#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>

#include "mitolina_types.hpp"
#include "api_simulate_helper.hpp"

using namespace Rcpp;


//' Simulate a geneology with varying population size.
//' 
//' This function simulates a geneology with varying population size specified
//' by a vector of population sizes, one for each generation. 
//' 
//' By the backwards simulating process of the Wright-Fisher model, 
//' individuals with no descendants in the end population are not simulated 
//' If for some reason additional full generations should be simulated, 
//' the number can be specified via the \code{extra_generations_full} parameter.
//' This can for example be useful if one wants to simulate the 
//' final 3 generations although some of these may not get (male) children.
//' 
//' Let \eqn{\alpha} be the parameter of a symmetric Dirichlet distribution 
//' specifying each man's probability to be the mother of an arbitrary 
//' male in the next generation. When \eqn{\alpha = 5}, a man's relative probability 
//' to be the mother has 95\% probability to lie between 0.32 and 2.05, compared with a 
//' constant 1 under the standard Wright-Fisher model and the standard deviation in 
//' the number of male offspring per man is 1.10 (standard Wright-Fisher = 1).
//' 
//' This symmetric Dirichlet distribution is implemented by drawing 
//' mother (unscaled) probabilities from a Gamma distribution with 
//' parameters \code{gamma_parameter_shape} and \code{gamma_parameter_scale} 
//' that are then normalised to sum to 1. 
//' To obtain a symmetric Dirichlet distribution with parameter \eqn{\alpha}, 
//' the following must be used:
//' \eqn{\code{gamma_parameter_shape} = \alpha}
//' and 
//' \eqn{\code{gamma_parameter_scale} = 1/\alpha}.
//' 
//' @param population_sizes The size of the population at each generation, g. 
//'        population_sizes[g] is the population size at generation g.
//'        The length of population_sizes is the number of generations being simulated.
//' @param extra_generations_full Additional full generations to be simulated.
//' @param gamma_parameter_shape Parameter related to symmetric Dirichlet distribution for each man's probability to be mother. Refer to details.
//' @param gamma_parameter_scale Parameter realted to symmetric Dirichlet distribution for each man's probability to be mother. Refer to details.
//' @param enable_gamma_variance_extension Enable symmetric Dirichlet (and disable standard Wright-Fisher).
//' @param progress Show progress.
//' @param individuals_generations_return How many generations back to return (pointers to) individuals for.
//' 
//' @return A mitolina_simulation / list with the following entries:
//' \itemize{
//'   \item \code{population}. An external pointer to the population.
//'   \item \code{generations}. Generations actually simulated, mostly useful when parameter \code{generations = -1}.
//'   \item \code{founders}. Number of founders after the simulated \code{generations}.
//'   \item \code{growth_type}. Growth type model.
//'   \item \code{sdo_type}. Standard deviation in a man's number of male offspring. StandardWF or GammaVariation depending on \code{enable_gamma_variance_extension}.
//'   \item \code{end_generation_individuals}. Pointers to individuals in end generation.
//'   \item \code{individuals_generations}. Pointers to individuals in end generation in addition to the previous \code{individuals_generations_return}.
//' }
//' 
//' @import Rcpp
//' @import RcppProgress
//' @import RcppArmadillo
//' @export
// [[Rcpp::export]]
List sample_mtdna_geneology_varying_size(
  IntegerVector population_sizes,
  int extra_generations_full = 0,  
  double gamma_parameter_shape = 7, double gamma_parameter_scale = 7, 
  bool enable_gamma_variance_extension = false,
  bool progress = true, 
  int individuals_generations_return = 2) {
  
  // boolean chosen like this to obey NA's
  bool all_gt_1 = is_true(all(population_sizes >= 1));
  if (!all_gt_1) {
    Rcpp::stop("Please specify only population_sizes >= 1");
  }
  
  int generations = population_sizes.length();
  
  if (generations < -1 || generations == 0) {
    Rcpp::stop("Please specify generations as -1 (for simulation to 1 founder) or > 0");
  }

  if (enable_gamma_variance_extension) {
    if (gamma_parameter_shape <= 0.0) {
      Rcpp::stop("gamma_parameter_shape must be > 0.0");
    }
    if (gamma_parameter_scale <= 0.0) {
      Rcpp::stop("gamma_parameter_scale must be > 0.0");
    }
  }
  
  //Rcpp::Rcout << "vary 1" << std::endl;

  
  Progress progress_bar(generations, progress);
  
  std::unordered_map<int, Individual*>* population_map = new std::unordered_map<int, Individual*>(); // pid's are garanteed to be unique
  Population* population = new Population(population_map);
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG);
  population_xptr.attr("class") = CharacterVector::create("mitolina_population", "externalptr");
  
  
  int individual_id = 1;
  std::vector<Individual*> end_generation(population_sizes[generations-1]);
  List end_generation_individuals(population_sizes[generations-1]);
  List last_k_generations_individuals;

  for (size_t i = 0; i < population_sizes[generations-1]; ++i) {
    Individual* indv = new Individual(individual_id++, 0);
    end_generation[i] = indv;    
    (*population_map)[indv->get_pid()] = indv;
    
    Rcpp::XPtr<Individual> indv_xptr(indv, RCPP_XPTR_2ND_ARG);
    end_generation_individuals[i] = indv_xptr;
    
    if (individuals_generations_return >= 0) {
      last_k_generations_individuals.push_back(indv_xptr);
    }
  }
  
  if (progress) {
    progress_bar.increment();
  }
  
  //Rcpp::Rcout << "vary 2" << std::endl;
  //Rcpp::Rcout << "  population_sizes[generations-1] = population_sizes[" << (generations-1) << "]" << std::endl;
  
  // Next generation  
  //std::vector<Individual*>* children_generation = &end_generation;
  std::vector<Individual*> children_generation(population_sizes[generations-1]);
  for (size_t i = 0; i < population_sizes[generations-1]; ++i) {
    children_generation[i] = end_generation[i];
  }
  std::vector<Individual*> mothers_generation;
  
  int founders_left = population_sizes[generations-1];
  
  //Rcpp::Rcout << "vary 3" << std::endl;
  //Rcpp::Rcout << "  generations = " << generations << std::endl;
  //Rcpp::Rcout << "  population_sizes.length() = " << population_sizes.length() << std::endl;
  
  // now, find out who the mothers to the children are
  for (size_t generation = 1; generation < generations; ++generation) {
    // Init ->
    //Rcpp::Rcout << "vary 4-" << generation << std::endl;

    //Rcpp::Rcout << "  population_size          = population_sizes[generations-(generation+1)] = population_sizes[" << generations << "-" << (generation+1) << "] = population_sizes[" << (generations-(generation+1)) << "]" << std::endl;
    //Rcpp::Rcout << "  children_population_size = population_sizes[generations-generation] = population_sizes[" << generations << "-" << generation << "] = population_sizes[" << (generations-generation) << "]" << std::endl;

    int population_size = population_sizes[generations-(generation+1)];    
    int children_population_size = population_sizes[generations-generation];

    
    WFRandomMother wf_random_mother(population_size);
    GammaVarianceRandomMother gamma_variance_mother(population_size, gamma_parameter_shape, gamma_parameter_scale);  
    SimulateChooseMother* choose_mother = &wf_random_mother;
    if (enable_gamma_variance_extension) {
      choose_mother = &gamma_variance_mother;
    }
    
    //Rcpp::Rcout << "vary 5-" << generation << std::endl;
    
    mothers_generation.clear();
    mothers_generation.resize(population_size);
    
    // <- Init
    
    int new_founders_left = 0;
    //Rcpp::Rcerr << "Generation " << generation << std::endl;
    
    // clear
    for (size_t i = 0; i < population_size; ++i) { // necessary?
      mothers_generation[i] = nullptr;
    }
    
    //Rcpp::Rcout << "vary 6-" << generation << std::endl;
    
    choose_mother->update_state_new_generation();
    
    //Rcpp::Rcout << "vary 7-" << generation << std::endl;
    
    /*
    FIXME: New, more explicit logic? Maybe copy create_mother_update_simulation_state_varying_size into here?
    */
    
    // now, run through children to pick each child's mother
    for (size_t i = 0; i < children_population_size; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (children_generation[i] == nullptr) {
        continue;
      }
      
      // child [i] in [generation-1]/children_generation has mother [mother_i] in [generation]/mothers_generation
      //int mother_i = sample_person_weighted(population_size, mothers_prob, mothers_prob_perm);
      int mother_i = choose_mother->get_mother_i();
      
      // if this is the mother's first child, create the mother
      if (mothers_generation[mother_i] == nullptr) {
        create_mother_update_simulation_state_varying_size(mother_i, &individual_id, generation, 
              individuals_generations_return, mothers_generation, population_map, 
              &new_founders_left, last_k_generations_individuals);
      }
      
      children_generation[i]->set_mother(mothers_generation[mother_i]);
      mothers_generation[mother_i]->add_child(children_generation[i]);
    }
    
    //Rcpp::Rcout << "vary 8-" << generation << std::endl;
    
    //Rcpp::Rcout << "generation = " << generation << "; extra_generations_full = " << extra_generations_full << std::endl;
    
    // create additional mothers (without children) if needed:
    if (generation <= extra_generations_full) {
      for (size_t mother_i = 0; mother_i < population_size; ++mother_i) {
        //Rcpp::Rcout << "vary 8-dong-" << mother_i << std::endl;
        
        if (mothers_generation[mother_i] != nullptr) {
          continue;
        }        
        
        // create mother, no children etc.
        create_mother_update_simulation_state_varying_size(mother_i, &individual_id, generation, 
              individuals_generations_return, mothers_generation, population_map, 
              &new_founders_left, last_k_generations_individuals);
      }      
    }
    
    //Rcpp::Rcout << "vary 9-" << generation << std::endl;
    
    // children_generation = &mothers_generation;
    // FIXME mikl 2017-06-26 19:09
    //children_generation = mothers_generation;
    //vs:
    children_generation.clear();
    children_generation.resize(population_size);
    for (size_t i = 0; i < population_size; ++i) {
      children_generation[i] = mothers_generation[i];
    }
    //<-
    
    if (Progress::check_abort()) {
      stop("Aborted");
    }
    
    if (progress) {
      progress_bar.increment();
    }
    
    founders_left = new_founders_left;
  }
  
  List res;
  res["population"] = population_xptr;
  res["generations"] = generations;
  res["founders"] = founders_left;
  res["growth_type"] = "VaryingPopulationSize";
  res["sdo_type"] = (enable_gamma_variance_extension) ? "GammaVariation" : "StandardWF";
  res["end_generation_individuals"] = end_generation_individuals;
  res["individuals_generations"] = last_k_generations_individuals;

  res.attr("class") = CharacterVector::create("mitolina_simulation", "list");
  
  return res;
}

