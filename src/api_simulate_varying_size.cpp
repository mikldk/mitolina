#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>

#include "mitolina_types.h"
#include "api_simulate_helper.h"

using namespace Rcpp;


//' Simulate a geneology with varying population size.
//' 
//' This function simulates a geneology with varying population size specified
//' by a vector of population sizes, one for each generation. 
//' 
//' By the backwards simulating process of the Wright-Fisher model, 
//' individuals with no descendants in the end population are not simulated 
//' If for some reason additional full generations should be simulated, 
//' the number can be specified via the `generations_full` parameter.
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
//' parameters `gamma_parameter_shape` and `gamma_parameter_scale`
//' that are then normalised to sum to 1. 
//' To obtain a symmetric Dirichlet distribution with parameter \eqn{\alpha}, 
//' the following must be used:
//' \eqn{`gamma_parameter_shape` = \alpha}
//' and 
//' \eqn{`gamma_parameter_scale` = 1/\alpha}.
//' 
//' @param population_sizes_females The size of the female population at each generation, `g`. All >= 1.
//'        `population_sizes_females[g]` is the population size at generation `g`.
//'        The length of population_sizes_females is the number of generations being simulated.
//' @param population_sizes_males The size of the male population at each generation, `g`. All >= 0.
//'        `population_sizes_males[g]` is the population size at generation `g`.
//' @param generations_full Number of full generations to be simulated.
//' @param generations_return How many generations to return (pointers to) individuals for.
//' @param enable_gamma_variance_extension Enable symmetric Dirichlet (and disable standard Wright-Fisher).
//' @param gamma_parameter_shape Parameter related to symmetric Dirichlet distribution for each man's probability to be mother. Refer to details.
//' @param gamma_parameter_scale Parameter realted to symmetric Dirichlet distribution for each man's probability to be mother. Refer to details.
//' @param progress Show progress.
//' 
//' @return A `mitolina_simulation` / list with the following entries:
//' \itemize{
//'   \item `population`. An external pointer to the population.
//'   \item `generations`. Generations actually simulated, mostly useful when parameter `generations = -1`.
//'   \item `founders`. Number of founders after the simulated `generations`.
//'   \item `growth_type`. Growth type model.
//'   \item `sdo_type`. Standard deviation in a man's number of male offspring. StandardWF or GammaVariation depending on `enable_gamma_variance_extension`.
//'   \item `end_generation_female_individuals`. Pointers to female individuals in end generation.
//'   \item `female_individuals_generations`. Pointers to female individuals in last `generations_return` generation (if `generations_return = 3`, then female individuals in the last three generations are returned).
//'   \item `end_generation_male_individuals`. Pointers to male individuals in end generation.
//'   \item `male_individuals_generations`. Pointers to male individuals in last `generations_return` generation (if `generations_return = 3`, then male individuals in the last three generations are returned).
//' }
//' @import Rcpp
//' @import RcppProgress
//' @import RcppArmadillo
//' @export
// [[Rcpp::export]]
List sample_mtdna_geneology_varying_size(
  IntegerVector population_sizes_females,
  IntegerVector population_sizes_males,
  int generations_full = 1,  
  int generations_return = 3,
  bool enable_gamma_variance_extension = false,
  double gamma_parameter_shape = 5.0, double gamma_parameter_scale = 1.0/5.0, 
  bool progress = true) {
  
  if (generations_full <= 0) {
    Rcpp::stop("generations_full must be at least 1");
  }
  // Always include full last generation, but how many additional?
  int extra_generations_full = generations_full - 1;
  
  if (generations_return <= 0) {
    Rcpp::stop("generations_return must be at least 1");
  }  
  // Always include full last generation, but how many additional?
  int extra_individuals_generations_return = generations_return - 1;
  
  
  // boolean chosen like this to obey NA's
  bool all_gte_1_females = is_true(all(population_sizes_females >= 1));
  if (!all_gte_1_females) {
    Rcpp::stop("Please specify only population_sizes_females >= 1");
  }
  
  // boolean chosen like this to obey NA's
  bool all_gte_0_males = is_true(all(population_sizes_males >= 0));
  if (!all_gte_0_males) {
    Rcpp::stop("Please specify only population_sizes_males >= 0");
  }
  
  int generations = population_sizes_females.length();
  
  if (generations == 0) {
    Rcpp::stop("Please specify at least 1 generation (the vector population_sizes must have length >= 1)");
  }
  
  if (population_sizes_males.length() != generations) {
    Rcpp::stop("length(population_sizes_males) != length(population_sizes_females)");
  }
  
  if (population_sizes_males[0] != 0) {
    Rcpp::stop("Please specify 0 in the first male generation (no male founders)");
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
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG_CLEANER);
  population_xptr.attr("class") = CharacterVector::create("mitolina_population", "externalptr");
  
  
  
  
  
  int individual_id = 1;
  
  
  
  // Create females
  std::vector<Individual*> end_generation_females(population_sizes_females[generations - 1]);
  List end_generation_female_individuals(population_sizes_females[generations - 1]);
  List last_k_generations_female_individuals;

  for (size_t i = 0; i < population_sizes_females[generations - 1]; ++i) {
    Individual* indv = new Individual(individual_id++, 0, true); // is_female = true: hence a female
    end_generation_females[i] = indv;    
    (*population_map)[indv->get_pid()] = indv;
    
    Rcpp::XPtr<Individual> indv_xptr(indv, RCPP_XPTR_2ND_ARG);
    end_generation_female_individuals[i] = indv_xptr;
    
    if (extra_individuals_generations_return >= 0) {
      last_k_generations_female_individuals.push_back(indv_xptr);
    }
  }
  
  // Create males
  std::vector<Individual*> end_generation_males(population_sizes_males[generations - 1]);
  List end_generation_male_individuals(population_sizes_males[generations - 1]);
  List last_k_generations_male_individuals;

  for (size_t i = 0; i < population_sizes_males[generations - 1]; ++i) {
    Individual* indv = new Individual(individual_id++, 0, false); // is_female = false: hence a male
    end_generation_males[i] = indv;    
    (*population_map)[indv->get_pid()] = indv;
    
    Rcpp::XPtr<Individual> indv_xptr(indv, RCPP_XPTR_2ND_ARG);
    end_generation_male_individuals[i] = indv_xptr;
    
    if (extra_individuals_generations_return >= 0) {
      last_k_generations_male_individuals.push_back(indv_xptr);
    }
  }
  
  
  if (progress) {
    progress_bar.increment();
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Next generation
  /////////////////////////////////////////////////////////////////////////////////////////////
  
  // Female init
  std::vector<Individual*> female_children_generation(population_sizes_females[generations-1]);
  for (size_t i = 0; i < population_sizes_females[generations-1]; ++i) {
    female_children_generation[i] = end_generation_females[i];
  }
  
  // Male init    
  std::vector<Individual*> male_children_generation(population_sizes_males[generations-1]);
  for (size_t i = 0; i < population_sizes_males[generations-1]; ++i) {
    male_children_generation[i] = end_generation_males[i];
  }  
  
  // Generations
  std::vector<Individual*> mothers_generation;  
  int founders_left = population_sizes_females[generations - 1];
  
  // now, find out who the mothers to the children (both females and males!) are
  for (size_t generation = 1; generation < generations; ++generation) {
    // Init ->
  
    // Mother choosing mechanism
    int mother_population_size = population_sizes_females[generations-(generation+1)];
    
    WFRandomMother wf_random_mother(mother_population_size);
    GammaVarianceRandomMother gamma_variance_mother(mother_population_size, gamma_parameter_shape, gamma_parameter_scale);  
    SimulateChooseMother* choose_mother = &wf_random_mother;
    if (enable_gamma_variance_extension) {
      choose_mother = &gamma_variance_mother;
    }
    // FIXME: Why not simply put this in constructor?
    choose_mother->update_state_new_generation();
    
    //--------------------------------------------------------------------------
    
    // Variables to keep tract of selected mothers
    
    mothers_generation.clear();
    mothers_generation.resize(mother_population_size);    
    // clear
    for (size_t i = 0; i < mother_population_size; ++i) { // necessary?
      mothers_generation[i] = nullptr;
    }
    int new_founders_left = 0;
    
    
    //--------------------------------------------------------------------------
    // now, run through **FEMALE** children to pick each child's mother
    int female_children_population_size = population_sizes_females[generations-generation];
    for (size_t i = 0; i < female_children_population_size; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (female_children_generation[i] == nullptr) {
        continue;
      }
      
      int mother_i = choose_mother->get_mother_i();
      
      // if this is the mother's first child, create the mother
      if (mothers_generation[mother_i] == nullptr) {
        create_mother_update_simulation_state_varying_size(mother_i, &individual_id, generation, 
              extra_individuals_generations_return, mothers_generation, population_map, 
              &new_founders_left, last_k_generations_female_individuals);
      }
      
      // mother set in add_child
      mothers_generation[mother_i]->add_child(female_children_generation[i]);
    }
    //--------------------------------------------------------------------------
    // now, run through **MALE** children to pick each child's mother
    int male_children_population_size = population_sizes_males[generations-generation];
    for (size_t i = 0; i < male_children_population_size; ++i) {
      int mother_i = choose_mother->get_mother_i();
      
      // if this is the mother's first child, create the mother
      if (mothers_generation[mother_i] == nullptr) {
        create_mother_update_simulation_state_varying_size(mother_i, &individual_id, generation, 
              extra_individuals_generations_return, mothers_generation, population_map, 
              &new_founders_left, last_k_generations_female_individuals);
      }

      // mother set in add_child
      mothers_generation[mother_i]->add_child(male_children_generation[i]);
    }    
    
    // Creates previous generation's males (i.e. for next iteration)
    male_children_population_size = population_sizes_males[generations-(generation+1)];
    male_children_generation.clear();
    male_children_generation.resize(male_children_population_size);
        
    for (size_t i = 0; i < male_children_population_size; ++i) {
      Individual* indv = new Individual(individual_id++, generation, false); // is_female = false: hence a male
      (*population_map)[indv->get_pid()] = indv;
      male_children_generation[i] = indv;
    }
    
    // possibly also add these new to last_k_generations_male_individuals
    if (generation <= extra_individuals_generations_return) {
      for (size_t i = 0; i < male_children_population_size; ++i) {
        Rcpp::XPtr<Individual> male_xptr(male_children_generation[i], RCPP_XPTR_2ND_ARG);
        last_k_generations_male_individuals.push_back(male_xptr);
      }
    }
    

    //--------------------------------------------------------------------------
    
    
    
    //Rcpp::Rcout << "vary 8-" << generation << std::endl;
    
    //Rcpp::Rcout << "generation = " << generation << "; extra_generations_full = " << extra_generations_full << std::endl;
    
    // create additional mothers (without children) if needed:
    if (generation <= extra_generations_full) {
      for (size_t mother_i = 0; mother_i < mother_population_size; ++mother_i) {
        //Rcpp::Rcout << "vary 8-dong-" << mother_i << std::endl;
        
        if (mothers_generation[mother_i] != nullptr) {
          continue;
        }        

        //Rcpp::Rcout << "   CREATE" << std::endl;
        
        //Rcpp::Rcout << "     " << individual_id << "->";
        
        // create mother, no children etc.
        create_mother_update_simulation_state_varying_size(mother_i, &individual_id, generation, 
              extra_individuals_generations_return, mothers_generation, population_map, 
              &new_founders_left, last_k_generations_female_individuals);
              
        // STOP: The created mother must have a mother next iteration!
              
        //Rcpp::Rcout << individual_id << std::endl;
      }      
    }
    
    //Rcpp::Rcout << "vary 9-" << generation << std::endl;
    
    female_children_generation.clear();
    female_children_generation.resize(mother_population_size);
    for (size_t i = 0; i < mother_population_size; ++i) {
      female_children_generation[i] = mothers_generation[i];
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
  res["end_generation_female_individuals"] = end_generation_female_individuals;
  res["female_individuals_generations"] = last_k_generations_female_individuals;
  res["end_generation_male_individuals"] = end_generation_male_individuals;
  res["male_individuals_generations"] = last_k_generations_male_individuals;

  res.attr("class") = CharacterVector::create("mitolina_simulation", "list");
  
  return res;
}

