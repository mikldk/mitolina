#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>

#include "mitolina_types.hpp"
#include "api_simulate.hpp"

using namespace Rcpp;

void create_mother_update_simulation_state_varying_size(
  int mother_i, 
  int* individual_id, 
  int generation, 
  int individuals_generations_return,
  std::vector<Individual*>& mothers_generation, 
  std::unordered_map<int, Individual*>* population_map, 
  int* new_founders_left,
  List& last_k_generations_individuals) {  
  
  Individual* mother = new Individual(*individual_id, generation);
  (*individual_id) = (*individual_id) + 1;
  
  mothers_generation[mother_i] = mother;
  (*population_map)[mother->get_pid()] = mother;

  (*new_founders_left) = (*new_founders_left) + 1;

  if (generation <= individuals_generations_return) {
    //Rcpp::Rcout << "create_mother_update_simulation_state_varying_size: generation = " << generation << "; individuals_generations_return = " << individuals_generations_return << std::endl;
    
    Rcpp::XPtr<Individual> mother_xptr(mother, RCPP_XPTR_2ND_ARG);
    last_k_generations_individuals.push_back(mother_xptr);
  }  
}

