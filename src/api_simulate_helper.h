#ifndef MITOLINE_API_SIMULATE_HELPER_H
#define MITOLINE_API_SIMULATE_HELPER_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>

#include "mitolina_types.h"

using namespace Rcpp;

void create_mother_update_simulation_state_varying_size(
  int mother_i, 
  int* individual_id, 
  int generation, 
  int individuals_generations_return,
  std::vector<Individual*>& mothers_generation, 
  std::unordered_map<int, Individual*>* population_map, 
  int* new_founders_left,
  List& last_k_generations_individuals);
  
#endif
