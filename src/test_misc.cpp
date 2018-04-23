/**
 test_misc.cpp
 Purpose: Code used in testing package.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */
 
#include <RcppArmadillo.h>

#include "mitolina_types.h"

using namespace Rcpp;

//' Generate test population
//' 
//' @return An external pointer to the population.
// [[Rcpp::export]]
Rcpp::XPtr<Population> test_create_population() {
  
  // pid's are garanteed to be unique
  std::unordered_map<int, Individual*>* population_map = 
    new std::unordered_map<int, Individual*>(); 
  
  Population* population = new Population(population_map);
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG);
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  //////////

  std::vector<Individual*> indvs;

  Individual* i1 = new Individual(1, 0, false); indvs.push_back(i1);
  Individual* i2 = new Individual(2, 0, true); indvs.push_back(i2);
  Individual* i3 = new Individual(3, 0, false); indvs.push_back(i3);
  Individual* i4 = new Individual(4, 0, false); indvs.push_back(i4);
  Individual* i5 = new Individual(5, 0, false); indvs.push_back(i5);
  
  Individual* i6 = new Individual(6, 1, true); indvs.push_back(i6);
  Individual* i7 = new Individual(7, 1, true); indvs.push_back(i7);
  Individual* i8 = new Individual(8, 1, true); indvs.push_back(i8);
  
  Individual* i9 = new Individual(9, 2, true); indvs.push_back(i9);
  Individual* i10 = new Individual(10, 2, true); indvs.push_back(i10);
  
  Individual* i11 = new Individual(11, 3, true); indvs.push_back(i11);
  
  Individual* i12 = new Individual(12, 3, true); indvs.push_back(i12);
  Individual* i13 = new Individual(13, 2, true); indvs.push_back(i13);
  Individual* i14 = new Individual(14, 1, true); indvs.push_back(i14);
  Individual* i15 = new Individual(15, 3, false); indvs.push_back(i15);
  Individual* i16 = new Individual(16, 3, false); indvs.push_back(i16);
  Individual* i17 = new Individual(17, 3, false); indvs.push_back(i17);
  
  Individual* i18 = new Individual(18, 2, false); indvs.push_back(i18);  
  Individual* i19 = new Individual(19, 2, false); indvs.push_back(i19);
  
  /*
   *     
   * G3          F11               F12
   *           /     \              |
   * G2     F9        F10          F13
   *        /  \       | \          |   \
   * G1   F6    7     F8  M19      F14   18M
   *       |   /\     /\         /  |  \
   * G0   M1  F2 M3  M4 M5     M15 M16 M17
   */
  i11->add_child(i9);
  i11->add_child(i10);
  i9->add_child(i6);
  i9->add_child(i7);
  i10->add_child(i8);
  i10->add_child(i19);
  
  i6->add_child(i1);
  i7->add_child(i2);
  i7->add_child(i3);
  
  i8->add_child(i4);
  i8->add_child(i5);
  
  //
  
  i12->add_child(i13);
  i13->add_child(i14);
  i13->add_child(i18);
  
  i14->add_child(i15);
  i14->add_child(i16);
  i14->add_child(i17);
  
  for (auto i : indvs) {
    (*population_map)[i->get_pid()] = i;
  }
  
  return population_xptr;
}

