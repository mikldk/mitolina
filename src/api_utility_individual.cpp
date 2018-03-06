#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "mitolina_types.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid) {  
  Population* pop = population;
  
  Individual* ind = population->get_individual(pid);
  //Rcpp::XPtr<Individual> res(ind, true);
  Rcpp::XPtr<Individual> res(ind, false); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
  res.attr("class") = Rcpp::CharacterVector::create("mitolina_individual", "externalptr");
  
  return res;
}


//' @export
// [[Rcpp::export]]
int get_pid(Rcpp::XPtr<Individual> individual) {  
  return individual->get_pid();
}

//' @export
// [[Rcpp::export]]
void print_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  int pid_f = (i->get_mother() != nullptr) ? i->get_mother()->get_pid() : -1;
  char gender = i->is_female() ? 'F' : 'M';
    
  std::vector<Individual*>* children = i->get_children();
  
  Rcpp::Rcout << "  pid = " << i->get_pid() << " [" << gender << "] with mother pid = " << pid_f << " and";
  
  if (children->size() == 0) {
    Rcpp::Rcout << " no children" << std::endl;
  } else {
    Rcpp::Rcout << " children (n = " << children->size() << "): " << std::endl;

    for (auto child : *children) {    
      std::vector<Individual*>* child_children = child->get_children();
      
      char child_gender = child->is_female() ? 'F' : 'M';
      
      Rcpp::Rcout << "    pid = " << child->get_pid() << " [" << gender << "] with mother pid = " << pid_f << " and " <<  child_children->size() << " children" << std::endl;
    }
  }
}

//' Get individual's generations
//' 
//' @export
// [[Rcpp::export]]
int get_generations_from_final(Rcpp::XPtr<Individual> individual) {  
  return individual->get_generations_from_final();
}

//' Get pedigree from individual
//' 
//' @export
//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;  
  Rcpp::XPtr<Pedigree> res(i->get_pedigree(), false); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = Rcpp::CharacterVector::create("mitolina_pedigree", "externalptr");
  
  return res;
}

//' Get pedigree id from pid
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids) {  
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  
  int N = pids.size();
  Rcpp::IntegerVector pedigree_ids(N);
  
  for (int i = 0; i < N; ++i) {
    Individual* ind = population->get_individual(pids[i]);
    pedigree_ids[i] = ind->get_pedigree_id();
  }
  
  return pedigree_ids;
}



//////////////////////////////////////

/*
//' @export
// [[Rcpp::export]]
int count_brothers(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_mother() == nullptr) {
    Rcpp::stop("Individual did not have a mother");
  }
  
  int mothers_girls = i->get_mother()->get_children_count();
  // exclude individual
  return (mothers_girls - 1);
}
*/
/*
//' @export
// [[Rcpp::export]]
int brothers_matching(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_mother() == nullptr) {
    Rcpp::stop("Individual did not have a mother");
  }
  
  if (!(i->is_haplotype_set())) {
    Rcpp::stop("Individual did not have a haplotype");
  }

  std::vector<bool> h = i->get_haplotype();  
  int loci = h.size();  
  
  std::vector<Individual*>* brothers = i->get_mother()->get_children();
  
  if (brothers->size() == 0) {
    return 0;
  }
  
  int matching = 0;

  for (auto brother : *brothers) {
    if (brother->get_pid() == i->get_pid()) {
      continue;
    }

    if (!(brother->is_haplotype_set())) {
      Rcpp::stop("Individual's brother did not have a haplotype");
    }  
  
    std::vector<bool> indv_h = brother->get_haplotype();    
    if (indv_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (indv_h == h) {
      matching += 1;
    }
  }
  
  return matching;
}

//' @export
// [[Rcpp::export]]
bool mother_matches(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;

  if (!(i->is_haplotype_set())) {
    Rcpp::stop("Individual did not have a haplotype");
  }
    
  if (i->get_mother() == nullptr) {
    Rcpp::stop("Individual did not have a mother");
  }
  
  Individual* mother = i->get_mother();
  
  if (!(mother->is_haplotype_set())) {
    Rcpp::stop("Individual's mother did not have a haplotype");
  }  
  
  std::vector<bool> h = i->get_haplotype();
  std::vector<bool> h_mother = mother->get_haplotype();    
  
  return (h.size() == h_mother.size() && h == h_mother);
}

//' @export
// [[Rcpp::export]]
bool grandmother_matches(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;

  if (!(i->is_haplotype_set())) {
    Rcpp::stop("Individual did not have a haplotype");
  }
  
  
  if (i->get_mother() == nullptr) {
    Rcpp::stop("Individual did not have a mother");
  }  
  Individual* mother = i->get_mother();
  if (!(mother->is_haplotype_set())) {
    Rcpp::stop("Individual's mother did not have a haplotype");
  }
  
  // It is not sufficient to calculate mother_matches(mother) as
  // mother and grandmother may not match

  if (mother->get_mother() == nullptr) {
    Rcpp::stop("Individual's mother did not have a mother");
  }
  Individual* grandmother = mother->get_mother();  
  if (!(grandmother->is_haplotype_set())) {
    Rcpp::stop("Individual's grandmother did not have a haplotype");
  }  
  
  std::vector<bool> h = i->get_haplotype();
  std::vector<bool> h_grandmother = grandmother->get_haplotype();    
  
  return (h.size() == h_grandmother.size() && h == h_grandmother);
}

*/


