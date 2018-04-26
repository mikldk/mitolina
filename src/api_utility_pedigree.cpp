#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "mitolina_types.h"

//' Get pedigree id
//' 
//' @export
// [[Rcpp::export]]
int get_pedigree_id(Rcpp::XPtr<Pedigree> ped) { 
  return ped->get_id();
}

bool is_pedigree_list(Rcpp::List pedigrees) {
  if (Rf_inherits(pedigrees, "mitolina_pedigreelist")) {
    return true;
  }
  
  return false;
}

void stopifnot_mitolina_pedigreelist(Rcpp::List pedigrees) {
  if (!is_pedigree_list(pedigrees)) {
    Rcpp::stop("pedigrees is not a mitolina_pedigreelist");
  }
}

//' Get number of pedigrees
//' 
//' @export
// [[Rcpp::export]]
int pedigrees_count(Rcpp::List pedigrees) {
  stopifnot_mitolina_pedigreelist(pedigrees);
  return pedigrees.size();
}

//' Get pedigree size
//' 
//' @export
// [[Rcpp::export]]
int pedigree_size(Rcpp::XPtr<Pedigree> ped) {  
  return ped->get_all_individuals()->size();
}

//' Get distribution of pedigree sizes
//' 
//' @export
//[[Rcpp::export]]
std::unordered_map<int, int> pedigrees_table(Rcpp::List pedigrees) {
  stopifnot_mitolina_pedigreelist(pedigrees);
  
  std::unordered_map<int, int> tab;
  
  for (int i = 0; i < pedigrees.size(); ++i) {
    Rcpp::XPtr<Pedigree> ped = pedigrees[i];
    tab[ped->get_all_individuals()->size()] += 1;
  }
  
  return tab;
}

/*
//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree_by_0index(Rcpp::List pedigrees, int index) {  
  stopifnot_mitolina_pedigreelist(pedigrees);
  
  if (index < 0 || index >= pedigrees.size()) {
    Rcpp::stop("pedigree index outside range");
  }
  
  Rcpp::XPtr<Pedigree> p = pedigrees.at(index);

  return p;
}
*/

/*
//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree_by_pedigree_id(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int index) {  
  std::vector<Pedigree*>* peds = pedigrees;
  Pedigree* p = peds->at(index);
  
  //Rcpp::XPtr<Pedigree> res(p, true);
  Rcpp::XPtr<Pedigree> res(p, false); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = Rcpp::CharacterVector::create("mitolina_pedigree", "externalptr");
  
  return res;
}
*/

//[[Rcpp::export]]
void print_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  Rcpp::Rcout << "Pedigree with " << p->get_all_individuals()->size() << " individuals:" << std::endl;
  
  for (auto i : *inds) {
    int pid_f = (i->get_mother() != NULL) ? i->get_mother()->get_pid() : -1;
    char gender = i->is_female() ? 'F' : 'M';
    
    Rcpp::Rcout << "  " << i->get_pid() << " [" << gender << "] with mother " << pid_f << std::endl;
  } 
}

//' get pids in pedigree
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  Rcpp::IntegerVector res(inds->size());
  int i = 0;
  for (auto ind : *inds) {   
    res(i) = ind->get_pid();
    ++i;
  } 
  
  return res;
}

//' get genders in pedigree
//' 
//' @export
// [[Rcpp::export]]
Rcpp::LogicalVector get_is_female_in_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  Rcpp::LogicalVector res(inds->size());
  int i = 0;
  for (auto ind : *inds) {   
    res(i) = ind->is_female();
    ++i;
  } 
  
  return res;
}

//' get pids in pedigree
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotypes_in_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
 
  size_t N = inds->size();
  Rcpp::List haps(N); 
  
  for (size_t i = 0; i < N; ++i) {
    Individual* indv = inds->at(i);
    haps(i) = indv->get_haplotype();
  }
  
  return haps;
}

//[[Rcpp::export]]
Rcpp::CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  Rcpp::CharacterMatrix edgelist(rels->size(), 2);
  int i = 0;
  
  for (auto pair: *rels) {
    edgelist(i, 0) = std::to_string(pair->first->get_pid());
    edgelist(i, 1) = std::to_string(pair->second->get_pid());
    ++i;
  }
  
  return edgelist;
}


//' Get pedigree information as graph (mainly intended for plotting)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  Rcpp::CharacterVector nodes(inds->size());
  
  int i = 0;
  for (auto individual : *inds) {
    nodes(i) = std::to_string(individual->get_pid());   
    ++i;
  }
  
  Rcpp::List ret;
  ret["nodes"] = nodes;
  ret["edgelist"] = get_pedigree_edgelist(ped);
  
  return ret;
}




















//' get pedigrees information in tidy format
//' 
// [[Rcpp::export]]
Rcpp::List get_pedigrees_tidy(Rcpp::List pedigrees) {  
  stopifnot_mitolina_pedigreelist(pedigrees);
  
  Rcpp::List ret_ped_ids;
  Rcpp::List ret_edgelists;
  Rcpp::List ret_haplotypes;
  Rcpp::List ret_pids;
  Rcpp::List ret_is_female;
  Rcpp::List ret_generation;  
  
  for (int ped_index = 0; ped_index < pedigrees.size(); ++ped_index) {
    Rcpp::XPtr<Pedigree> ped = pedigrees.at(ped_index);

    ret_ped_ids.push_back(ped->get_id());    
    
    std::vector< std::pair<Individual*, Individual*>* >* rels = ped->get_relations();
    
    Rcpp::IntegerMatrix edgelist(rels->size(), 2);
    int i = 0;
    
    for (auto pair: *rels) {
      edgelist(i, 0) = pair->first->get_pid();
      edgelist(i, 1) = pair->second->get_pid();
      ++i;
    }
    
    ret_edgelists.push_back(edgelist);

    
    std::vector<Individual*>* inds = ped->get_all_individuals();
    
    size_t N = inds->size();
    Rcpp::List haps(N);
    Rcpp::IntegerVector pids(N);
    Rcpp::LogicalVector is_female(N);
    Rcpp::IntegerVector generation(N);
    
    for (size_t i = 0; i < N; ++i) {
      Individual* indv = inds->at(i);
      haps(i) = indv->get_haplotype();
      pids(i) = indv->get_pid();
      is_female(i) = indv->is_female();
      generation(i) = indv->get_generations_from_final();
    }
    
    ret_haplotypes.push_back(haps);
    ret_pids.push_back(pids);
    ret_is_female.push_back(is_female);
    ret_generation.push_back(generation);
  }
  
  Rcpp::List ret;
  
  ret["ped_ids"] = ret_ped_ids;
  ret["edgelists"] = ret_edgelists;
  ret["haplotypes"] = ret_haplotypes;
  ret["pids"] = ret_pids;
  ret["is_female"] = ret_is_female;
  ret["generations_from_final"] = ret_generation;
  
  return ret;
}


