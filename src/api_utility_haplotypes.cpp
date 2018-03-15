#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "mitolina_types.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List indices_in_mixture(Rcpp::IntegerMatrix haplotypes, Rcpp::IntegerVector H1, Rcpp::IntegerVector H2) { 
  size_t N = haplotypes.nrow();
  
  Rcpp::List res;
  
  if (N == 0) {
    return res;
  }

  // mainly count wanted, but indices are good for debuggin
  Rcpp::IntegerVector res_in_mixture;
  Rcpp::IntegerVector res_H1;
  Rcpp::IntegerVector res_H2;
  
  size_t loci = haplotypes.ncol();

  for (size_t i = 0; i < N; ++i) {
    Rcpp::IntegerVector h = haplotypes(i, Rcpp::_);
    
    bool in_mixture = true;
    bool match_H1 = true; // faster than Rcpp equal/all sugar 
    bool match_H2 = true;
    
    for (size_t locus = 0; locus < loci; ++locus) {
      if (in_mixture && (h[locus] != H1[locus]) && (h[locus] != H2[locus])) {
        in_mixture = false;
      }
      
      if (match_H1 && (h[locus] != H1[locus])) {
        match_H1 = false;
      }
      
      if (match_H2 && (h[locus] != H2[locus])) {
        match_H2 = false;
      }
      
      // if neither have a chance, just stop
      if (!in_mixture && !match_H1 && !match_H2) {
        break;
      }
    }
    
    if (in_mixture) {
      res_in_mixture.push_back(i + 1); // R indexing
    }
    
    if (match_H1) {
      res_H1.push_back(i + 1); // R indexing
    }
    
    if (match_H2) {
      res_H2.push_back(i + 1); // R indexing
    }
  }
  
  res["in_mixture"] = res_in_mixture;
  res["match_H1"] = res_H1;
  res["match_H2"] = res_H2;

  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::List pedigree_get_haplotypes_pids(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids) {  
  size_t N = pids.size();
  Rcpp::List haps(N);

  for (size_t i = 0; i < N; ++i) {
    Individual* indv = population->get_individual(pids[i]);
    haps(i) = indv->get_haplotype();
  }

  return haps;
}
 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix individuals_get_haplotypes(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals) {   
  size_t n = individuals.size();
 
  if (n <= 0) {
    Rcpp::IntegerMatrix empty_haps(0, 0);
    return empty_haps;
  }
 
  size_t loci = individuals[0]->get_haplotype().size();

  if (loci <= 0) {
    Rcpp::stop("Expected > 0 loci");
    Rcpp::IntegerMatrix empty_haps(0, 0);
    return empty_haps;
  }

  Rcpp::IntegerMatrix haps(n, loci);

  for (size_t i = 0; i < n; ++i) {
    std::vector<bool> hap = individuals[i]->get_haplotype();

    if (hap.size() != loci) {
      Rcpp::stop("Expected > 0 loci for all haplotypes");
      Rcpp::IntegerMatrix empty_haps(0, 0);
      return empty_haps;
    }
    
    Rcpp::IntegerVector h = Rcpp::wrap(hap);
    haps(i, Rcpp::_) = h;
  }

  return haps;
}

//' @export
// [[Rcpp::export]]
void pedigree_populate_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, Rcpp::NumericVector mutation_rates) {  
  Pedigree* p = ped;
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);

  if (loci != mut_rates.size()) {
    Rcpp::stop("Number of loci specified in haplotype must equal number of mutation rates specified");
  }
    
  ped->populate_haplotypes(loci, mut_rates);
}

//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, Rcpp::NumericVector mutation_rates, bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  
  if (loci != mut_rates.size()) {
    Rcpp::stop("Number of loci specified in haplotype must equal number of mutation rates specified");
  }
  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_haplotypes(loci, mut_rates);
    
     if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}

//' Custom founders
//' 
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes_custom_founders(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, 
                                       Rcpp::NumericVector mutation_rates,
                                       Rcpp::Nullable<Rcpp::Function> get_founder_haplotype = R_NilValue,
                                       bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  
  if (get_founder_haplotype.isNull()) {
    Rcpp::stop("get_founder_haplotype must not be NULL");
  }  
  
  Rcpp::Function g_founder_hap = Rcpp::as<Rcpp::Function>(get_founder_haplotype);

  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_haplotypes_custom_founders(mut_rates, g_founder_hap);
    
     if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}

//' @export
// [[Rcpp::export]]
std::vector<bool> get_haplotype(Rcpp::XPtr<Individual> individual) {
  return individual->get_haplotype();
}

//' @export
// [[Rcpp::export]]
int get_haplotype_no_variants(Rcpp::XPtr<Individual> individual) {
  return individual->get_haplotype_total_no_variants();
}

//' @export
// [[Rcpp::export]]
int count_haplotype_occurrences_individuals(const Rcpp::List individuals, 
    const Rcpp::LogicalVector haplotype, 
    const int haplotype_no_variants) {
    
  int n = individuals.size();
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotype);
  
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];

    if (haplotype_no_variants != indv->get_haplotype_total_no_variants()) {
      continue;
    }
    
    std::vector<bool> indv_h = indv->get_haplotype();
    
    if (indv_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (indv_h == h) {
      count += 1;
    }
  }
  
  return count;
}

//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotype_matching_individuals(const Rcpp::List individuals, 
    const Rcpp::LogicalVector haplotype, 
    const int haplotype_no_variants) {
    
  int n = individuals.size();
  int loci = haplotype.size();
  Rcpp::List ret_indv;
  
  std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotype);
  
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];

    if (haplotype_no_variants != indv->get_haplotype_total_no_variants()) {
      continue;
    }
    
    std::vector<bool> indv_h = indv->get_haplotype();
    
    if (indv_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (indv_h == h) {
      ret_indv.push_back(indv);
    }
  }
  
  return ret_indv;
}

// There are count_haplotype_occurrences_pedigree matches. 
// This gives details on meiotic distance and the max L0 distance of the haplotypes on the path between the suspect and the matching individual in the pedigree. Some of these matches may have (back)mutations between in between them.
//
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix pedigree_haplotype_matches_in_pedigree_meiosis_L0_dists(const Rcpp::XPtr<Individual> suspect, bool matches_are_female = true, int generation_upper_bound_in_result = -1) {
  Rcpp::stop("Extract matches, then count them and then find this information");
  
  const std::vector<bool> h = suspect->get_haplotype();
  const int h_no_variants = suspect->get_haplotype_total_no_variants();

  const Pedigree* pedigree = suspect->get_pedigree();
  const int suspect_pedigree_id = suspect->get_pedigree_id();
  const std::vector<Individual*>* family = pedigree->get_all_individuals();
  
  std::vector<int> meiosis_dists;
  std::vector<int> max_L0_dists;
  std::vector<int> pids;
  
  // includes suspect by purpose
  for (auto dest : *family) { 
    int generation = dest->get_generations_from_final();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    if (dest->is_female() != matches_are_female) {
      continue;
    }
    
    // only considering within pedigree matches
    if (dest->get_pedigree_id() != suspect_pedigree_id) {
      continue;
    }
    
    if (h_no_variants != dest->get_haplotype_total_no_variants()) {
      continue;
    }
    
    std::vector<bool> dest_h = dest->get_haplotype();
    
    if (dest_h.size() != h.size()) {
      Rcpp::stop("haplotype and dest_h did not have same number of loci");
    }
    
    if (dest_h == h) {
      std::vector<Individual*> path = suspect->calculate_path_to(dest);  
      int meiosis_dist = suspect->meiosis_dist_tree(dest);
      
      int meiosis_dist_from_path = path.size() - 1; // n vertices means n-1 edges (tree)
      //Rcpp::Rcout << ">> path from " << suspect->get_pid() << " to " << dest->get_pid() << " has length = " << meiosis_dist_from_path << " and meioses = " << meiosis_dist << (meiosis_dist_from_path == meiosis_dist ? " ok" : " ERROR") << ": " << std::endl;
      
      int max_L0 = 0;
      
      //Rcpp::Rcout << "  ";
      
      for (auto intermediate_node : path) { 
        //Rcpp::Rcout << intermediate_node->get_pid();
        
        int d = suspect->get_haplotype_L0(intermediate_node);
        
        if (d > max_L0) {
          max_L0 = d;
          //Rcpp::Rcout << "!";
        }
        
        //Rcpp::Rcout << " ";
      }
      
      //Rcpp::Rcout << std::endl;      
      
      if (meiosis_dist == -1) {
        Rcpp::stop("Cannot occur in pedigree!");
      }
      
      meiosis_dists.push_back(meiosis_dist);
      max_L0_dists.push_back(max_L0);
      pids.push_back(dest->get_pid());
    }
  }
  
  size_t n = meiosis_dists.size();
  
  Rcpp::IntegerMatrix matches(n, 3);
  colnames(matches) = Rcpp::CharacterVector::create("meioses", "max_L0", "pid");
  
  for (size_t i = 0; i < n; ++i) {
    matches(i, 0) = meiosis_dists[i];
    matches(i, 1) = max_L0_dists[i];
    matches(i, 2) = pids[i];
  }
  
  return matches;
}
  

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_matches_info(const Rcpp::XPtr<Individual> suspect, const Rcpp::List matching_indv) {
  const std::vector<bool> h = suspect->get_haplotype();
  const int suspect_pedigree_id = suspect->get_pedigree_id();

  std::vector<int> meiosis_dists;
  std::vector<int> max_L0_dists;
  std::vector<int> pids;
  
  // includes suspect by purpose
  for (int i = 0; i < matching_indv.size(); ++i) {
    Rcpp::XPtr<Individual> dest = matching_indv[i];
    
    // only considering within pedigree matches
    if (dest->get_pedigree_id() != suspect_pedigree_id) {
      continue;
    }
    
  
    std::vector<Individual*> path = suspect->calculate_path_to(dest);  
    int meiosis_dist = path.size() - 1;
    //int meiosis_dist = suspect->meiosis_dist_tree(dest);    
    //int meiosis_dist_from_path = path.size() - 1; // n vertices means n-1 edges (tree)
    //Rcpp::Rcout << ">> path from " << suspect->get_pid() << " to " << dest->get_pid() << " has length = " << meiosis_dist_from_path << " and meioses = " << meiosis_dist << (meiosis_dist_from_path == meiosis_dist ? " ok" : " ERROR") << ": " << std::endl;
    
    int max_L0 = 0;
    
    //Rcpp::Rcout << "  ";
    
    for (auto intermediate_node : path) { 
      //Rcpp::Rcout << intermediate_node->get_pid();
      
      int d = suspect->get_haplotype_L0(intermediate_node);
      
      if (d > max_L0) {
        max_L0 = d;
      }
    }
    
    if (meiosis_dist == -1) { // <=> path.size() == 0
      Rcpp::stop("Cannot occur in pedigree!");
    }
    
    meiosis_dists.push_back(meiosis_dist);
    max_L0_dists.push_back(max_L0);
    pids.push_back(dest->get_pid());
  }
  
  size_t n = meiosis_dists.size();
  
  Rcpp::IntegerMatrix matches(n, 3);
  colnames(matches) = Rcpp::CharacterVector::create("meioses", "max_L0", "pid");
  
  for (size_t i = 0; i < n; ++i) {
    matches(i, 0) = meiosis_dists[i];
    matches(i, 1) = max_L0_dists[i];
    matches(i, 2) = pids[i];
  }
  
  return matches;
}
  
  
//' @export
// [[Rcpp::export]]
int meiotic_dist(Rcpp::XPtr<Individual> ind1, Rcpp::XPtr<Individual> ind2) {
  return ind1->meiosis_dist_tree(ind2);
}

//' @export
// [[Rcpp::export]]
int count_haplotype_occurrences_pedigree(Rcpp::XPtr<Pedigree> pedigree, 
  const Rcpp::LogicalVector haplotype, 
  const int haplotype_total_no_variants,
  int generation_upper_bound_in_result = -1) {
  
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotype);

  std::vector<Individual*>* family = pedigree->get_all_individuals();

  for (auto dest : *family) {    
    int generation = dest->get_generations_from_final();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    std::vector<bool> dest_h = dest->get_haplotype();
    
    if (dest_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (dest_h == h) {
      count += 1;
    }    
  }
  
  return count;
}


//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector haplotypes_to_hashes(Rcpp::ListOf< Rcpp::LogicalVector > haplotypes) {   
  size_t n = haplotypes.size();
  std::unordered_map< std::vector<bool>, std::vector<int> > hashtable;
  
  for (size_t i = 0; i < n; ++i) {
    std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotypes[i]);
    hashtable[h].push_back(i);
  }  
  
  Rcpp::IntegerVector hap_ids(n);
  int id = 1;
  
  for (auto it : hashtable) {
    std::vector<int> indices = it.second;
    
    for (size_t j = 0; j < indices.size(); ++j) {
      hap_ids[ indices[j] ] = id;
    }
    
    ++id;
  }

  return hap_ids;
}






//' @export
// [[Rcpp::export]]
Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > build_haplotypes_hashmap(const Rcpp::List individuals) {
  int n = individuals.size();
  
  std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > >* hashtable = new std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > >();
  Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > res(hashtable, RCPP_XPTR_2ND_ARG);
  res.attr("class") = Rcpp::CharacterVector::create("mitolina_haplotype_hashmap", "externalptr");
  
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    (*hashtable)[indv->get_haplotype()].push_back(indv);
  }
  
  return res;
}


//' @export
// [[Rcpp::export]]
void print_haplotypes_hashmap(const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > hashmap) {
  std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > map = *hashmap;

  auto fn = map.hash_function();
  int n = 0;
  
  for (auto elem : map) {
    size_t hash = fn(elem.first);
    std::vector< Rcpp::XPtr<Individual> > indvs = elem.second;
    
    //Rcpp::Rcout << "Key = " << elem.first << std::endl;
    if (indvs.size() > 50) {
      Rcpp::Rcout << "Hash = " << hash << " (size = " << indvs.size() << ")" << std::endl;
    }
    
    int m = 0;
    
    for (auto indv: indvs) {
      if (indvs.size() > 50) {
        if (m == 0) {
          Rcpp::Rcout << "   ";
        }
        
        if (m <= 6) {
          Rcpp::Rcout << indv->get_pid() << " ";
        } else if (m == 7) { 
          Rcpp::Rcout << "..." << std::endl;
        }
      }
      
      m += 1;
      n += 1;      
    }
    
    //Individual* indv = elem.second;
    //Rcpp::Rcout << "Key = " << elem.first << ", pid = " << indv->get_pid() << std::endl;
    //Rcpp::Rcout << "Key = " << elem.first << std::endl;
  }
  
  Rcpp::Rcout << "Total = " << n << std::endl;
}


//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotype_matching_individuals_from_hashmap(
    const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > hashmap,
    const Rcpp::LogicalVector haplotype) {
    
  std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > map = *hashmap;  
  Rcpp::List ret_indv_empty;
  std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotype);
  
  std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > >::const_iterator got = map.find(h);

  if (got == map.end()) {
    return ret_indv_empty;
  } else {
    std::vector< Rcpp::XPtr<Individual> > indvs = got->second;
    
    //Rcpp::List ret_indv(indvs.size());
    Rcpp::List ret_indv;
    
    for (auto indv: indvs) {
      ret_indv.push_back(indv);
    }
    
    return ret_indv;
  }

  return ret_indv_empty;
}

