#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "mitolina_types.h"

//' Get haplotypes from a vector of pids.
//' 
//' Requires that haplotypes are first populated, e.g. 
//' with [pedigrees_all_populate_haplotypes()] or 
//' [pedigrees_all_populate_haplotypes_custom_founders()].
//' 
//' @param population Population
//' @param pids Vector of pids to get haplotypes for.
//' 
//' @return List of haplotypes where row `i` is the haplotype of `individuals[[i]]`.
//' 
//' @seealso [get_haplotypes_individuals()].
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotypes_pids(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids) { 
  size_t N = pids.size();
  Rcpp::List haps(N);

  // FIXME: LogicalMatrix?
  
  for (size_t i = 0; i < N; ++i) {
    Individual* indv = population->get_individual(pids[i]);
    haps(i) = indv->get_haplotype();
  }

  return haps;
}
 
//' Get haplotype matrix from list of individuals
//' 
//' Requires that haplotypes are first populated, e.g. 
//' with [pedigrees_all_populate_haplotypes()] or 
//' [pedigrees_all_populate_haplotypes_custom_founders()].
//' 
//' @param individuals Individuals to get haplotypes for.
//' @return Matrix of haplotypes where row `i` is the haplotype of `individuals[[i]]`.
//' 
//' @seealso [get_haplotypes_pids()].
//' 
//' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix get_haplotypes_individuals(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals) {  
  size_t n = individuals.size();
 
  if (n <= 0) {
    Rcpp::LogicalMatrix empty_haps(0, 0);
    return empty_haps;
  }
 
  size_t loci = individuals[0]->get_haplotype().size();

  if (loci <= 0) {
    Rcpp::stop("Expected > 0 loci");
    Rcpp::LogicalMatrix empty_haps(0, 0);
    return empty_haps;
  }

  Rcpp::LogicalMatrix haps(n, loci);

  for (size_t i = 0; i < n; ++i) {
    std::vector<bool> hap = individuals[i]->get_haplotype();

    if (hap.size() != loci) {
      Rcpp::stop("Expected > 0 loci for all haplotypes");
      Rcpp::LogicalMatrix empty_haps(0, 0);
      return empty_haps;
    }
    
    Rcpp::LogicalVector h = Rcpp::wrap(hap);
    haps(i, Rcpp::_) = h;
  }

  return haps;
}

//' Is individuals females (or males)
//' 
//' @param individuals Individuals to get haplotypes for.
//' @return Logical vector: true for female, false for male
//' 
//' @export
// [[Rcpp::export]]
Rcpp::LogicalVector get_individuals_is_female(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals) {  
  size_t n = individuals.size();
 
  if (n <= 0) {
    Rcpp::LogicalVector empty(0);
    return empty;
  }
 
  Rcpp::LogicalVector sexes(n);

  for (size_t i = 0; i < n; ++i) {
    sexes[i] = individuals[i]->is_female();
  }

  return sexes;
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


//' Populate haplotypes in pedigrees (founder types same for all).
//' 
//' Populate haplotypes from founder and down in all pedigrees.
//' Note, that haplotypes are binary (TRUE/FALSE) and 
//' that all founders get haplotype `rep(FALSE, loci)`.
//' 
//' Note, that pedigrees must first have been inferred by [build_pedigrees()].
//' 
//' @param pedigrees Pedigree list in which to populate haplotypes
//' @param loci Number of loci
//' @param mutation_rates Vector with mutation rates, length `loci`
//' @param progress Show progress
//'
//' @seealso [pedigrees_all_populate_haplotypes_custom_founders()].
//' 
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, Rcpp::NumericVector mutation_rates, bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  int loci = mut_rates.size();
  
  if (loci <= 0) {
    Rcpp::stop("At least one mutation rate / locus required");
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

//' Populate haplotypes in pedigrees (custom founder types).
//' 
//' Populate haplotypes from founder and down in all pedigrees.
//' All founders get a haplotype from calling the user 
//' provided function `get_founder_haplotype()` that must return a vector of TRUE/FALSE values.
//' 
//' Note, that pedigrees must first have been inferred by [build_pedigrees()].
//' 
//' @param pedigrees Pedigree list in which to populate haplotypes
//' @param mutation_rates Vector with mutation rates
//' @param get_founder_haplotype Function taking no arguments returning a haplotype, i.e. a logical vector (TRUE/FALSE values) of length `length(mutation_rates)`
//' @param progress Show progress
//'
//' @seealso [pedigrees_all_populate_haplotypes()].
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
    const Rcpp::LogicalVector haplotype) {
    
  int n = individuals.size();
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotype);
  int haplotype_no_variants = std::count(h.begin(), h.end(), true);
  
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

//' Get individuals matching from list of individuals
//' 
//' Get the indvididuals that matches `haplotype` in `individuals`.
//' 
//' @param individuals List of individuals to count occurrences in.
//' @param haplotype Haplotype to count occurrences of.
//' @return List of individuals that matches `haplotype` amongst `individuals`.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotype_matching_individuals(const Rcpp::List individuals, 
    const Rcpp::LogicalVector haplotype) {

  int n = individuals.size();
  int loci = haplotype.size();
  Rcpp::List ret_indv;
  
  std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotype);
  int haplotype_no_variants = std::count(h.begin(), h.end(), true);
  
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
  int generation_upper_bound_in_result = -1) {
  
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<bool> h = Rcpp::as< std::vector<bool> >(haplotype);
  int haplotype_no_variants = std::count(h.begin(), h.end(), true);
  
  std::vector<Individual*>* family = pedigree->get_all_individuals();

  for (auto dest : *family) {    
    int generation = dest->get_generations_from_final();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    std::vector<bool> dest_h = dest->get_haplotype();
    
    if (haplotype_no_variants != dest->get_haplotype_total_no_variants()) {
      continue;
    }
    
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

