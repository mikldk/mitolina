#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]

#include <progress.hpp>

#include <string>

#include "mitolina_types.h"
#include "api_utility_pedigree.h"

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
//' @return List of haplotypes where element `i` is the haplotype of `individuals[[i]]`.
//' 
//' @seealso [get_haplotypes_pids()].
//' 
//' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix get_haplotypes_individuals(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals) {  
  size_t n = individuals.size();
 
  if (n <= 0) {
    Rcpp::LogicalMatrix empty_haps;
    return empty_haps;
  }
 
  size_t loci = individuals[0]->get_haplotype().size();

  if (loci <= 0) {
    Rcpp::stop("Expected > 0 loci");
    Rcpp::LogicalMatrix empty_haps;
    return empty_haps;
  }

  //Rcpp::List haps(n);
  //Rcpp::LogicalMatrix haps(n, loci);
  Rcpp::LogicalMatrix haps_trans(loci, n);

  for (size_t i = 0; i < n; ++i) {
    std::vector<bool> hap = individuals[i]->get_haplotype();

    if (hap.size() != loci) {
      Rcpp::stop("Expected > 0 loci for all haplotypes");
      Rcpp::LogicalMatrix empty_haps;
      return empty_haps;
    }
    
    Rcpp::LogicalVector h = Rcpp::wrap(hap);
    //haps[i] = h;
    haps_trans(Rcpp::_, i) = h;
  }
  
  //Rcpp::LogicalMatrix haps = Rcpp::transpose(haps_trans);
  Rcpp::LogicalMatrix haps = Rcpp::tranpose_impl<LGLSXP, Rcpp::PreserveStorage>(haps_trans);

  return haps;
}



//' Is individuals females (or males)
//' 
//' @param individuals Individuals to get sex for.
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


//' Is individual female (or male)
//' 
//' @param individual Individual to get sex for.
//' @return Logical: true for female, false for male
//' 
//' @export
// [[Rcpp::export]]
Rcpp::LogicalVector get_individual_is_female(Rcpp::XPtr<Individual> individual) {  
  return individual->is_female();
}

/*
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
*/

//' Populate haplotypes in pedigrees (founder types same for all).
//' 
//' Populate haplotypes from founder and down in all pedigrees.
//' Note, that haplotypes are binary (TRUE/FALSE) and 
//' that all founders get haplotype `rep(FALSE, loci)`.
//' 
//' Note, that pedigrees must first have been inferred by [build_pedigrees()].
//' 
//' @param pedigrees Pedigree list in which to populate haplotypes
//' @param mutation_rates Vector with mutation rates, length `loci`
//' @param progress Show progress?
//'
//' @seealso [pedigrees_all_populate_haplotypes_custom_founders()].
//' 
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes(Rcpp::List pedigrees, 
                                       Rcpp::NumericVector mutation_rates, 
                                       bool progress = true) {
  
  stopifnot_mitolina_pedigreelist(pedigrees);

  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  int loci = mut_rates.size();
  
  if (loci <= 0) {
    Rcpp::stop("At least one mutation rate / locus required");
  }
  
  size_t N = pedigrees.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    Rcpp::XPtr<Pedigree> ped = pedigrees.at(i);
    ped->populate_haplotypes(loci, mut_rates);
    
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
void pedigrees_all_populate_haplotypes_custom_founders(Rcpp::List pedigrees, 
                                       Rcpp::NumericVector mutation_rates,
                                       Rcpp::Nullable<Rcpp::Function> get_founder_haplotype = R_NilValue,
                                       bool progress = true) {

  stopifnot_mitolina_pedigreelist(pedigrees);
  
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  
  if (get_founder_haplotype.isNull()) {
    Rcpp::stop("get_founder_haplotype must not be NULL");
  }  
  
  Rcpp::Function g_founder_hap = Rcpp::as<Rcpp::Function>(get_founder_haplotype);

  size_t N = pedigrees.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    Rcpp::XPtr<Pedigree> ped = pedigrees.at(i);
    ped->populate_haplotypes_custom_founders(mut_rates, g_founder_hap);
    
     if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}

//' Get haplotype from an individual
//' 
//' Requires that haplotypes are first populated, e.g. 
//' with [pedigrees_all_populate_haplotypes()] or 
//' [pedigrees_all_populate_haplotypes_custom_founders()].
//' 
//' @param individual Individual to get haplotypes for.
//' @return Haplotype for `individual`.
//' 
//' @seealso [get_haplotypes_individuals()] and [get_haplotypes_pids()].
//' 
//' @export
// [[Rcpp::export]]
std::vector<bool> get_haplotype(Rcpp::XPtr<Individual> individual) {
  return individual->get_haplotype();
}

//' Get number of variants in haplotype
//'
//' Number of variants is for example faster when checking for equality or when 
//' summarising (like plotting). If haplotypes do not have same number of variants, they cannot be equal.
//'
//' @param individual Individual to get number of variants for
//'
//' @export
// [[Rcpp::export]]
int get_haplotype_no_variants(Rcpp::XPtr<Individual> individual) {
  return individual->get_haplotype_total_no_variants();
}

//' Count haplotypes occurrences in list of individuals
//' 
//' Counts the number of types `haplotype` appears in `individuals`.
//' 
//' @param individuals List of individuals to count occurrences in.
//' @param haplotype Haplotype to count occurrences of.
//' 
//' @return Number of times that `haplotype` occurred amongst `individuals`.
//' 
//' @seealso [get_matches_info()], [count_haplotype_occurrences_individuals()].
//' 
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

//' Information about matching individuals
//' 
//' Note: This function does not check that individuals in 
//' `matching_indv` actually match.
//'
//' Note: only considering individuals within same pedigree!
//'  
//' This gives detailed information about matching individuals in the pedigree, 
//' e.g. meiotic distances and maximum L0 distance (number of sites they differ) on the path as some of these 
//' matches may have (back)mutations between in between them (but often this will be 0).
//' 
//' @param suspect Individual that others must match the profile of.
//' @param matching_indv List of matching individuals to get information for.
//' 
//' @return Matrix with information about matching individuals. 
//' Columns in order: meioses (meiotic distance to `suspect`), 
//' max_L0 (on the path between the matching individual and `suspect`, 
//' what is the maximum L0 distance between the `suspect`'s profile and the 
//' profiles of the individuals on the path), 
//' pid (pid of matching individual)
//' 
//' @seealso [count_haplotype_occurrences_individuals()] and [count_haplotype_occurrences_pedigree()].
//'
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
  
  
//' Meiotic distance between two individuals
//' 
//' Get the number of meioses between two individuals.
//' Note, that pedigrees must first have been inferred by [build_pedigrees()].
//' 
//' @param ind1 Individual 1
//' @param ind2 Individual 2
//' 
//' @return Number of meioses between `ind1` and `ind2` if they are in the same pedigree, else -1.
//' 
//' @export
// [[Rcpp::export]]
int meiotic_dist(Rcpp::XPtr<Individual> ind1, Rcpp::XPtr<Individual> ind2) {
  return ind1->meiosis_dist_tree(ind2);
}

//' Count haplotypes occurrences in pedigree
//' 
//' Counts the number of types `haplotype` appears in `pedigree`.
//' 
//' @param pedigree Pedigree to count occurrences in.
//' @param haplotype Haplotype to count occurrences of.
//' @param generation_upper_bound_in_result Only consider matches in 
//' generation 0, 1, ... generation_upper_bound_in_result.
//' -1 means disabled, consider all generations.
//' End generation is generation 0.
//' Second last generation is 1. 
//' And so on.
//' 
//' @return Number of times that `haplotype` occurred in `pedigree`.
//' 
//' @seealso [get_matches_info()], [count_haplotype_occurrences_individuals()].
//' 
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

//' Convert haplotypes to hashes (integers)
//' 
//' Individuals with the same haplotype will have the same hash (integer)
//' and individuals with different haplotypes will have different hashes (integers).
//' 
//' This can be useful if for example using haplotypes to define groups 
//' and the haplotype itself is not of interest.
//' 
//' @param haplotypes List of haplotypes (list of logical vectors)
//' 
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




/*

//' Build hashmap of haplotypes to individuals
//' 
//' Makes it possible to find all individuals with a certain haplotype.
//' Must be used with e.g. [get_haplotype_matching_individuals_from_hashmap()] 
//' or [print_haplotypes_hashmap()].
//' 
//' @param individuals List of individuals to build hashmap of
//' @param max_load_factor Tuning parameter for hash table
//' @param verbose_interval 0 for no verbose output, e.g. 1,000 for output for every 1000 individual added
//' @return Hashmap with haplotypes as keys and vector of individuals as value
//' 
//' @seealso [get_haplotype_matching_individuals_from_hashmap()] 
//' and [print_haplotypes_hashmap()].
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > build_haplotypes_hashmap(
    const Rcpp::List& individuals, 
    const float max_load_factor = 10,
    const int verbose_interval = 0) {
  
  int n = individuals.size();
  
  std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > >* hashtable = new std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > >();
  hashtable->reserve(n);
  hashtable->max_load_factor(max_load_factor);
  
  for (int i = 0; i < n; ++i) {
    if (verbose_interval > 0 && i % verbose_interval == 0) {
      Rcpp::Rcout << 100.0*((double)(i + 1)/(double)n) << "%:" << std::endl;
      Rcpp::Rcout << "  current max_load_factor: " << hashtable->max_load_factor() << std::endl;
      Rcpp::Rcout << "  current max_load_factor: " << hashtable->max_load_factor() << std::endl;
      Rcpp::Rcout << "  current size: " << hashtable->size() << std::endl;
      Rcpp::Rcout << "  current bucket_count: " << hashtable->bucket_count() << std::endl;
      Rcpp::Rcout << "  current load_factor: " << hashtable->load_factor() << std::endl;    
    }
    
    Rcpp::XPtr<Individual> indv = individuals[i];
    (*hashtable)[indv->get_haplotype()].push_back(indv);
  }

  Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > res(hashtable, RCPP_XPTR_2ND_ARG);
  res.attr("class") = Rcpp::CharacterVector::create("mitolina_haplotype_hashmap", "externalptr");
  
  return res;
}


//' Print haplotype hashmap
//' 
//' Print hashmap a haplotypes to individuals made by [build_haplotypes_hashmap()].
//' 
//' @param hashmap Hashmap made by [build_haplotypes_hashmap()]
//' 
//' @seealso [get_haplotype_matching_individuals_from_hashmap()] 
//' and [build_haplotypes_hashmap()].
//' 
//' @export
// [[Rcpp::export]]
void print_haplotypes_hashmap(const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > hashmap) {
  std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > >& map = *hashmap;

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


//' Get individuals with a certain haplotype by hashmap lookup
//' 
//' By using hashmap made by [build_haplotypes_hashmap()], 
//' it is easy to get all individuals with a certain haplotype.
//' 
//' @param hashmap Hashmap to make lookup in, made by [build_haplotypes_hashmap()]
//' @param haplotype to get individuals that has this haplotype
//' 
//' @return List of individuals with a given haplotype
//' 
//' @seealso [print_haplotypes_hashmap()] 
//' and [build_haplotypes_hashmap()].
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotype_matching_individuals_from_hashmap(
    const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > >& hashmap,
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




//' Delete haplotype hashmap
//' 
//' Delete hashmap made by [build_haplotypes_hashmap()].
//' 
//' @param hashmap Hashmap made by [build_haplotypes_hashmap()]
//' 
//' @seealso [get_haplotype_matching_individuals_from_hashmap()] 
//' and [build_haplotypes_hashmap()].
//' 
//' @export
// [[Rcpp::export]]
void delete_haplotypes_hashmap(Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > hashmap) {
  std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > >* map = hashmap;
  delete map;
}


*/












//' Infer haplotype ids
//' 
//' Makes it faster to compare haplotypes by an id instead of the entire mitogenome.
//' 
//' @param individuals List of individuals to infer haplotype ids for
//' @param progress Show progress?
//' 
//' @export
// [[Rcpp::export]]
void infer_haplotype_ids(const Rcpp::List& individuals, bool progress = true) {
  int n = individuals.size();
  Progress p(n, progress);
  
  std::unordered_map< std::vector<bool>, int > haplo_to_id;
  int hap_id = 1;
  
  for (size_t i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    
    if (!(indv->is_haplotype_set())) {
      Rcpp::stop("Haplotype not yet set");
    }
    
    std::vector<bool> h = indv->get_haplotype();
    
    std::unordered_map<std::vector<bool>, int>::const_iterator got = haplo_to_id.find(h);
    
    if (got == haplo_to_id.end()) {
      // not found
      indv->set_haplotype_id(hap_id);
      haplo_to_id[h] = hap_id;
      ++hap_id;
    } else {
      // found
      indv->set_haplotype_id(got->second);
    }
    
    if (progress) {
      p.increment();
    }
  }
}



//' Get haplotype id from individual
//' 
//' @param individual Individual to get haplotypes for.
//' 
//' @return Haplotype id
//' 
//' @export
// [[Rcpp::export]]
int get_haplotype_id_individual(Rcpp::XPtr<Individual> individual) {  
  return individual->get_haplotype_id();
}


//' Get haplotype ids from list of individuals
//' 
//' @param individuals Individuals to get haplotypes for.
//' 
//' @return List of haplotype ids where element `i` is the haplotype of `individuals[[i]]`.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_haplotype_ids_individuals(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals) {  
  size_t n = individuals.size();
 
  if (n <= 0) {
    Rcpp::IntegerVector empty_hapids;
    return empty_hapids;
  }
 
  Rcpp::IntegerVector hapids(n);

  for (size_t i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    hapids[i] = indv->get_haplotype_id();
  }  

  return hapids;
}




//' Build hashmap of haplotype ids to individuals
//' 
//' Makes it possible to find all individuals with a certain haplotype id.
//' Must be used with e.g. [get_haplotypeid_matching_individuals_from_hashmap()].
//' 
//' @param individuals List of individuals to build hashmap of
//' @param progress Show progress?
//' 
//' @return Hashmap with haplotype id as keys and vector of individuals as value
//' 
//' @seealso [get_haplotypeid_matching_individuals_from_hashmap()].
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< std::vector< std::vector<Individual*>* > > build_haplotypeids_hashmap(
    const Rcpp::List& individuals, bool progress = true) {
  
  
  int n = individuals.size();
  Progress p(2*n, progress);
  
  int max_id = 0;
  
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    int id = indv->get_haplotype_id();
    
    if (id > max_id) {
      max_id = id;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  // size + 1: max_id will be attained, 0 index not used
  std::vector< std::vector<Individual*>* >* hashtable = new std::vector< std::vector<Individual*>* >(max_id + 1);
  for (int id = 0; id <= max_id; ++id) {
    std::vector<Individual*>* vec = new std::vector<Individual*>();
    hashtable->at(id) = vec;
  }
    
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv_xptr = individuals[i];
    Individual* indv = indv_xptr;
    int id = indv->get_haplotype_id();
    
    hashtable->at(id)->push_back(indv);
    
    if (progress) {
      p.increment();
    }
  }

  Rcpp::XPtr< std::vector< std::vector<Individual*>* > > res(hashtable, RCPP_XPTR_2ND_ARG);
  res.attr("class") = Rcpp::CharacterVector::create("mitolina_haplotypeid_hashmap", "externalptr");
  
  return res;
}


//' Get individuals with a certain haplotype id by hashmap lookup
//' 
//' By using hashmap made by [build_haplotypeids_hashmap()], 
//' it is easy to get all individuals with a certain haplotype id.
//' 
//' @param hashmap Hashmap to make lookup in, made by [build_haplotypeids_hashmap()]
//' @param haplotype_id to get individuals that has this haplotype id
//' 
//' @return List of individuals with a given haplotype id
//' 
//' @seealso [build_haplotypeids_hashmap()].
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotypeid_matching_individuals_from_hashmap(
    const Rcpp::XPtr< std::vector< std::vector<Individual*>* > >& hashmap,
    const int haplotype_id) {
  
  Rcpp::List ret_indv_empty;  
  
  if (haplotype_id >= hashmap->size()) {
    return ret_indv_empty;
  }
  
  std::vector<Individual*>* indvs_with_id = hashmap->at(haplotype_id);
  int n = indvs_with_id->size();  
  Rcpp::List ret_indv(n);
  
  for (int i = 0; i < n; ++i) {
    Individual* indv = indvs_with_id->at(i);    
    Rcpp::XPtr<Individual> indv_xptr(indv, RCPP_XPTR_2ND_ARG);
    ret_indv[i] = indv_xptr;
  }
  
  return ret_indv;  
}




//' Delete haplotype id hashmap
//' 
//' Delete hashmap made by [build_haplotypeids_hashmap()].
//' 
//' @param hashmap Hashmap made by [build_haplotypeids_hashmap()]
//' 
//' @seealso [get_haplotypeid_matching_individuals_from_hashmap()] 
//' and [build_haplotypeids_hashmap()].
//' 
//' @export
// [[Rcpp::export]]
void delete_haplotypeids_hashmap(Rcpp::XPtr< std::vector< std::vector<Individual*>* > > hashmap) {
  std::vector< std::vector<Individual*>* >* map = hashmap;
  
  int n = map->size();  
  
  for (int i = 0; i < n; ++i) {
    std::vector<Individual*>* map_i = map->at(i);
    delete map_i;
  }
  
  delete map;
}

