// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "mitolina_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// wipe_pedigrees
void wipe_pedigrees(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP _mitolina_wipe_pedigrees(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    wipe_pedigrees(pedigrees);
    return R_NilValue;
END_RCPP
}
// build_pedigrees
Rcpp::XPtr< std::vector<Pedigree*> > build_pedigrees(Rcpp::XPtr<Population> population, bool progress);
RcppExport SEXP _mitolina_build_pedigrees(SEXP populationSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(build_pedigrees(population, progress));
    return rcpp_result_gen;
END_RCPP
}
// sample_mtdna_geneology_varying_size
List sample_mtdna_geneology_varying_size(IntegerVector population_sizes_females, IntegerVector population_sizes_males, int extra_generations_full, double gamma_parameter_shape, double gamma_parameter_scale, bool enable_gamma_variance_extension, bool progress, int extra_individuals_generations_return);
RcppExport SEXP _mitolina_sample_mtdna_geneology_varying_size(SEXP population_sizes_femalesSEXP, SEXP population_sizes_malesSEXP, SEXP extra_generations_fullSEXP, SEXP gamma_parameter_shapeSEXP, SEXP gamma_parameter_scaleSEXP, SEXP enable_gamma_variance_extensionSEXP, SEXP progressSEXP, SEXP extra_individuals_generations_returnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type population_sizes_females(population_sizes_femalesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type population_sizes_males(population_sizes_malesSEXP);
    Rcpp::traits::input_parameter< int >::type extra_generations_full(extra_generations_fullSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_parameter_shape(gamma_parameter_shapeSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_parameter_scale(gamma_parameter_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type enable_gamma_variance_extension(enable_gamma_variance_extensionSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< int >::type extra_individuals_generations_return(extra_individuals_generations_returnSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_mtdna_geneology_varying_size(population_sizes_females, population_sizes_males, extra_generations_full, gamma_parameter_shape, gamma_parameter_scale, enable_gamma_variance_extension, progress, extra_individuals_generations_return));
    return rcpp_result_gen;
END_RCPP
}
// get_haplotypes_pids
Rcpp::List get_haplotypes_pids(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids);
RcppExport SEXP _mitolina_get_haplotypes_pids(SEXP populationSEXP, SEXP pidsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pids(pidsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotypes_pids(population, pids));
    return rcpp_result_gen;
END_RCPP
}
// get_haplotypes_individuals
Rcpp::LogicalMatrix get_haplotypes_individuals(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals);
RcppExport SEXP _mitolina_get_haplotypes_individuals(SEXP individualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ListOf< Rcpp::XPtr<Individual> > >::type individuals(individualsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotypes_individuals(individuals));
    return rcpp_result_gen;
END_RCPP
}
// get_individuals_is_female
Rcpp::LogicalVector get_individuals_is_female(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals);
RcppExport SEXP _mitolina_get_individuals_is_female(SEXP individualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ListOf< Rcpp::XPtr<Individual> > >::type individuals(individualsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_individuals_is_female(individuals));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_populate_haplotypes
void pedigree_populate_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, Rcpp::NumericVector mutation_rates);
RcppExport SEXP _mitolina_pedigree_populate_haplotypes(SEXP pedSEXP, SEXP lociSEXP, SEXP mutation_ratesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mutation_rates(mutation_ratesSEXP);
    pedigree_populate_haplotypes(ped, loci, mutation_rates);
    return R_NilValue;
END_RCPP
}
// pedigrees_all_populate_haplotypes
void pedigrees_all_populate_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, Rcpp::NumericVector mutation_rates, bool progress);
RcppExport SEXP _mitolina_pedigrees_all_populate_haplotypes(SEXP pedigreesSEXP, SEXP lociSEXP, SEXP mutation_ratesSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mutation_rates(mutation_ratesSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    pedigrees_all_populate_haplotypes(pedigrees, loci, mutation_rates, progress);
    return R_NilValue;
END_RCPP
}
// pedigrees_all_populate_haplotypes_custom_founders
void pedigrees_all_populate_haplotypes_custom_founders(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, Rcpp::NumericVector mutation_rates, Rcpp::Nullable<Rcpp::Function> get_founder_haplotype, bool progress);
RcppExport SEXP _mitolina_pedigrees_all_populate_haplotypes_custom_founders(SEXP pedigreesSEXP, SEXP mutation_ratesSEXP, SEXP get_founder_haplotypeSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mutation_rates(mutation_ratesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::Function> >::type get_founder_haplotype(get_founder_haplotypeSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    pedigrees_all_populate_haplotypes_custom_founders(pedigrees, mutation_rates, get_founder_haplotype, progress);
    return R_NilValue;
END_RCPP
}
// get_haplotype
std::vector<bool> get_haplotype(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _mitolina_get_haplotype(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotype(individual));
    return rcpp_result_gen;
END_RCPP
}
// get_haplotype_no_variants
int get_haplotype_no_variants(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _mitolina_get_haplotype_no_variants(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotype_no_variants(individual));
    return rcpp_result_gen;
END_RCPP
}
// count_haplotype_occurrences_individuals
int count_haplotype_occurrences_individuals(const Rcpp::List individuals, const Rcpp::LogicalVector haplotype);
RcppExport SEXP _mitolina_count_haplotype_occurrences_individuals(SEXP individualsSEXP, SEXP haplotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type haplotype(haplotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(count_haplotype_occurrences_individuals(individuals, haplotype));
    return rcpp_result_gen;
END_RCPP
}
// get_haplotype_matching_individuals
Rcpp::List get_haplotype_matching_individuals(const Rcpp::List individuals, const Rcpp::LogicalVector haplotype);
RcppExport SEXP _mitolina_get_haplotype_matching_individuals(SEXP individualsSEXP, SEXP haplotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type haplotype(haplotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotype_matching_individuals(individuals, haplotype));
    return rcpp_result_gen;
END_RCPP
}
// get_matches_info
Rcpp::IntegerMatrix get_matches_info(const Rcpp::XPtr<Individual> suspect, const Rcpp::List matching_indv);
RcppExport SEXP _mitolina_get_matches_info(SEXP suspectSEXP, SEXP matching_indvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::XPtr<Individual> >::type suspect(suspectSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type matching_indv(matching_indvSEXP);
    rcpp_result_gen = Rcpp::wrap(get_matches_info(suspect, matching_indv));
    return rcpp_result_gen;
END_RCPP
}
// meiotic_dist
int meiotic_dist(Rcpp::XPtr<Individual> ind1, Rcpp::XPtr<Individual> ind2);
RcppExport SEXP _mitolina_meiotic_dist(SEXP ind1SEXP, SEXP ind2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type ind2(ind2SEXP);
    rcpp_result_gen = Rcpp::wrap(meiotic_dist(ind1, ind2));
    return rcpp_result_gen;
END_RCPP
}
// count_haplotype_occurrences_pedigree
int count_haplotype_occurrences_pedigree(Rcpp::XPtr<Pedigree> pedigree, const Rcpp::LogicalVector haplotype, int generation_upper_bound_in_result);
RcppExport SEXP _mitolina_count_haplotype_occurrences_pedigree(SEXP pedigreeSEXP, SEXP haplotypeSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type pedigree(pedigreeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type haplotype(haplotypeSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(count_haplotype_occurrences_pedigree(pedigree, haplotype, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// haplotypes_to_hashes
Rcpp::IntegerVector haplotypes_to_hashes(Rcpp::ListOf< Rcpp::LogicalVector > haplotypes);
RcppExport SEXP _mitolina_haplotypes_to_hashes(SEXP haplotypesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ListOf< Rcpp::LogicalVector > >::type haplotypes(haplotypesSEXP);
    rcpp_result_gen = Rcpp::wrap(haplotypes_to_hashes(haplotypes));
    return rcpp_result_gen;
END_RCPP
}
// build_haplotypes_hashmap
Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > build_haplotypes_hashmap(const Rcpp::List individuals);
RcppExport SEXP _mitolina_build_haplotypes_hashmap(SEXP individualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type individuals(individualsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_haplotypes_hashmap(individuals));
    return rcpp_result_gen;
END_RCPP
}
// print_haplotypes_hashmap
void print_haplotypes_hashmap(const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > hashmap);
RcppExport SEXP _mitolina_print_haplotypes_hashmap(SEXP hashmapSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > >::type hashmap(hashmapSEXP);
    print_haplotypes_hashmap(hashmap);
    return R_NilValue;
END_RCPP
}
// get_haplotype_matching_individuals_from_hashmap
Rcpp::List get_haplotype_matching_individuals_from_hashmap(const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > hashmap, const Rcpp::LogicalVector haplotype);
RcppExport SEXP _mitolina_get_haplotype_matching_individuals_from_hashmap(SEXP hashmapSEXP, SEXP haplotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::XPtr< std::unordered_map< std::vector<bool>, std::vector< Rcpp::XPtr<Individual> > > > >::type hashmap(hashmapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type haplotype(haplotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotype_matching_individuals_from_hashmap(hashmap, haplotype));
    return rcpp_result_gen;
END_RCPP
}
// get_individual
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid);
RcppExport SEXP _mitolina_get_individual(SEXP populationSEXP, SEXP pidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type pid(pidSEXP);
    rcpp_result_gen = Rcpp::wrap(get_individual(population, pid));
    return rcpp_result_gen;
END_RCPP
}
// get_pid
int get_pid(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _mitolina_get_pid(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pid(individual));
    return rcpp_result_gen;
END_RCPP
}
// print_individual
void print_individual(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _mitolina_print_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    print_individual(individual);
    return R_NilValue;
END_RCPP
}
// get_generations_from_final
int get_generations_from_final(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _mitolina_get_generations_from_final(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_generations_from_final(individual));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_from_individual
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _mitolina_get_pedigree_from_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_from_individual(individual));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_id_from_pid
Rcpp::IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids);
RcppExport SEXP _mitolina_get_pedigree_id_from_pid(SEXP populationSEXP, SEXP pidsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pids(pidsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_id_from_pid(population, pids));
    return rcpp_result_gen;
END_RCPP
}
// pop_size
int pop_size(Rcpp::XPtr<Population> population);
RcppExport SEXP _mitolina_pop_size(SEXP populationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    rcpp_result_gen = Rcpp::wrap(pop_size(population));
    return rcpp_result_gen;
END_RCPP
}
// get_individuals
Rcpp::ListOf< Rcpp::XPtr<Individual> > get_individuals(Rcpp::XPtr<Population> population);
RcppExport SEXP _mitolina_get_individuals(SEXP populationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    rcpp_result_gen = Rcpp::wrap(get_individuals(population));
    return rcpp_result_gen;
END_RCPP
}
// meioses_generation_distribution
Rcpp::IntegerMatrix meioses_generation_distribution(Rcpp::XPtr<Individual> individual, int generation_upper_bound_in_result);
RcppExport SEXP _mitolina_meioses_generation_distribution(SEXP individualSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(meioses_generation_distribution(individual, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// population_size_generation
int population_size_generation(Rcpp::XPtr<Population> population, int generation_upper_bound_in_result);
RcppExport SEXP _mitolina_population_size_generation(SEXP populationSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(population_size_generation(population, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_size_generation
int pedigree_size_generation(Rcpp::XPtr<Pedigree> pedigree, int generation_upper_bound_in_result);
RcppExport SEXP _mitolina_pedigree_size_generation(SEXP pedigreeSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type pedigree(pedigreeSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_size_generation(pedigree, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_id
int get_pedigree_id(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_get_pedigree_id(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_id(ped));
    return rcpp_result_gen;
END_RCPP
}
// pedigrees_count
int pedigrees_count(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP _mitolina_pedigrees_count(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigrees_count(pedigrees));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_size
int pedigree_size(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_pedigree_size(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_size(ped));
    return rcpp_result_gen;
END_RCPP
}
// pedigrees_table
std::unordered_map<int, int> pedigrees_table(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP _mitolina_pedigrees_table(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigrees_table(pedigrees));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree
Rcpp::XPtr<Pedigree> get_pedigree(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int index);
RcppExport SEXP _mitolina_get_pedigree(SEXP pedigreesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree(pedigrees, index));
    return rcpp_result_gen;
END_RCPP
}
// print_pedigree
void print_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_print_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    print_pedigree(ped);
    return R_NilValue;
END_RCPP
}
// get_pids_in_pedigree
Rcpp::IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_get_pids_in_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pids_in_pedigree(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_is_female_in_pedigree
Rcpp::LogicalVector get_is_female_in_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_get_is_female_in_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_is_female_in_pedigree(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_haplotypes_in_pedigree
Rcpp::List get_haplotypes_in_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_get_haplotypes_in_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotypes_in_pedigree(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_edgelist
Rcpp::CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_get_pedigree_edgelist(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_edgelist(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_as_graph
Rcpp::List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _mitolina_get_pedigree_as_graph(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_as_graph(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigrees_tidy
Rcpp::List get_pedigrees_tidy(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP _mitolina_get_pedigrees_tidy(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigrees_tidy(pedigrees));
    return rcpp_result_gen;
END_RCPP
}
// test_create_population
Rcpp::XPtr<Population> test_create_population();
RcppExport SEXP _mitolina_test_create_population() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(test_create_population());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mitolina_wipe_pedigrees", (DL_FUNC) &_mitolina_wipe_pedigrees, 1},
    {"_mitolina_build_pedigrees", (DL_FUNC) &_mitolina_build_pedigrees, 2},
    {"_mitolina_sample_mtdna_geneology_varying_size", (DL_FUNC) &_mitolina_sample_mtdna_geneology_varying_size, 8},
    {"_mitolina_get_haplotypes_pids", (DL_FUNC) &_mitolina_get_haplotypes_pids, 2},
    {"_mitolina_get_haplotypes_individuals", (DL_FUNC) &_mitolina_get_haplotypes_individuals, 1},
    {"_mitolina_get_individuals_is_female", (DL_FUNC) &_mitolina_get_individuals_is_female, 1},
    {"_mitolina_pedigree_populate_haplotypes", (DL_FUNC) &_mitolina_pedigree_populate_haplotypes, 3},
    {"_mitolina_pedigrees_all_populate_haplotypes", (DL_FUNC) &_mitolina_pedigrees_all_populate_haplotypes, 4},
    {"_mitolina_pedigrees_all_populate_haplotypes_custom_founders", (DL_FUNC) &_mitolina_pedigrees_all_populate_haplotypes_custom_founders, 4},
    {"_mitolina_get_haplotype", (DL_FUNC) &_mitolina_get_haplotype, 1},
    {"_mitolina_get_haplotype_no_variants", (DL_FUNC) &_mitolina_get_haplotype_no_variants, 1},
    {"_mitolina_count_haplotype_occurrences_individuals", (DL_FUNC) &_mitolina_count_haplotype_occurrences_individuals, 2},
    {"_mitolina_get_haplotype_matching_individuals", (DL_FUNC) &_mitolina_get_haplotype_matching_individuals, 2},
    {"_mitolina_get_matches_info", (DL_FUNC) &_mitolina_get_matches_info, 2},
    {"_mitolina_meiotic_dist", (DL_FUNC) &_mitolina_meiotic_dist, 2},
    {"_mitolina_count_haplotype_occurrences_pedigree", (DL_FUNC) &_mitolina_count_haplotype_occurrences_pedigree, 3},
    {"_mitolina_haplotypes_to_hashes", (DL_FUNC) &_mitolina_haplotypes_to_hashes, 1},
    {"_mitolina_build_haplotypes_hashmap", (DL_FUNC) &_mitolina_build_haplotypes_hashmap, 1},
    {"_mitolina_print_haplotypes_hashmap", (DL_FUNC) &_mitolina_print_haplotypes_hashmap, 1},
    {"_mitolina_get_haplotype_matching_individuals_from_hashmap", (DL_FUNC) &_mitolina_get_haplotype_matching_individuals_from_hashmap, 2},
    {"_mitolina_get_individual", (DL_FUNC) &_mitolina_get_individual, 2},
    {"_mitolina_get_pid", (DL_FUNC) &_mitolina_get_pid, 1},
    {"_mitolina_print_individual", (DL_FUNC) &_mitolina_print_individual, 1},
    {"_mitolina_get_generations_from_final", (DL_FUNC) &_mitolina_get_generations_from_final, 1},
    {"_mitolina_get_pedigree_from_individual", (DL_FUNC) &_mitolina_get_pedigree_from_individual, 1},
    {"_mitolina_get_pedigree_id_from_pid", (DL_FUNC) &_mitolina_get_pedigree_id_from_pid, 2},
    {"_mitolina_pop_size", (DL_FUNC) &_mitolina_pop_size, 1},
    {"_mitolina_get_individuals", (DL_FUNC) &_mitolina_get_individuals, 1},
    {"_mitolina_meioses_generation_distribution", (DL_FUNC) &_mitolina_meioses_generation_distribution, 2},
    {"_mitolina_population_size_generation", (DL_FUNC) &_mitolina_population_size_generation, 2},
    {"_mitolina_pedigree_size_generation", (DL_FUNC) &_mitolina_pedigree_size_generation, 2},
    {"_mitolina_get_pedigree_id", (DL_FUNC) &_mitolina_get_pedigree_id, 1},
    {"_mitolina_pedigrees_count", (DL_FUNC) &_mitolina_pedigrees_count, 1},
    {"_mitolina_pedigree_size", (DL_FUNC) &_mitolina_pedigree_size, 1},
    {"_mitolina_pedigrees_table", (DL_FUNC) &_mitolina_pedigrees_table, 1},
    {"_mitolina_get_pedigree", (DL_FUNC) &_mitolina_get_pedigree, 2},
    {"_mitolina_print_pedigree", (DL_FUNC) &_mitolina_print_pedigree, 1},
    {"_mitolina_get_pids_in_pedigree", (DL_FUNC) &_mitolina_get_pids_in_pedigree, 1},
    {"_mitolina_get_is_female_in_pedigree", (DL_FUNC) &_mitolina_get_is_female_in_pedigree, 1},
    {"_mitolina_get_haplotypes_in_pedigree", (DL_FUNC) &_mitolina_get_haplotypes_in_pedigree, 1},
    {"_mitolina_get_pedigree_edgelist", (DL_FUNC) &_mitolina_get_pedigree_edgelist, 1},
    {"_mitolina_get_pedigree_as_graph", (DL_FUNC) &_mitolina_get_pedigree_as_graph, 1},
    {"_mitolina_get_pedigrees_tidy", (DL_FUNC) &_mitolina_get_pedigrees_tidy, 1},
    {"_mitolina_test_create_population", (DL_FUNC) &_mitolina_test_create_population, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_mitolina(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
