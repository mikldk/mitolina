context("Simulate population")

set.seed(1)
sim_res_growth <- sample_mtdna_geneology_varying_size(
  population_sizes_females = rep(1e3, 20),
  population_sizes_males = c(0L, rep(1e3, 19)),
  generations_full = 3,
  enable_gamma_variance_extension = TRUE,
  gamma_parameter_shape = 5,
  gamma_parameter_scale = 1/5,
  progress = FALSE)

test_that("sample_geneology works", {
  expect_failure(expect_null(sim_res_growth))
  expect_output(print(sim_res_growth$population), regexp = "^Population with .* individuals$")
  expect_equal(length(sim_res_growth$end_generation_female_individuals), 1000L)
  expect_equal(length(sim_res_growth$end_generation_male_individuals), 1000L)
  expect_equal(length(sim_res_growth$female_individuals_generations), 3L*1000L)
  expect_equal(length(sim_res_growth$male_individuals_generations), 3L*1000L)
  expect_equal(sim_res_growth$sdo_type, "GammaVariation")
})


peds <- build_pedigrees(sim_res_growth$population, progress = FALSE)
test_that("build_pedigrees", {
  expect_true(pedigrees_count(peds) > 0)
})


sites <- 1000L
get_founder_mito <- function() {
  sample(c(FALSE, TRUE), sites, replace = TRUE)
}
mu <- rep(1e-2, sites)
pedigrees_all_populate_haplotypes_custom_founders(
  pedigrees = peds,
  mutation_rates = mu,
  get_founder_haplotype = get_founder_mito, 
  progress = FALSE)

haps_3gen <- get_haplotypes_individuals(sim_res_growth$male_individuals_generations)
pids_3gen <- unlist(lapply(sim_res_growth$male_individuals_generations, get_pid))
haps_pids_3gen <- get_haplotypes_pids(sim_res_growth$population, pids_3gen)
test_that("male_individuals_generations: get_haplotypes_individuals / get_haplotypes_pids haps", {
  expect_equal(nrow(haps_3gen), length(pids_3gen))
  expect_equal(nrow(haps_3gen), length(haps_pids_3gen))
  expect_equal(haps_3gen, do.call(rbind, haps_pids_3gen))
})

ped <- peds[[1]]
haps_ped <- get_haplotypes_in_pedigree(ped)
pids_ped <- get_pids_in_pedigree(ped)
haps_pids_ped <- get_haplotypes_pids(sim_res_growth$population, pids_ped)
test_that("ped: get_haplotypes_in_pedigree / get_haplotypes_pids haps", {
  expect_equal(length(haps_ped), length(pids_ped))
  expect_equal(length(haps_ped), length(haps_pids_ped))
  expect_equal(haps_ped, haps_pids_ped)
})

indv_ped <- lapply(pids_ped, get_individual, population = sim_res_growth$population)
haps_indv_ped <- get_haplotypes_individuals(indv_ped)
test_that("ped: get_individual / get_haplotypes_individuals / get_haplotypes_in_pedigree haps", {
  expect_equal(nrow(haps_indv_ped), length(haps_ped))
  expect_equal(nrow(haps_indv_ped), length(pids_ped))
  expect_equal(haps_indv_ped, do.call(rbind, haps_ped))
})

