context("Simulate population")

set.seed(1)
sim_res_growth <- sample_mtdna_geneology_varying_size(
  population_sizes_females = rep(1e3, 20),
  population_sizes_males = c(0L, rep(1e3, 19)),
  enable_gamma_variance_extension = TRUE,
  gamma_parameter_shape = 5,
  gamma_parameter_scale = 1/5,
  extra_generations_full = 2,
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
