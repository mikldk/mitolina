context("Pedigrees and haplotypes")

test_pop <- test_create_population()

test_that("test_create_population works", {
  expect_failure(expect_null(test_pop))
  expect_output(print(test_pop), regexp = "^Population with 19 individuals$")
  expect_equal(pop_size(test_pop), 19L)
})

indvs <- get_individuals(test_pop)
test_that("get_individuals works", {
  expect_failure(expect_null(indvs))
  expect_equal(length(indvs), 19L)
})

peds <- build_pedigrees(test_pop, progress = FALSE)
test_that("build_pedigrees works", {
  expect_output(print(peds), regexp = "^List of 2 pedigrees \\(of size 12, 7\\)$")
  expect_equal(pedigrees_count(peds), 2L)
})

if (FALSE) {
  plot(peds[[1]])
  plot(peds[[2]])
}

ped <- peds[[1]]
pids <- sort(get_pids_in_pedigree(ped))
test_that("pedigree pids works", {
  expect_equal(length(pids), 12L)
  expect_true(all(pids == c(1L:11L, 19L)))
  expect_equal(length(get_pids_in_pedigree(peds[[2]])), 7L)
})

test_that("meiotic_dist works", {
  expect_equal(0L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 1L)))
  
  expect_equal(1L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 6L)))
  
  expect_equal(4L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 10L)))
  
  expect_equal(-1L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                 get_individual(test_pop, pid = 12L)))
  
  
  expect_equal(2L, meiotic_dist(get_individual(test_pop, pid = 15L), 
                                get_individual(test_pop, pid = 16L)))
})


LOCI <- 100L
pedigrees_all_populate_haplotypes(peds, mutation_rates = rep(0L, LOCI), progress = FALSE)
test_that("pedigrees_all_populate_haplotypes works", {
  expect_output(print(peds), regexp = "^List of 2 pedigrees \\(of size 12, 7\\)$")
})

test_that("get_haplotype_matching_individuals works", {
  expect_equal(length(indvs), length(get_haplotype_matching_individuals(indvs, rep(FALSE, LOCI))))
  expect_equal(0L, length(get_haplotype_matching_individuals(indvs, rep(TRUE, LOCI))))
  expect_equal(lapply(indvs, get_pid), lapply(get_haplotype_matching_individuals(indvs, rep(FALSE, LOCI)), get_pid))
})

test_that("count_haplotype_occurrences_individuals works", {
  expect_equal(19L, count_haplotype_occurrences_individuals(indvs, rep(FALSE, LOCI)))
  expect_equal(0L, count_haplotype_occurrences_individuals(indvs, rep(TRUE, LOCI)))
  
  expect_equal(12L, count_haplotype_occurrences_pedigree(ped, rep(FALSE, LOCI)))
  expect_equal(5L, count_haplotype_occurrences_pedigree(ped, 
                                                        rep(FALSE, LOCI), 
                                                        generation_upper_bound_in_result = 0L))
  expect_equal(5L+4L, count_haplotype_occurrences_pedigree(ped, 
                                                           rep(FALSE, LOCI), 
                                                           generation_upper_bound_in_result = 1L))
  expect_equal(5L+4L+2L, count_haplotype_occurrences_pedigree(ped, 
                                                              rep(FALSE, LOCI), 
                                                              generation_upper_bound_in_result = 2L))
})

indvs_is_female <- get_individuals_is_female(indvs)
indv_females <- indvs[indvs_is_female]
indv_males <- indvs[!indvs_is_female]

test_that("get_individuals_is_female works", {
  expect_equal(10L, length(indv_females))
  expect_equal(9L, length(indv_males))
})

suspect <- get_individual(test_pop, pid = 1L)
match_females <- get_haplotype_matching_individuals(individuals = indv_females, 
                                                    haplotype = get_haplotype(suspect))
match_males <- get_haplotype_matching_individuals(individuals = indv_males, 
                                                  haplotype = get_haplotype(suspect))

test_that("get_haplotype_matching_individuals works", {
  expect_equal(length(match_females), length(indv_females))
  expect_equal(length(match_males), length(indv_males))
})

match_info_females <- get_matches_info(suspect = suspect, matching_indv = match_females)
match_info_males <- get_matches_info(suspect = suspect, matching_indv = match_males)

mei_res_females <- match_info_females[order(match_info_females[, 3L]), ] # order by pid
meioses_females <- mei_res_females[, 1L]
mei_res_males <- match_info_males[order(match_info_males[, 3L]), ] # order by pid
meioses_males <- mei_res_males[, 1L]

test_that("get_matches_info works", {
  expect_equal(mei_res_females[, 3L], c(2L, 6L:11L))
  expect_true(all(mei_res_females[, 2L] == 0L)) # max L0 == 0
  
  expect_equal(mei_res_males[, 3L], c(1L, 3L, 4L, 5L, 19L))
  expect_true(all(mei_res_males[, 2L] == 0L)) # max L0 == 0
  
  expect_equal(length(match_females), length(indv_females))
  expect_equal(length(match_males), length(indv_males))
  
  expect_equal(4L, meiotic_dist(suspect, get_individual(test_pop, pid = 3L))) 
  expect_equal(meiotic_dist(suspect, get_individual(test_pop, pid = 3L)), 
               meioses_males[mei_res_males[, 3L] == 3L])
  
  expect_equal(4L, meiotic_dist(suspect, get_individual(test_pop, pid = 10L))) 
  expect_equal(meiotic_dist(suspect, get_individual(test_pop, pid = 10L)), 
               meioses_females[mei_res_females[, 3L] == 10L])
})


haps_from_ped <- get_haplotypes_in_pedigree(ped)
haps_from_pids <- get_haplotypes_pids(test_pop, pids)
haps_from_indvs <- get_haplotypes_individuals(indvs)
hap_from_indv <- lapply(pids, function(pid) get_haplotype(get_individual(test_pop, pid)))

test_that("pedigrees_all_populate_haplotypes haplotypes works", {
  #haps_from_ped
  expect_true(is.list(haps_from_ped))
  expect_equal(length(haps_from_ped), 12L)
  expect_equal(length(haps_from_ped[[1L]]), LOCI)
  expect_true(all(unlist(haps_from_ped) == 0L))
  
  #haps_from_pids
  expect_true(is.list(haps_from_pids))
  expect_equal(length(haps_from_pids), 12L)
  expect_equal(length(haps_from_pids[[1L]]), LOCI)
  expect_true(all(unlist(haps_from_pids) == 0L))
  
  #haps_from_indvs
  expect_equal(length(haps_from_indvs), 19L)
  expect_true(all(unique(unlist(haps_from_indvs))) == FALSE)
  expect_equal(as.integer(unique(unlist(haps_from_indvs))), 0L)
  
  #hap_from_indv
  expect_equal(haps_from_ped, hap_from_indv)
})

no_vars <- unlist(lapply(indvs, get_haplotype_no_variants))
test_that("get_haplotype_no_variants works", {
  expect_true(all(no_vars == 0L))
  expect_equal(no_vars, unlist(lapply(haps_from_indvs, sum)))
})

f_hap <- c(TRUE, rep(FALSE, LOCI-1L))
pedigrees_all_populate_haplotypes_custom_founders(peds, 
                                                  mutation_rates = rep(0, LOCI), 
                                                  get_founder_haplotype = function() f_hap,
                                                  progress = FALSE)

haps_from_ped <- get_haplotypes_in_pedigree(ped)
haps_from_pids <- get_haplotypes_pids(test_pop, pids)
haps_from_indvs <- get_haplotypes_individuals(indvs)
hap_from_indv <- lapply(pids, function(pid) get_haplotype(get_individual(test_pop, pid)))

test_that("pedigrees_all_populate_haplotypes_custom_founders works", {
  #haps_from_ped
  expect_true(is.list(haps_from_ped))
  expect_equal(length(haps_from_ped), 12L)
  expect_equal(length(haps_from_ped[[1L]]), LOCI)
  expect_true(length(unique(unlist(haps_from_ped))) == 2L)
  
  expect_equal(haps_from_ped, haps_from_pids)
  expect_equal(haps_from_ped, hap_from_indv)
  expect_equal(haps_from_indvs, lapply(seq_along(indvs), function(j) f_hap))
})

no_vars <- unlist(lapply(indvs, get_haplotype_no_variants))
test_that("get_haplotype_no_variants works", {
  expect_true(all(no_vars == 1L))
  expect_true(all(no_vars == sum(f_hap)))
  expect_equal(no_vars, unlist(lapply(haps_from_indvs, sum)))
})

##########################################################

hashmap <- build_haplotypes_hashmap(indvs)

test_that("hashmap works", {
  expect_output(print_haplotypes_hashmap(hashmap), "Total = 19", fixed = TRUE)
  expect_equal(length(indvs), length(get_haplotype_matching_individuals_from_hashmap(hashmap, f_hap)))
  expect_equal(0L, length(get_haplotype_matching_individuals_from_hashmap(hashmap, !f_hap)))
})


#######################################################

test_that("get_generations_from_final", {
  expect_equal(0L, get_generations_from_final(get_individual(test_pop, 1)))
  expect_equal(1L, get_generations_from_final(get_individual(test_pop, 6)))
})

test_that("print_individual", {
  expect_output(print_individual(indvs[[1L]]), 
                "  pid = 19 [M] in generation 1 with mother pid = 10 and no children",
                fixed = TRUE)
})

test_that("get_pedigree_id_from_pid", {
  expect_equal(get_pedigree_id(peds[[1]]), 1)
  expect_equal(get_pedigree_id(peds[[2]]), 2)
})

test_that("get_pedigree_id_from_pid", {
  expect_equal(get_pedigree_id_from_pid(test_pop, 1), get_pedigree_id_from_pid(test_pop, 6))
  expect_equal(get_pedigree_id_from_pid(test_pop, 15), get_pedigree_id_from_pid(test_pop, 16))
})

test_that("pedigrees_table", {
  expect_equal(pedigrees_table(peds), c("7" = 1, "12" = 1))
})

test_that("print_pedigree", {
  expect_output(print_pedigree(peds[[2]]), 
                "Pedigree with 7 individuals:
  18 [M] with mother 13
  13 [F] with mother 12
  12 [F] with mother -1
  14 [F] with mother 13
  15 [M] with mother 14
  16 [M] with mother 14
  17 [M] with mother 14",
                fixed = TRUE)
})

test_that("get_is_female_in_pedigree", {
  expect_equal(sum(get_is_female_in_pedigree(ped = peds[[1]])), 7)
  expect_equal(sum(!get_is_female_in_pedigree(ped = peds[[1]])), 5)
  
  expect_equal(sum(get_is_female_in_pedigree(ped = peds[[2]])), 3)
  expect_equal(sum(!get_is_female_in_pedigree(ped = peds[[2]])), 4)
})

peds_tidy <- get_pedigrees_tidy(pedigrees = peds)
test_that("get_pedigrees_tidy", {
  expect_equal(unlist(peds_tidy$ped_ids), c(1, 2))
})


#######################################################
# destructor
rm(test_pop)
gc()
#######################################################

