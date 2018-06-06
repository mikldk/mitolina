## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(fig.width = 7, cache = FALSE)

## ---- message=FALSE------------------------------------------------------
library(mitolina)

## ------------------------------------------------------------------------
set.seed(1)

## ---- message=FALSE------------------------------------------------------
library(dplyr)

## ------------------------------------------------------------------------
sim <- sample_mtdna_geneology_varying_size(
  population_sizes_females = rep(10, 10),
  population_sizes_males = c(0L, rep(10, 9)),
  progress = FALSE)
pop <- sim$population
pop

## ------------------------------------------------------------------------
peds <- build_pedigrees(pop, progress = FALSE)
peds

## ------------------------------------------------------------------------
end_gen_peds <- lapply(sim$end_generation_male_individuals, get_pedigree_from_individual)
pedids <- sapply(end_gen_peds, get_pedigree_id)
id_tab <- sort(table(pedids))
ped_id_endgen <- as.integer(names(id_tab)[1]) # smallest
indv_index <- which(pedids == ped_id_endgen)[1]
ped <- get_pedigree_from_individual(sim$end_generation_male_individuals[[indv_index]])
pedigree_size(ped)

## ------------------------------------------------------------------------
plot(ped)

## ------------------------------------------------------------------------
sites <- 100

## ------------------------------------------------------------------------
get_founder_mito <- function() {
  sample(c(FALSE, TRUE), sites, replace = TRUE)
}
x <- get_founder_mito()
head(x)
sum(x)

## ------------------------------------------------------------------------
mu <- rep(1e-2, sites)

## ------------------------------------------------------------------------
pedigrees_all_populate_haplotypes_custom_founders(
  pedigrees = peds,
  mutation_rates = mu,
  get_founder_haplotype = get_founder_mito, 
  progress = FALSE)

## ------------------------------------------------------------------------
plot(ped, num_vars = TRUE)

## ---- fig.height=8, fig.width=12-----------------------------------------
library(tidygraph)
library(ggraph)

g <- as_tbl_graph(peds) %>% 
  activate(nodes) %>% 
  filter(ped_id == ped_id_endgen) %>% 
  mutate(num_vars = sapply(haplotype, sum))

p <- ggraph(g, layout = 'tree') +
  geom_edge_link() +
  geom_node_point(aes(color = sex, fill = num_vars, shape = sex), size = 12, stroke = 2) +
  geom_node_text(aes(label = paste0(name, "\n", num_vars)), size = 4, color = "white") +
  theme_graph(base_family = "") +
  scale_shape_manual(NULL, values = c("Female" = 21, "Male" = 22)) +
  scale_color_manual(NULL, values = c("Female" = "black", "Male" = "red"))
p

## ------------------------------------------------------------------------
haps <- get_haplotypes_in_pedigree(ped)
length(haps)
mean(sapply(haps, sum))

## ------------------------------------------------------------------------
haps <- get_haplotypes_individuals(sim$end_generation_male_individuals)
length(haps)

## ------------------------------------------------------------------------
haps <- get_haplotypes_individuals(sim$male_individuals_generations)
length(haps)

## ------------------------------------------------------------------------
data(mtdna_partitions)
mtdna_partitions %>% print(n = 2)

## ------------------------------------------------------------------------
data(mtdna_mut_Rieux)
mtdna_mut_Rieux

## ------------------------------------------------------------------------
data(mtdna_mut_Oversti)
mtdna_mut_Oversti

## ------------------------------------------------------------------------
mu_mean_year <- mtdna_partitions %>% 
  left_join(mtdna_mut_Rieux, by = "PartitionRieux") %>% pull(MutRateYearly)
mu_mean_gen <- 25*mu_mean_year

## ------------------------------------------------------------------------
1 - prod(1 - mu_mean_gen)

## ------------------------------------------------------------------------
d_mu_sample <- mtdna_mut_Rieux %>% 
  rowwise() %>% 
  mutate(mu = 25*rnorm(1, mean = MutRateYearly, sd = MutRateSEYearly))

mu_sample <- mtdna_partitions %>% 
  left_join(d_mu_sample, by = "PartitionRieux") %>% pull(mu)

## ------------------------------------------------------------------------
get_overall_mu <- function() {
  d_mu_sample <- mtdna_mut_Rieux %>% 
    rowwise() %>% 
    mutate(mu = 25*rnorm(1, mean = MutRateYearly, sd = MutRateSEYearly))
  
  mu_sample <- mtdna_partitions %>% 
    left_join(d_mu_sample, by = "PartitionRieux") %>% pull(mu)
  
  1 - prod(1 - mu_sample)
}

set.seed(10)
mu_overall <- replicate(100, get_overall_mu())

summary(mu_overall)
mean(mu_overall)
qnorm(0.975)*sd(mu_overall)
quantile(mu_overall, c(0.025, 0.975))

## ------------------------------------------------------------------------
get_founder_mito <- function() {
  sample(c(FALSE, TRUE), length(mu_mean_gen), replace = TRUE)
}

pedigrees_all_populate_haplotypes_custom_founders(
  pedigrees = peds,
  mutation_rates = mu_mean_gen,
  get_founder_haplotype = get_founder_mito, 
  progress = FALSE)

## ------------------------------------------------------------------------
haps <- get_haplotypes_individuals(sim$male_individuals_generations)
table(sapply(haps, sum))

## ------------------------------------------------------------------------
haps <- get_haplotypes_individuals(sim$end_generation_male_individuals)
table(sapply(haps, sum))

## ---- fig.height=8, fig.width=12-----------------------------------------
g <- as_tbl_graph(peds) %>% 
  activate(nodes) %>% 
  filter(ped_id == ped_id_endgen) %>% 
  mutate(num_vars = sapply(haplotype, sum))

p <- ggraph(g, layout = 'tree') +
  geom_edge_link() +
  geom_node_point(aes(color = sex, fill = factor(num_vars), shape = sex), size = 12, stroke = 2, show.legend = FALSE) +
  geom_node_text(aes(label = paste0(name, "\n", num_vars)), size = 4, color = "white") +
  theme_graph(base_family = "") +
  scale_shape_manual(NULL, values = c("Female" = 21, "Male" = 22)) +
  scale_color_manual(NULL, values = c("Female" = "black", "Male" = "red"))
p

