## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(fig.width = 7)

## ---- message=FALSE------------------------------------------------------
library(mitolina)

## ------------------------------------------------------------------------
set.seed(1)

## ------------------------------------------------------------------------
data(mtdna_partitions)
mtdna_partitions %>% print(n = 2)

## ------------------------------------------------------------------------
data(mtdna_mut_SoaresK8)
mtdna_mut_SoaresK8

## ------------------------------------------------------------------------
data(mtdna_mut_RieuxK4)
mtdna_mut_RieuxK4

data(mtdna_mut_RieuxK1)
mtdna_mut_RieuxK1

## ------------------------------------------------------------------------
data(mtdna_mut_OverstiK4)
mtdna_mut_OverstiK4

## ------------------------------------------------------------------------
data(mtdna_mut_schemes)
mtdna_mut_schemes %>% print(n = 2)
mtdna_mut_schemes %>% count(MutationScheme)

## ------------------------------------------------------------------------
d_mu <- mtdna_mut_schemes %>% 
  filter(MutationScheme == "Rieux 2014 (K = 4)", MutRateYearly > 0)
mu_mean_year <- d_mu %>% pull(MutRateYearly)
mu_mean_gen <- 25*mu_mean_year

## ------------------------------------------------------------------------
1 - prod(1 - mu_mean_gen)

