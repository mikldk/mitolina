#' Mitogenome annotated
#'
#' A dataset containing information about the mitogenome.
#' 
#' Note, that `PartitionRieux` and `PartitionOversti` have been 
#' modified to obtain the same number of positions as were given in the papers.
#'
#' @format A tibble / data frame with 16,569 rows and 5 variables:
#' \describe{
#'   \item{Position}{Position in the mitogenome}
#'   \item{rCRS}{The base in the revised Cambridge Reference Sequence (rCRS)}
#'   \item{PartitionPhyloTree}{Partition according to \url{http://www.phylotree.org}}
#'   \item{PartitionRieux}{Partition according to Rieux et al. (2014) (see references)}
#'   \item{PartitionOversti}{Partition according to Rieux et al. (2014) and Översti et al. 2017 (see references)}
#' }
#' @source \url{http://www.phylotree.org/resources/rCRS_annotated.htm}
#' @seealso mtdna_mut_Rieux mtdna_mut_Oversti
#' @references 
#' * PhyloTree: \url{http://www.phylotree.org/resources/rCRS_annotated.htm}; 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_partitions"

#' Mitogenome mutation rates according to Rieux et al. (2014)
#'
#' A dataset containing mutation information about the mitogenome 
#' partitioned into 4 partitions 
#' according to Rieux et al. (2014).
#' 
#' Note the value `[Not included?]` in the `PartitionRieux` column.
#'
#' @format A tibble / data frame with 5 rows and 6 variables:
#' \describe{
#'   \item{PartitionRieux}{Partition according to Rieux et al. (2014) (see references)}
#'   \item{MutRateYearly}{Mutation rate per site per year}
#'   \item{MutHPDLower}{Mutation rate highest posterior density lower limit}
#'   \item{MutHPDUpper}{Mutation rate highest posterior density upper limit}
#'   \item{MutRateSEYearly}{Standard error of `MutRateYearly` by assuming that the highest posterior region is a 95\% normal}
#'   \item{MutationScheme}{Name of the mutation scheme}
#' }
#' @source \url{https://doi.org/10.1093/molbev/msu222}
#' @seealso mtdna_partitions mtdna_mut_Oversti
#' @references 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
"mtdna_mut_Rieux"

#' Mitogenome mutation rates according to Översti et al. (2017)
#'
#' A dataset containing mutation information about the mitogenome 
#' partitioned into 4 partitions 
#' according to Rieux et al. (2014).
#'
#' Note the value \code{[Not included?]} in the \code{PartitionOversti} column.
#'
#' @format A tibble / data frame with 5 rows and 6 variables:
#' \describe{
#'   \item{PartitionOversti}{Partition according to Översti et al. (2017) (see references)}
#'   \item{MutRateYearly}{Mutation rate per site per year}
#'   \item{MutHPDLower}{Mutation rate highest posterior density lower limit}
#'   \item{MutHPDUpper}{Mutation rate highest posterior density upper limit}
#'   \item{MutRateSEYearly}{Standard error of `MutRateYearly` by assuming that the highest posterior region is a 95\% normal}
#'   \item{MutationScheme}{Name of the mutation scheme}
#' }
#' @source \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}
#' @seealso mtdna_partitions mtdna_mut_Rieux
#' @references 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_mut_Oversti"

#' 588 forensic-quality haplotypes representing three U.S. populations
#'
#' Full mtGenome reference data: Development and characterization of 588 forensic-quality 
#' haplotypes representing three U.S. populations
#'
#' @format A tibble / data frame with 588 rows and 5 variables:
#' \describe{
#'   \item{id}{Sample id}
#'   \item{HG}{Haplogroup}
#'   \item{Variants}{List column with variants for each individual}
#'   \item{Origin}{Origin}
#'   \item{State}{US state}
#' }
#' @source \url{https://doi.org/10.1016/j.fsigen.2014.09.021}
#' @references Just et al. "Full mtGenome reference data: Development and characterization of 588 
#' forensic-quality haplotypes representing three U.S. populations", Forensic Science International: 
#' Genetics, 14, 2015.
"mtdna_vars_Just"
