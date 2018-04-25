#' Mitogenome annotated
#'
#' A dataset containing information about the mitogenome.
#'
#' @format A tibble / data frame with 16,569 rows and 5 variables:
#' \describe{
#'   \item{Position}{Position in the mitogenome}
#'   \item{rCRS}{The base in the revised Cambridge Reference Sequence (rCRS)}
#'   \item{PartitionPhyloTree}{Partition according to \url{http://www.phylotree.org}}
#'   \item{PartitionRieux}{Partition according to Rieux et al. (2014) (see references)}
#'   \item{PartitionSoares}{Partition according to Soares et al. (2009) (see references)}
#' }
#' @source \url{http://www.phylotree.org/resources/rCRS_annotated.htm}
#' @seealso mtdna_mut_schemes mtdna_mut_SoaresK8 mtdna_mut_RieuxK4 mtdna_mut_RieuxK1 mtdna_mut_OverstiK4
#' @references 
#' * PhyloTree: \url{http://www.phylotree.org/resources/rCRS_annotated.htm}; 
#' * Soares et al. (2009): \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}; 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_partitions"

#' Mitogenome mutation rates according to Soares et al. (2009)
#'
#' A dataset containing mutation information about the mitogenome according to
#' Soares et al. (2009).
#' Soares et al. (2009) does not give standard errors par partition, but they
#' write "The substitution rate for the entire molecule was 
#' \eqn{1.665 \times 10^{-8}} (\eqn{\pm 1.479\times 10^{-9}}).". 
#' We assume that the error margin is a 95\% confidence interval for a
#' normal distribution. This gives a standard error of
#' \eqn{1.479\times 10^{-9}/1.96 = 7.5459184\times 10^{-10}}. 
#'
#' @format A tibble / data frame with 8 rows and 4 variables:
#' \describe{
#'   \item{PartitionSoares}{Partition according to Soares et al. (2009) (see references)}
#'   \item{MutRateYearly}{Mutation rate per site per year}
#'   \item{MutRateSEYearly}{Standard error of `MutRateYearly`}
#'   \item{MutationScheme}{Name of the mutation scheme}
#' }
#' @source \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}
#' @seealso mtdna_mut_schemes mtdna_partitions mtdna_mut_RieuxK4 mtdna_mut_RieuxK1 mtdna_mut_OverstiK4
#' @references 
#' * PhyloTree: \url{http://www.phylotree.org/resources/rCRS_annotated.htm}; 
#' * Soares et al. (2009): \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}; 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_mut_SoaresK8"

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
#' @seealso mtdna_mut_schemes mtdna_partitions mtdna_mut_SoaresK8 mtdna_mut_RieuxK1 mtdna_mut_OverstiK4
#' @references 
#' * PhyloTree: \url{http://www.phylotree.org/resources/rCRS_annotated.htm}; 
#' * Soares et al. (2009): \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}; 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_mut_RieuxK4"

#' Mitogenome mutation rates according to Rieux et al. (2014)
#'
#' A dataset containing mutation information about the mitogenome 
#' partitioned into 1 partitions 
#' according to Rieux et al. (2014).
#' 
#' @format A tibble / data frame with 1 row and 5 variables:
#' \describe{
#'   \item{MutRateYearly}{Mutation rate per site per year}
#'   \item{MutHPDLower}{Mutation rate highest posterior density lower limit}
#'   \item{MutHPDUpper}{Mutation rate highest posterior density upper limit}
#'   \item{MutRateSEYearly}{Standard error of `MutRateYearly` by assuming that the highest posterior region is a 95\% normal}
#'   \item{MutationScheme}{Name of the mutation scheme}
#' }
#' @source \url{https://doi.org/10.1093/molbev/msu222}
#' @seealso mtdna_mut_schemes mtdna_partitions mtdna_mut_RieuxK4 mtdna_mut_SoaresK8 mtdna_mut_OverstiK4
#' @references 
#' * PhyloTree: \url{http://www.phylotree.org/resources/rCRS_annotated.htm}; 
#' * Soares et al. (2009): \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}; 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_mut_RieuxK1"

#' Mitogenome mutation rates according to Översti et al. (2017)
#'
#' A dataset containing mutation information about the mitogenome 
#' partitioned into 4 partitions 
#' according to Rieux et al. (2014).
#'
#' Note the value \code{[Not included?]} in the \code{PartitionRieux} column.
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
#' @source \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}
#' @seealso mtdna_mut_schemes mtdna_partitions mtdna_mut_RieuxK4 mtdna_mut_RieuxK1 mtdna_mut_SoaresK8
#' @references 
#' * PhyloTree: \url{http://www.phylotree.org/resources/rCRS_annotated.htm}; 
#' * Soares et al. (2009): \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}; 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_mut_OverstiK4"

#' Mitogenome with mutation information
#'
#' A dataset containing information about the mitogenome.
#'
#' @format A tibble / data frame with 66,276 (4*16,569) rows and 4 variables:
#' \describe{
#'   \item{Position}{Position in the mitogenome}
#'   \item{MutationScheme}{Mutation scheme: 
#'   "Soares 2009 (K = 8)" [mtdna_mut_SoaresK8], 
#'   "Rieux 2014 (K = 1)" [mtdna_mut_RieuxK1], 
#'   "Rieux 2014 (K = 4)" [mtdna_mut_RieuxK4], or
#'   "Översti 2017 (K = 4) 16569" [mtdna_mut_OverstiK4]}.
#'   \item{MutRateYearly}{Mutation rate per site per year}
#'   \item{MutRateSEYearly}{Standard error of `MutRateYearly`; please refer to 
#'      [mtdna_mut_SoaresK8], [mtdna_mut_RieuxK1], 
#'      [mtdna_mut_RieuxK4], or [mtdna_mut_OverstiK4] for details.}
#' }
#' @source \url{http://www.phylotree.org/resources/rCRS_annotated.htm}
#' @seealso mtdna_partitions mtdna_mut_SoaresK8 mtdna_mut_RieuxK4 mtdna_mut_RieuxK1 mtdna_mut_OverstiK4
#' @references 
#' * PhyloTree: \url{http://www.phylotree.org/resources/rCRS_annotated.htm}; 
#' * Soares et al. (2009): \url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}; 
#' * Rieux et al. (2014): \url{https://doi.org/10.1093/molbev/msu222}; 
#' * Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
"mtdna_mut_schemes"
