% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mtdna_mut_Oversti}
\alias{mtdna_mut_Oversti}
\title{Mitogenome mutation rates according to Översti et al. (2017)}
\format{A tibble / data frame with 5 rows and 6 variables:
\describe{
\item{PartitionOversti}{Partition according to Översti et al. (2017) (see references)}
\item{MutRateYearly}{Mutation rate per site per year}
\item{MutHPDLower}{Mutation rate highest posterior density lower limit}
\item{MutHPDUpper}{Mutation rate highest posterior density upper limit}
\item{MutRateSEYearly}{Standard error of \code{MutRateYearly} by assuming that the highest posterior region is a 95\% normal}
\item{MutationScheme}{Name of the mutation scheme}
}}
\source{
\url{https://dx.doi.org/10.1016/j.ajhg.2009.05.001}
}
\usage{
mtdna_mut_Oversti
}
\description{
A dataset containing mutation information about the mitogenome
partitioned into 4 partitions
according to Rieux et al. (2014).
}
\details{
Note the value \code{[Not included?]} in the \code{PartitionOversti} column.
}
\references{
\itemize{
\item Översti et al. (2017): \url{https://doi.org/10.1038/s41598-017-05673-7}
}
}
\seealso{
mtdna_partitions mtdna_mut_Rieux
}
\keyword{datasets}
