\name{SPLS}
\alias{SPLS}
\title{Sparse PLS for RV Tests}
\description{Use SPLS for selecting significant variants and testing the variants associated with disease traits.}
\usage{
SPLS(x, y, scale = TRUE, ncomp, eta.grid, size.max, a = 2, 
	npermutation = 0, npermutation.max, min.nonsignificant.counts)
}
\arguments{
 \item{x}{ Genotype matrix, each row as an individual and each column as a snp}
  \item{y}{ Phenotype vector}
  \item{scale}{ see spls}
  \item{ncomp}{ Number of components}
  \item{eta.grid}{ see spls}
  \item{size.max}{ Maximum number of variants included}
  \item{a}{ Penalty parameter for information criterion, a=2 for AIC.}
  \item{npermutation}{ Number of permutation, if less than 1, the permutation will not be run.}
  \item{npermutation.max}{ Maximum permutation}
  \item{min.nonsignificant.counts}{ Minimum nonsignificant counts}
}
\details{
Use spls package to implement SPLS and an information criterion (AIC, BIC, GIC) to select a set of variants.
}
\value{
\item{nonsignificant.counts }{Counts of permuted data that have a higher score than unpermuted data.}
\item{pvalue.empirical }{Empirical pvalue via permutation}
\item{pvalue.nominal }{Not availabe}
\item{vs }{The selected variants}
\item{total.permutation }{Total permutation}
}
\references{
Xu C, Ladouceur M, Dastani Z, Richards JB, Ciampi A, Greenwood CMT. (2012) Multiple Regression Methods Show Great Potential for Rare Variant Association Tests. PLoS ONE 7(8): e41694. doi:10.1371/journal.pone.0041694
}
\author{C. Xu}
\seealso{
\code{\link{spls}, \link{LASSO}}
}
\keyword{ models }
\keyword{ multivariate }