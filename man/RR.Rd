\name{RR}
\alias{RR}
\title{Ridge Regression for RV Tests}
\description{Use ridge regression for testing rare variants association with disease traits.}
\usage{
RR(x, y, scale = FALSE, lambda = 1, npermutation = 100, npermutation.max, min.nonsignificant.counts)
}
\arguments{
  \item{x}{ Genotype matrix}
  \item{y}{ Phenotype vector}
  \item{scale}{ If TRUE, scale x and y.}
  \item{lambda}{ Regularization parameter}
  \item{npermutation}{ Number of permutation, if less than 1, the permutation will not be run.}
  \item{npermutation.max}{ Maximum permutation, if missing, equal to npermutation.}
  \item{min.nonsignificant.counts}{ Minimum nonsignificant counts, if missing, equal to 10.}
}
\value{
\item{score }{Correlation between y and y_est}
\item{nonsignificant.counts }{Counts of permuted data that have a higher score than unpermuted data.}
\item{pvalue.empirical }{Empirical pvalue via permutation}
\item{pvalue.nominal }{Theoretical pvalue, not available now.}
\item{total.permutation }{Total permutation}
}
\references{
C. Xu, M. Ladouceur, Z. Dastani, J. B. Richards, A. Ciampi, C. M.T. Greenwood (2012), Multiple regression methods show great potential for rare variant association tests, PLoSONE.
}
\author{C. Xu}
\seealso{
\code{\link{PCR}, \link{PLS}}
}
\keyword{ models }
\keyword{ multivariate }