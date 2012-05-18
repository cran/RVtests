\name{PCR}
\alias{PCR}
\title{Principal Components Regression for RV tests}
\description{Use principal components for testing rare variants association with disease traits.}
\usage{
PCR(x, y, scale = FALSE, ncomp, varpercent, npermutation = 100, npermutation.max, min.nonsignificant.counts)
}
\arguments{
  \item{x}{ Genotype matrix}
  \item{y}{ Phenotype vector}
  \item{scale}{ If TRUE, scale x and y.}
  \item{ncomp}{ Number of components, which could be a vector containing a set of numbers.}
  \item{varpercent}{ Explained variance percentage}
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
\item{ncomp.varp }{Number of components required for specified variance percentage}
}
\references{
C. Xu, M. Ladouceur, Z. Dastani, J. B. Richards, A. Ciampi, C. M.T. Greenwood (2012), Multiple regression methods show great potential for rare variant association tests, PLoSONE.
}
\author{C. Xu}
\seealso{
\code{\link{PLS}, \link{RR}}
}
\keyword{ models }
\keyword{ multivariate }
