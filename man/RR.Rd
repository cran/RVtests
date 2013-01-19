\name{RR}
\alias{RR}
\title{Ridge Regression for RV Tests}
\description{Use ridge regression for testing rare variants association with disease traits.}
\usage{
RR(x, y, z = NULL, scale = FALSE, weights = 1, lambda = 1, npermutation = 1000, npermutation.max, min.nonsignificant.counts = 0.1*npermutation)
}
\arguments{
  \item{x}{ Genotype matrix}
  \item{y}{ Phenotype vector}
  \item{z}{ Covariate matrix}	
  \item{scale}{ If TRUE, scale x and y.}
  \item{weights}{ Genotype weights}
	
  \item{lambda}{ Regularization parameter}
  \item{npermutation}{ Number of permutation}
  \item{npermutation.max}{ Maximum permutation}
  \item{min.nonsignificant.counts}{ Minimum nonsignificant counts}
}
\value{
\item{nonsignificant.counts }{Counts of permuted data that have a higher score than unpermuted data.}
\item{total.permutation }{Total permutation}
\item{score }{Correlation between y and y_est if z=NULL. }
\item{pvalue.empirical }{Empirical pvalue via permutation}
\item{pvalue.nominal }{Theoretical pvalue, not available.}
}
\author{C. Xu}
\seealso{
\code{\link{PCR}, \link{PLS}}
}
\keyword{ models }
\keyword{ multivariate }