\name{PLS}
\alias{PLS}
\title{Partial Least Squares Regression for RV tests}
\description{Use PLS components for testing rare variants association with disease traits.}
\usage{
PLS(x, y, scale = FALSE, ncomp, varpercent, npermutation = 100, npermutation.max, min.nonsignificant.counts)
}
\arguments{
  \item{x}{ Genotype matrix}
  \item{y}{ Phenotype vector}
  \item{scale}{ If TRUE, scale x and y.}
  \item{ncomp}{ Number of components, which could be a vector containing a set of numbers.}
  \item{varpercent}{ Explained variance percentage}
  \item{npermutation}{ Number of permutation, if less than 1, the permutation will not be run.}
  \item{npermutation.max}{ Maximum permutation}
  \item{min.nonsignificant.counts}{ Minimum nonsignificant counts}
}
\value{
\item{score }{Correlation between y and y_est}
\item{nonsignificant.counts }{Counts of permuted data that have a higher score than unpermuted data.}
\item{pvalue.empirical }{Empirical pvalue via permutation}
\item{pvalue.nominal }{Theoretical pvalue, not available now.}
\item{total.permutation }{Total permutation}
\item{ncomp.varp }{Number of components required for specified variance percentage}
}
\author{C. Xu}
\seealso{
\code{\link{PCR}, \link{SPLS}}
}
\keyword{ models }
\keyword{ multivariate }