\name{PCR}
\alias{PCR}
\title{Principal Components Regression for RV tests}
\description{Use principal components for testing rare variants association with disease traits.}
\usage{
PCR(x, y, scale = FALSE, ncomp, varpercent, 
	npermutation = 100, npermutation.max, min.nonsignificant.counts)
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
\references{
Xu C, Ladouceur M, Dastani Z, Richards JB, Ciampi A, Greenwood CMT. (2012) Multiple Regression Methods Show Great Potential for Rare Variant Association Tests. PLoS ONE 7(8): e41694. doi:10.1371/journal.pone.0041694
}
\author{C. Xu}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{PLS}, \link{RR}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }
