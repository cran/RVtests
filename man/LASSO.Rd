\name{LASSO}
\alias{LASSO}
\title{LASSO for Rare Variant Tests}
\description{
Use LASSO for selecting significant variants and testing the variants associated with disease traits.
}
\usage{
LASSO(x, y, family = c("gaussian", "binomial", "poisson", "multinomial", "cox"), 
	alpha = 1, nlambda = 100, lambda.min.ratio, standardize = TRUE, 
	size.max, a = 2, npermutation = 0, npermutation.max, min.nonsignificant.counts)
}
\arguments{
  \item{x}{ Genotype matrix}
  \item{y}{ Phenotype vector}
  \item{family}{ Family: gaussian, binomial, poisson, multinomial, and cox}
  \item{alpha}{ alpha = 1 for LASSO, see R package glmnet}
  \item{nlambda}{ see glmnet}
  \item{lambda.min.ratio}{ see glmnet}
  \item{standardize}{ see glmnet}
  \item{size.max}{ Maximum number of variants included}
  \item{a}{ Penalty parameter for information criterion, a=2 for AIC.}
  \item{npermutation}{ Number of permutation, if less than 1, the permutation will not be run.}
  \item{npermutation.max}{ Maximum permutation, if missing, equal to npermutation.}
  \item{min.nonsignificant.counts}{ Minimum nonsignificant counts, if missing, equal to 10.}
}
\details{
Use glmnet package to implement LASSO and an information criterion, AIC, BIC, or GIC, to select the best subset of variants.
}
\value{
\item{nonsignificant.counts }{Counts of permuted data that have a higher score than unpermuted data.}
\item{pvalue.empirical }{Empirical pvalue via permutation}
\item{pvalue.nominal }{Theoretical pvalue for the selected variants}
\item{vs }{The selected variants}
\item{total.permutation }{Total permutation}
\item{family }{Family}
}
\references{
C. Xu, M. Ladouceur, Z. Dastani, J. B. Richards, A. Ciampi, C. M.T. Greenwood (2012), Multiple regression methods show great potential for rare variant association tests, PLoSONE.
}
\author{C. Xu}
\seealso{
\code{\link{SPLS}, \link{glmnet}}
}
\keyword{ models }
\keyword{ multivariate }