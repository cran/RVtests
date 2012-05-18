\name{RVtests-package}
\alias{RVtests-package}
\alias{RVtests}
\docType{package}
\title{Rare Variant Tests Using Multiple Regression Methods}
\description{
Use regularized multiple regression methods to test rare variants association with disease traits.
Include 5 methods: ridge regression (RR), principal components regression (PCR), partial least 
squares regression (PLS), sparse PLS (SPLS), and the LASSO.
The empirical pvalue is obtained by an adaptive permutation. 
The permutation parameters should be properly larger to get a reliable pvalue.
Usually let the minimum permutation, npermutation = 1000,
and the minimum nonsignificant counts, min.nonsignificant.counts = 50.
If the maximum permutation equals the minimum permutation, a full permutation is run.
}
\details{
\tabular{ll}{
Package: \tab RVtests\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-05-18\cr
License: \tab GLP 2.0 or greater\cr
LazyLoad: \tab yes\cr
}
}
\author{
C. Xu,
C. M. Greenwood.
Maintainer: C. Xu <changjiang.xu@mail.mcgill.ca>
}
\references{
C. Xu, M. Ladouceur, Z. Dastani, J. B. Richards, A. Ciampi, C. M.T. Greenwood (2012), Multiple regression methods show great potential for rare variant association tests, PLoSONE.
}
\keyword{ package }
\seealso{
\code{\link[SKAT]{SKAT}, \link[thgenetics]{rareGeneTest}} 
}