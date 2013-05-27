\name{RVtests-package}
\alias{RVtests-package}
\alias{RVtests}
\docType{package}
\title{Rare Variants Tests}
\description{
Use multiple regression methods to test rare variants association with disease traits.
}
\details{
\tabular{ll}{
Package: \tab RVtests\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2013-05-27\cr
License: \tab GLP 2.0 or greater\cr
LazyLoad: \tab yes\cr
}
An overview of how to use the package, including the most important functions 
}
\author{
C. Xu,
C. M. Greenwood,
Maintainer: <changjiang.xu@mail.mcgill.ca>
}
\references{
Xu C, Ladouceur M, Dastani Z, Richards JB, Ciampi A, Greenwood CMT. (2012) Multiple Regression Methods Show Great Potential for Rare Variant Association Tests. PLoS ONE 7(8): e41694. doi:10.1371/journal.pone.0041694
}
\examples{
data(sample.cgeno)
str(sample.cgeno) 
x=count2geno(sample.cgeno$cgeno)
dim(x)

set.seed(31018)
y<- rowSums(x[,2:4]*rep(rnorm(3,1,0.1), each=nrow(x))) + 0.4*rnorm(nrow(x))

tmp<- proc.time();RR(x,y,lambda=0:5); proc.time()-tmp
tmp<- proc.time();RR(x,y,weights=c(rep(2,10), rep(1, ncol(x)-10)), lambda=0:5); proc.time()-tmp
tmp<- proc.time();RR(x,y,weights=c(rep(1,10), rep(0, ncol(x)-10)), lambda=0:5); proc.time()-tmp
}
