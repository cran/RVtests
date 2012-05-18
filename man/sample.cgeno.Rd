\name{sample.cgeno}
\alias{sample.cgeno}
\docType{data}
\title{Genotype counts dataset}
\description{A list of genotype counts, phenotype, and polyphen weight}
\usage{data(sample.cgeno)}
\format{
List of 3

 cgeno: data.frame with 960 observations of 3 variables, indid (individual ID), snpid (SNP ID), and count.

 phen: data.frame with 262 observations of  2 variables, indid (individual ID), and trait.

 polyphen.weight: data.frame with 71 obs. of  2 variables snpid (SNP ID), and weight.
}
\details{The dataset was used in comparing VT and WOD methods, see reference below.}
\references{
Ladouceur M, Dastani Z, Aulchenko YS, Greenwood CM, Richards JB (2012) 
The empirical power of rare variant association methods: Results from sanger sequencing in 1,998 individuals.
PloS Genetics 8: e1002496.
}
\examples{
data(sample.cgeno)
str(sample.cgeno) 
}
\keyword{datasets}
