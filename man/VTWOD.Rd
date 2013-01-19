\name{VTWOD}
\alias{VTWOD}
\title{VT and WOD for RV Tests}
\description{Include methods: T1, T5, WE, VT, and WOD.}
\usage{
VTWOD(x, y, polyphen.weight, flipPhenotype = 0, npermutation = 1000, npermutation.max, min.nonsignificant.counts)
}
\arguments{
 \item{x}{ Genotype matrix}
  \item{y}{ Phenotype vector}
  \item{polyphen.weight}{ Polyphen weight} 
  \item{flipPhenotype}{ Logical, if TRUE, flip phenotype to opposite by multipling -1} 
  \item{npermutation}{ Number of permutation, if less than 1, the permutation will not be run.}
  \item{npermutation.max}{ Maximum permutation}
  \item{min.nonsignificant.counts}{ Minimum nonsignificant counts}
}
\value{
\item{score }{Scores of T1, T5, WE, VT, and WOD}
\item{nonsignificant.counts }{Counts of permuted data that have a higher score than unpermuted data.}
\item{pvalue.empirical }{Empirical pvalue via permutation}
\item{pvalue.nominal }{Theoretical pvalue, not available now.}
\item{total.permutation }{Total permutation}
}
\references{
Ladouceur M, Dastani Z, Aulchenko YS, Greenwood CM, Richards JB (2012) The empirical power of rare variant association methods: Results from sanger sequencing in 1,998 individuals. PloS Genetics 8: e1002496.

Price AL, Kryukov GV, de Bakker PI, Purcell SM, Staples J, et al. (2010) Pooled association tests for rare variants in exon-resequencing studies. Am J Hum Genet 86: 832 - 838.
}

\author{ C. Xu}

\note{
This R implementation by Adam Kiezun, based on reference implementation in C by Alkes Price.
Added WOD tests to the program in 2011 by Celia Greenwood
}
\seealso{
\code{\link{RR}, \link{PCR}, \link{PLS}}
}
\keyword{ models }
\keyword{ multivariate }