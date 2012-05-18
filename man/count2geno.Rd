\name{count2geno}
\alias{count2geno}
\title{Transforming genotype alleles counts to genotype matrix}
\description{Genotype data format transformation: from alleles counts to genotype coding matrix.}
\usage{count2geno(cgeno, indid)}
\arguments{
  \item{cgeno}{ Genotype counts data format: a matrix with 3 columns: indid (individual ID), snpid (SNP ID), and count}
  \item{indid}{ Individual ID, must include indid in cgeno. If missing, it is the same with indid in cgeno.}
}
\value{A matrix of genotypes}
\author{C. Xu}
\seealso{
\code{\link{geno2count}}
}
\examples{
data(sample.cgeno)
str(sample.cgeno) 
#count2geno(sample.cgeno$cgeno)
 }
\keyword{datasets}