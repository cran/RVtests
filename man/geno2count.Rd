\name{geno2count}
\alias{geno2count}
\title{Transforming genotype codes matrix to genotype counts}
\description{Genotype counts }
\usage{
geno2count(genotype)
}
\arguments{
  \item{genotype}{ Genotype matrix or data frame with row and column names, each row as an individual and each column as a snp}
}
\value{
Data frame of genotype counts with 3 columns: indid (individual IDs), snpid (SNP IDs), and count
}
\author{C. Xu}
\seealso{
\code{\link{count2geno}}
}
\keyword{datasets}