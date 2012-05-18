\name{geno2count}
\alias{geno2count}
\title{Transforming genotype codes matrix to counts}
\description{Genotype data format transformation: from genotype codes matrix to counts. The '0' codes will be removed.}
\usage{geno2count(genotype)}
\arguments{
  \item{genotype}{ Genotype codes matrix}
}
\value{
Data frame of genotype counts with 3 columns: indid, snpid, and count
}
\author{C. Xu}
\seealso{
\code{\link{count2geno}}
}
\keyword{datasets}