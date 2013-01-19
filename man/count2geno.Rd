\name{count2geno}
\alias{count2geno}
\title{Transforming genotype counts to genotype codes}
\description{
Transform genotype counts data format to genotype codes format.
}
\usage{
count2geno(cgeno, indid)
}
\arguments{
  \item{cgeno}{ A matrix with 3 columns: indid, snpid, and count}
  \item{indid}{ Individual ID, including indid in cgeno}
}
\value{
  A matrix of genotypes 
}
\author{C. Xu}

\seealso{
\code{\link{geno2count}}
}
\keyword{datasets}