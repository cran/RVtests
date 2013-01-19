\name{sample.cgeno}
\alias{sample.cgeno}
\docType{data}
\title{Genotype Counts dataset}
\description{A list of genotype counts, phenotype, and polyphen weight}
\usage{data(sample.cgeno)}
\format{
  The format is:
List of 3
 $ cgeno          :'data.frame':	960 obs. of  3 variables:
  ..$ indid: int [1:960] 16929 18167 18168 28671 31308 32037 33182 49716 53138 13206 ...
  ..$ snpid: int [1:960] 57548466 57548466 57548466 57548466 57548466 57548466 57548466 57548466 57548466 57550649 ...
  ..$ count: int [1:960] 1 1 1 1 1 1 1 1 1 1 ...
 $ phen           :'data.frame':	262 obs. of  2 variables:
  ..$ indid: int [1:262] 32 90 101 109 129 133 225 236 253 349 ...
  ..$ trait: num [1:262] 0.128 0.166 0.884 0.929 0.195 ...
 $ polyphen.weight:'data.frame':	71 obs. of  2 variables:
  ..$ snpid : int [1:71] 57548364 57548466 57550649 57550666 57556205 57556220 57556236 57567762 57569339 57569466 ...
  ..$ weight: num [1:71] 0.5 0.5 0.5 0.5 0.055 0 0.706 0.5 0.995 0.5 ...
}
\details{The dataset was used in comparing VT and WOD methods.}
\examples{
data(sample.cgeno)
str(sample.cgeno) 
}
\keyword{datasets}
