\name{simGeno}
\alias{simGeno}
\title{
Function to simulate SNP matrix
}
\description{
Simulate SNP matrix coded 0, 1 and 2; with random genotypes.
}
\usage{
simGeno(Nind, Nmarkers)
}
\arguments{
  \item{Nind}{
number of individuals to simulate.
}
  \item{Nmarkers}{
Nmarkers number of SNP markers to generate.
}
}
\value{
a matrix of dimensions Nind x Nmarkers.
}
\references{
Wu, R., Ma, C., & Casella, G. (2007). Statistical genetics of quantitative traits: linkage, maps and QTL. Springer Science & Business Media.
}
\author{
Martin Nahuel Garcia <orcid:0000-0001-5760-986X>
}
\examples{
#simulate 100 individuals and 1000 SNPs
set.seed(123)
simGeno(100, 1000)
#[1] "simG was generated"
dim(simG);simG[1:5,1:5]
#[1]  100 1000
#[,1] [,2] [,3] [,4] [,5]
#[1,]    0    1    0    2    2
#[2,]    2    0    2    0    0
#[3,]    1    1    1    2    2
#[4,]    2    2    1    2    1
#[5,]    2    1    1    1    1
}
