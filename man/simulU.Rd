\name{simulU}
\alias{simulU}
\title{
Function to simulate a random SNP matrix, phenotype and QTLs with their effects
}
\description{
This function simulate a SNP matrix (coded as 0, 1, 2) and traits with a selected number of QTLs and their effects that will be sampled from a Uniform distribution.
}
\usage{
simulU(Nind, Nmarkers, Nqtl, Pmean, Perror)
}
\arguments{
  \item{Nind}{
number of individuals to simulate.
}
  \item{Nmarkers}{
number of SNP markers to generate.
}
  \item{Nqtl}{
number of QTLs controlling the trait.
}
  \item{Pmean}{
phenotype mean.
}
  \item{Perror}{
standard deviation of error (portion of phenotype not explained by genomic information).
}
}
\value{
An object of class list containing the SNP matrix, the trait, the markers associated and their effects.
  \item{geno }{SNP matrix generated.}
  \item{pheno }{vector with the trait values simulated.}
  \item{QTN }{column in the SNP matrix with the SNP associated.}
  \item{Meffects }{effects of the associated SNPs.}
}
\references{
Wu, R., Ma, C., & Casella, G. (2007). Statistical genetics of quantitative traits: linkage, maps and QTL. Springer Science & Business Media.
}
\author{
Martin Nahuel Garcia <orcid:0000-0001-5760-986X>
}
\seealso{
simGeno, simulN
}
\examples{
set.seed(123)
simulU(100, 1000, 50, 12, .5)
#[1] "usimout was generated"
str(usimout)
#List of 4
#$ geno    : num [1:100, 1:1000] 0 2 1 2 2 0 1 2 1 1 ...
#$ pheno   : num [1:100, 1] 10.3 14.7 11.8 10.2 13.1 ...
#$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#$ Meffects: num [1:50] 0.2355 0.0158 -0.1369 -0.1246 0.7426 ...
}
