\name{simPheno}
\alias{simPheno}
\title{
Function to simulate phenotypes
}
\description{
Simulate a phenotype from a genotype matrix with QTLs with random effects sampled from a Normal distribution.
}
\usage{
simPheno(x, Nqtl, Esigma, Pmean, Perror)
}
\arguments{
  \item{x}{
SNP matrix coded like 0 homozygous; 1 heterozygous; 2 homozygous
}
  \item{Nqtl}{
number of QTLs to simulate
}
  \item{Esigma}{
standard deviation of effects with distribution N~(0,Esigma^2)
}
  \item{Pmean}{
phenotype mean
}
  \item{Perror}{
standard deviation of error (portion of phenotype not explained by genomic information)
}
}
\value{
An object of class list containing the trait, the markers associated and their effects.
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
\examples{
set.seed(123)
simGeno(100, 1000)
#' #[1] "simG was generated"
simPheno(simG, 50, .8, 12, .5)
#[1] "simP was generated"
str(simP)
#List of 3
#$ pheno   : num [1:100, 1] 24 20.5 15.6 13.6 18.5 ...
#$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#$ Meffects: num [1:50] 0.2396 -0.138 0.906 0.0186 1.0687 ...
}
