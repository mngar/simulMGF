\name{simulN}
\alias{simulN}
\title{
Function to simulate a SNP matrix, a phenotypic trait and the effects of associated SNPs.
}
\description{
This function simulate a SNP matrix (coded as 0, 1, 2) and traits with a selected number of QTLs and their effects that will be sampled from a Normal distribution.
}
\usage{
simulN(Nind, Nmarkers, Nqtl, Esigma, Pmean, Perror)
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
  \item{Esigma}{
standard deviation of effects with distribution N~(0,Esigma^2).
}
  \item{Pmean}{
phenotype mean.
}
  \item{Perror}{
standard deviation of error (portion of phenotype not explained by genomic information).
}
}
\details{
Genotypic data is simulated as the round value sampled from an uniform distribution with interval (-.5,2.5). Phenotypic data are obtained as a linear function defined by:
\deqn{y = Pmean + \sum QTN*Meffects + Perror }

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
\note{
  The genotype is simulated in the same way of simGeno function.
  The trait, QTLs and their effects are simulated in the same way of simPheno function.
}
\seealso{
simGeno, simPheno, simulU
}
\examples{
set.seed(123)
simulN(100, 1000, 50, .9, 12, .5)
#[1] "nsimout was generated"
str(nsimout)
#List of 4
#$ geno    : num [1:100, 1:1000] 0 2 1 2 2 0 1 2 1 1 ...
#$ pheno   : num [1:100, 1] 25.4 21.6 16 13.8 19.4 ...
#$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#$ Meffects: num [1:50] 0.2696 -0.1552 1.0192 0.0209 1.2023 ...
}
