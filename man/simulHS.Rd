\name{simulHS}
\alias{simulHS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate half sib progeny genotypes
}
\description{
Simulate half sib progeny from one genotyped parent assuming a random genotype for the other parental. We assume that these are diploid organisms.
}
\usage{
simulHS(x, Nprogeny)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
genotype matrix of a set of moms
}
  \item{Nprogeny}{
number of progeny's genotypes to simulate for each mom
}
}
\details{
The function assume: a diploid organism; mendelian segregation of alleles; and independent segregation.
}
\value{
a matrix of dimensions (nrow(x)*Nprogeny) x ncol(x)
}
\references{
Wu, R., Casella, G., & Ma, C. X. (2007). Linkage Analysis for Distorted and Misclassified Markers. Statistical Genetics of Quantitative Traits: Linkage, Maps, and QTL, 123-143.
}
\author{
Martin Nahuel Garcia <orcid:0000-0001-5760-986X>
}
\examples{
#' #simulate 100 individuals and 1000 SNPs
set.seed(123)
simGeno(100, 1000)
#[1] "simG was generated"
#simulate the genotype of 3 sets 5 HS (one set by mom)
simulHS(simG[1:3,],5)
#[1] "simulatedHS was generated"
dim(simulatedHS)
#[1]   15 1000
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
