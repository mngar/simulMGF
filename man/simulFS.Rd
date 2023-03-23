\name{simulFS}
\alias{simulFS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate full sib progeny genotype
}
\description{
Simulate full sib progeny genotypes from the genotype of the parents (matrixes with the same dimensions). Pair of parents mating will be in the order of the matrixes. We assume that these are diploid organisms.
}
\usage{
simulFS(x, y, Nprogeny)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
genotype matrix of a set of moms
}
  \item{y}{
genotype matrix of a set of dads
}
  \item{Nprogeny}{
Nprogeny number of progeny's genotypes to generate from each pair of parents
}
}
\value{
a matrix of dimensions (nrow(x)*Nprogeny) x ncol(x)
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
#simulate the genotype of 5 FS from 3 pairs of parents
simulFS(simG[1:3,],simG[4:6,],5)
#[1] "simulatedFS was generated"
 dim(simulatedFS)
#[1]   15 1000
# The first 5 individuals are progeny of mom 1 and dad 1, the second 5 individuals
# are progeny of mom 2 and dad 2, and so on.
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
