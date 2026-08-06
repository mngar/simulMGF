# Function to simulate half sib progeny genotypes

Simulate half sib progeny from one genotyped parent assuming a random
genotype for the other parental. We assume that these are diploid
organisms.

## Usage

``` r
simulHS(x, Nprogeny)
```

## Arguments

- x:

  genotype matrix of a set of moms

- Nprogeny:

  number of progeny's genotypes to simulate for each mom

## Details

The function assume: a diploid organism; mendelian segregation of
alleles; and independent segregation.

## Value

a matrix of dimensions (nrow(x)\*Nprogeny) x ncol(x)

## References

Wu, R., Ma, C., & Casella, G. (2007). Statistical genetics of
quantitative traits: linkage, maps and QTL. Springer Science & Business
Media.

## Author

Martin Nahuel Garcia \<orcid:0000-0001-5760-986X\>

## Examples

``` r
#' #simulate 100 individuals and 1000 SNPs
set.seed(123)
simGeno(100, 1000)
#> [1] "simG was generated"
#[1] "simG was generated"
#simulate the genotype of 3 sets 5 HS (one set by mom)
simulHS(simG[1:3,],5)
#> [1] "simulatedHS was generated"
#[1] "simulatedHS was generated"
dim(simulatedHS)
#> [1]   15 1000
#[1]   15 1000
```
