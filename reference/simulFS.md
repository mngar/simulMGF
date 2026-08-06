# Function to simulate full sib progeny genotype

Simulate full sib progeny genotypes from the genotype of the parents
(matrixes with the same dimensions). Pair of parents mating will be in
the order of the matrixes. We assume that these are diploid organisms.

## Usage

``` r
simulFS(x, y, Nprogeny)
```

## Arguments

- x:

  genotype matrix of a set of moms

- y:

  genotype matrix of a set of dads

- Nprogeny:

  Nprogeny number of progeny's genotypes to generate from each pair of
  parents

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
#simulate 100 individuals and 1000 SNPs
set.seed(123)
simGeno(100, 1000)
#> [1] "simG was generated"
#[1] "simG was generated"
#simulate the genotype of 5 FS from 3 pairs of parents
simulFS(simG[1:3,],simG[4:6,],5)
#> [1] "simulatedFS was generated"
#[1] "simulatedFS was generated"
 dim(simulatedFS)
#> [1]   15 1000
#[1]   15 1000
# The first 5 individuals are progeny of mom 1 and dad 1, the second 5 individuals
# are progeny of mom 2 and dad 2, and so on.
```
