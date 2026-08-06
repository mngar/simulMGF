# Function to simulate SNP matrix

Simulate SNP matrix coded 0, 1 and 2; with random genotypes.

## Usage

``` r
simGeno(Nind, Nmarkers)
```

## Arguments

- Nind:

  number of individuals to simulate.

- Nmarkers:

  Nmarkers number of SNP markers to generate.

## Value

a matrix of dimensions Nind x Nmarkers.

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
dim(simG);simG[1:5,1:5]
#> [1]  100 1000
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    2    2
#> [2,]    2    0    2    0    0
#> [3,]    1    1    1    2    2
#> [4,]    2    2    1    2    1
#> [5,]    2    1    1    1    1
#[1]  100 1000
#[,1] [,2] [,3] [,4] [,5]
#[1,]    0    1    0    2    2
#[2,]    2    0    2    0    0
#[3,]    1    1    1    2    2
#[4,]    2    2    1    2    1
#[5,]    2    1    1    1    1
```
