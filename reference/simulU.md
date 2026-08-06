# Function to simulate a random SNP matrix, phenotype and QTLs with their effects

This function simulate a SNP matrix (coded as 0, 1, 2) and traits with a
selected number of QTLs and their effects that will be sampled from a
Uniform distribution.

## Usage

``` r
simulU(Nind, Nmarkers, Nqtl, Pmean, Perror)
```

## Arguments

- Nind:

  number of individuals to simulate.

- Nmarkers:

  number of SNP markers to generate.

- Nqtl:

  number of QTLs controlling the trait.

- Pmean:

  phenotype mean.

- Perror:

  standard deviation of error (portion of phenotype not explained by
  genomic information).

## Value

An object of class list containing the SNP matrix, the trait, the
markers associated and their effects.

- geno :

  SNP matrix generated.

- pheno :

  vector with the trait values simulated.

- QTN :

  column in the SNP matrix with the SNP associated.

- Meffects :

  effects of the associated SNPs.

## References

Wu, R., Ma, C., & Casella, G. (2007). Statistical genetics of
quantitative traits: linkage, maps and QTL. Springer Science & Business
Media.

## Author

Martin Nahuel Garcia \<orcid:0000-0001-5760-986X\>

## See also

simGeno, simulN

## Examples

``` r
set.seed(123)
simulU(100, 1000, 50, 12, .5)
#> [1] "usimout was generated"
#[1] "usimout was generated"
str(usimout)
#> List of 4
#>  $ geno    : num [1:100, 1:1000] 0 2 1 2 2 0 1 2 1 1 ...
#>  $ pheno   : num [1:100, 1] 10.3 14.7 11.8 10.2 13.1 ...
#>  $ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#>  $ Meffects: num [1:50] 0.2355 0.0158 -0.1369 -0.1246 0.7426 ...
#List of 4
#$ geno    : num [1:100, 1:1000] 0 2 1 2 2 0 1 2 1 1 ...
#$ pheno   : num [1:100, 1] 10.3 14.7 11.8 10.2 13.1 ...
#$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#$ Meffects: num [1:50] 0.2355 0.0158 -0.1369 -0.1246 0.7426 ...
```
