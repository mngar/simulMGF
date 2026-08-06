# Function to simulate phenotypes

Simulate a phenotype from a genotype matrix with QTLs with random
effects sampled from a Normal distribution.

## Usage

``` r
simPheno(x, Nqtl, Esigma, Pmean, Perror)
```

## Arguments

- x:

  SNP matrix coded like 0 homozygous; 1 heterozygous; 2 homozygous

- Nqtl:

  number of QTLs to simulate

- Esigma:

  standard deviation of effects with distribution N~(0,Esigma^2)

- Pmean:

  phenotype mean

- Perror:

  standard deviation of error (portion of phenotype not explained by
  genomic information)

## Value

An object of class list containing the trait, the markers associated and
their effects.

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

## Examples

``` r
set.seed(123)
simGeno(100, 1000)
#> [1] "simG was generated"
#' #[1] "simG was generated"
simPheno(simG, 50, .8, 12, .5)
#> [1] "simP was generated"
#[1] "simP was generated"
str(simP)
#> List of 3
#>  $ pheno   : num [1:100, 1] 24 20.5 15.6 13.6 18.5 ...
#>  $ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#>  $ Meffects: num [1:50] 0.2396 -0.138 0.906 0.0186 1.0687 ...
#List of 3
#$ pheno   : num [1:100, 1] 24 20.5 15.6 13.6 18.5 ...
#$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#$ Meffects: num [1:50] 0.2396 -0.138 0.906 0.0186 1.0687 ...
```
