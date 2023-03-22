#simG <- simP <- simulatedFS <- simulatedHS <- nsimout <- usimout <- NULL
#' Function to simulate a random SNP matrix, phenotype and QTLs with effects sampled from a Normal distribution
#'
#' @param Nind number of individuals to simulate
#' @param Nmarkers number of SNP markers to generate
#' @param Nqtl number of QTLs to simulate
#' @param Esigma standard deviation of effects with distribution N~(0,Esigma^2)
#' @param Pmean phenotype mean
#' @param Perror standard deviation of error (portion of phenotype not explained by genomic information)
#'
#' @return a list with 4 elements
#' @export
#'
#' @examples
#' set.seed(123)
#' simulN(100, 1000, 50, .9, 12, .5)
#' #[1] "nsimout was generated"
#' str(nsimout)
#' #List of 4
#' #$ geno    : num [1:100, 1:1000] 0 2 1 2 2 0 1 2 1 1 ...
#' #$ pheno   : num [1:100, 1] 25.4 21.6 16 13.8 19.4 ...
#' #$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#' #$ Meffects: num [1:50] 0.2696 -0.1552 1.0192 0.0209 1.2023 ...
simulN <- function(Nind, Nmarkers, Nqtl, Esigma, Pmean, Perror){
  #genotype simulation
  x <- matrix(NA, Nind, Nmarkers)
  for(i in 1:Nmarkers) {
    x[,i] <- round(runif(Nind,-0.5,2.5))
  }
  #SNPs associated (random)
  QTN <- sample(1:Nmarkers, Nqtl, replace=F)
  #SNP effects from uniform distribution
  Neffects <- rnorm(Nqtl,0,Esigma)
  #portion of phenotype not explained by markers
  error <- rnorm(Nind,0,Perror)
  #phenotype simulation
  for(i in 1:Nqtl)
  {Term <- x[,QTN[i]]*Neffects[i]
  assign(paste0("Term", i), Term)
  }
  t2 <- 0
  for(i in 1:Nqtl)
  {t1 <- get(paste0("Term", i))
  t2 <- t1+t2
  }
  y <- Pmean+t2+error
  dim(y) <- c(length(y), 1)
  genphen <- list(geno = x, pheno = y, QTN = QTN, Meffects = Neffects)
  nsimout <- NULL
  nsimout <<- genphen
  return("nsimout was generated")
}

#' Function to simulate half sib progeny genotypes from the genotype of one parent
#'
#' Simulate half sib progeny from one genotyped parent assuming a random genotype for the other parental. We assume that these are diploid organisms.
#'
#' @param x genotype matrix of a set of moms
#' @param Nprogeny number of progeny's genotypes to simulate from each mom
#'
#' @return a matrix of dimensions (nrow(x)*Nprogeny) x ncol(x)
#' @export
#'
#' @examples
#' #' #simulate 100 individuals and 1000 SNPs
#' set.seed(123)
#' simGeno(100, 1000)
#' #[1] "simG was generated"
#' #simulate the genotype of 3 sets 5 HS (one set by mom)
#' simulHS(simG[1:3,],5)
#' #[1] "simulatedHS was generated"
#' dim(simulatedHS)
#' #[1]   15 1000
simulHS <- function(x, Nprogeny){
  #progeny genomic matrix based on mendelian segregation
  newx <- matrix(NA,nrow(x)*Nprogeny,ncol(x))
  g0 <- c(0,1)
  g1 <- c(0,1,1,2)
  g2 <- c(1,2)
  for(i in 1:nrow(x)){
    for(j in 1:Nprogeny){
      for(m in 1:ncol(x)){
        if(x[i,m] == 0){newx[i+nrow(x)*(j-1),m]=sample(g0,1,replace=T)
        } else if(x[i,m] == 1){newx[i+nrow(x)*(j-1),m]=sample(g1,1,replace=T)
        } else if(x[i,m] == 2){newx[i+nrow(x)*(j-1),m]=sample(g2,1,replace=T)
        } else if(is.na(x[i,m])){newx[i+nrow(x)*(j-1),m]=="NA"
        }
      }
    }
  }
  simulatedHS <- NULL
  simulatedHS <<- newx
  return("simulatedHS was generated")
}

#' Function to simulate full sib
#'
#' Simulate full sib progeny genotypes from the genotype of parentals (matrixes with the same dimensions). Pair of parents mating will be in the order of the matrixes. We assume that these are diploid organisms.
#'
#' @param x genotype matrix of a set of moms
#' @param y genotype matrix of a set of dads
#' @param Nprogeny number of progeny's genotypes to generate from each pair of parents
#'
#' @return a matrix of dimensions (nrow(x)*Nprogeny) x ncol(x)
#' @export
#'
#' @examples
#' #simulate 100 individuals and 1000 SNPs
#' set.seed(123)
#' simGeno(100, 1000)
#' #[1] "simG was generated"
#' #simulate the genotype of 5 FS from 3 pairs of parents
#' simulFS(simG[1:3,],simG[4:6,],5)
#' #[1] "simulatedFS was generated"
#'  dim(simulatedFS)
#' #[1]   15 1000
#' # The first 5 individuals are progeny of mom 1 and dad 1, the second 5 individuals are progeny of mom 2 and dad 2, and so on.
simulFS <- function(x, y, Nprogeny){
  #progeny genomic matrix based on mendelian segregation
  newy <- matrix(NA,nrow(x)*Nprogeny,ncol(x))
  g0 <- c(0,1)
  g1 <- c(0,1,1,2)
  g2 <- c(1,2)
  for(i in 1:nrow(x)){
    for(j in 1:Nprogeny){
      for(m in 1:ncol(x)){
        if(x[i,m] == 0 &&  y[i,m] == 0){newy[i+nrow(x)*(j-1),m]=0
        } else if(x[i,m] == 1 &&  y[i,m] == 0){newy[i+nrow(x)*(j-1),m]=sample(g0,1,replace=T)
        } else if(x[i,m] == 2 &&  y[i,m] == 0){newy[i+nrow(x)*(j-1),m]=1
        } else if(x[i,m] == 0 &&  y[i,m] == 1){newy[i+nrow(x)*(j-1),m]=sample(g0,1,replace=T)
        } else if(x[i,m] == 1 &&  y[i,m] == 1){newy[i+nrow(x)*(j-1),m]=sample(g1,1,replace=T)
        } else if(x[i,m] == 2 &&  y[i,m] == 1){newy[i+nrow(x)*(j-1),m]=sample(g2,1,replace=T)
        } else if(x[i,m] == 0 &&  y[i,m] == 2){newy[i+nrow(x)*(j-1),m]=1
        } else if(x[i,m] == 1 &&  y[i,m] == 2){newy[i+nrow(x)*(j-1),m]=sample(g2,1,replace=T)
        } else if(x[i,m] == 2 &&  y[i,m] == 2){newy[i+nrow(x)*(j-1),m]=2
        } else if(is.na(x[i,m]) ||  is.na(y[i,m])){newy[i+nrow(x)*(j-1),m]=="NA"
        }
      }
    }
  }
  simulatedFS <- NULL
  simulatedFS <<- newy
  return("simulatedFS was generated")
}

#' Function to simulate phenotypes
#'
#' Simulate a phenotype from a genotype matrix with QTLs with random effects sampled from a Normal distribution
#'
#' @param x SNP matrix coded like 0 homozygote; 1 heterozygote; 2 homozygote
#' @param Nqtl number of QTLs to simulate
#' @param Esigma standard deviation of effects with distribution N~(0,Esigma^2)
#' @param Pmean phenotype mean
#' @param Perror standard deviation of error (portion of phenotype not explained by genomic information)
#'
#' @return An object of class list containing the trait, the markers associated and their effects.
#' @export
#'
#' @examples
#' set.seed(123)
#' simGeno(100, 1000)
#' #' #[1] "simG was generated"
#' simPheno(simG, 50, .8, 12, .5)
#' #[1] "simP was generated"
#' str(simP)
#' #List of 3
#' #$ pheno   : num [1:100, 1] 24 20.5 15.6 13.6 18.5 ...
#' #$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#' #$ Meffects: num [1:50] 0.2396 -0.138 0.906 0.0186 1.0687 ...
simPheno <- function(x, Nqtl, Esigma, Pmean, Perror){
  #genotype data
  Nind <- nrow(x)
  Nmarkers <- ncol(x)
  #SNPs associated (random)
  QTN <- sample(1:Nmarkers, Nqtl, replace=F)
  #SNP effects from uniform distribution
  Neffects <- rnorm(Nqtl,0,Esigma)
  #portion of phenotype not explained by markers
  error <- rnorm(Nind,0,Perror)
  #phenotype simulation
  for(i in 1:Nqtl)
  {Term <- x[,QTN[i]]*Neffects[i]
  assign(paste0("Term", i), Term)
  }
  t2 <- 0
  for(i in 1:Nqtl)
  {t1 <- get(paste0("Term", i))
  t2 <- t1+t2
  }
  y <- Pmean+t2+error
  dim(y) <- c(length(y), 1)
  pheno <- list(pheno = y, QTN = QTN, Meffects = Neffects)
  simP <- NULL
  simP <<- pheno
  return("simP was generated")
}

#' Function to simulate SNP matrix
#'
#' Simulate SNP matrix coded 0, 1 and 2; with random genotypes
#'
#' @param Nind number of individuals to simulate
#' @param Nmarkers number of SNP markers to generate
#'
#' @return a matrix of dimensions Nind x Nmarkers
#' @export
#'
#' @examples
#' #simulate 100 individuals and 1000 SNPs
#' set.seed(123)
#' simGeno(100, 1000)
#' #[1] "simG was generated"
#' dim(simG);simG[1:5,1:5]
#' #[1]  100 1000
#' #[,1] [,2] [,3] [,4] [,5]
#' #[1,]    0    1    0    2    2
#' #[2,]    2    0    2    0    0
#' #[3,]    1    1    1    2    2
#' #[4,]    2    2    1    2    1
#' #[5,]    2    1    1    1    1
simGeno <- function(Nind, Nmarkers){
  #genotype simulation
  x <- matrix(NA, Nind, Nmarkers)
  for(i in 1:Nmarkers) {
    x[,i] <- round(runif(Nind,-0.5,2.5))
  }
  geno <- x
  simG <- NULL
  simG <<- geno
  return("simG was generated")
}

#' Function to simulate a random SNP matrix, phenotype and QTLs with effects sampled from a Uniform distribution
#'
#' @param Nind number of individuals to simulate
#' @param Nmarkers number of SNP markers to generate
#' @param Nqtl number of QTLs to simulate
#' @param Pmean phenotype mean
#' @param Perror standard deviation of error (portion of phenotype not explained by genomic information)
#'
#' @return a list with 4 elements
#' @export
#'
#' @examples
#' set.seed(123)
#' simulU(100, 1000, 50, 12, .5)
#' #[1] "usimout was generated"
#' str(usimout)
#' #List of 4
#' #$ geno    : num [1:100, 1:1000] 0 2 1 2 2 0 1 2 1 1 ...
#' #$ pheno   : num [1:100, 1] 10.3 14.7 11.8 10.2 13.1 ...
#' #$ QTN     : int [1:50] 568 474 529 349 45 732 416 51 413 514 ...
#' #$ Meffects: num [1:50] 0.2355 0.0158 -0.1369 -0.1246 0.7426 ...
simulU <- function(Nind, Nmarkers, Nqtl, Pmean, Perror){
  #genotype simulation
  x <- matrix(NA, Nind, Nmarkers)
  for(i in 1:Nmarkers) {
    x[,i] <- round(runif(Nind,-0.5,2.5))
  }
  #SNPs associated (random)
  QTN <- sample(1:Nmarkers, Nqtl, replace=F)
  #SNP effects from uniform distribution
  Ueffects <- runif(Nqtl, -1, 1)
  #portion of phenotype not explained by markers
  error <- rnorm(Nind,0,Perror)
  #phenotype simulation
  for(i in 1:Nqtl)
  {Term <- x[,QTN[i]]*Ueffects[i]
  assign(paste0("Term", i), Term)
  }
  t2 <- 0
  for(i in 1:Nqtl)
  {t1 <- get(paste0("Term", i))
  t2 <- t1+t2
  }
  y <- Pmean+t2+error
  dim(y) <- c(length(y), 1)
  genphen <- list(geno = x, pheno = y, QTN = QTN, Meffects = Ueffects)
  usimout <- NULL
  usimout <<- genphen
  return("usimout was generated")
}


