
<br>

## Impute z-statistics for missing variants using observed z-statistics and LD matrix

## Installation
```r
devtools::install_github("GabrielHoffman/imputez")
```
Pasaniuc, et al. (2014).

Uses implicit covariance and emprical Bayes shrinkage of the eigen-values to accelerate computations. \code{imputez()} is cubic in the number of features, p, but \code{imputezDecorr()}  is the minimum of O(n p^2) and O(n^2p).  For large number of features, this is can be a dramatic speedup.

```r
library(imputez)
library(GenomicDataStream)
library(mvtnorm)
library(dplyr)

# VCF file for reference
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# initialize data stream
gds = GenomicDataStream(file, "DS", initialize=TRUE)

# read genotype data from reference
dat = getNextChunk(gds)

# simulate z-statistics with correlation structure
# from the LD of the reference panel
z = c(rmvnorm(1, rep(0, 10), cor(dat$X)))

# Combine z-statistics with variant ID, position, etc
df = dat$info %>%
    mutate(z = z, GWAS_A1 = A1, GWAS_A2 = A2) %>%
    rename(REF_A1 = A1, REF_A2 = A2)

# Given observed z-statistics and 
# GenomicDataStream of reference panel,
# Impute z-statistics from variants missing z-statistics.
# Here drop variant 2, and then impute its z-statistic
res = run_imputez(df[-2,], gds, 10000, 1000)

# Results of imputed z-statistics
res
#> # A tibble: 1 × 9
#>   ID          A1    A2       z.stat        se  r2.pred lambda    maf nVariants
#>   <chr>       <chr> <chr>     <dbl>     <dbl>    <dbl>  <dbl>  <dbl>     <int>
#> 1 1:11000:T:C T     C     0.0000165 0.0000291 8.45e-10   1.00 0.0159         9
```
