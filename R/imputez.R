#' Impute z-score based correlations
#'
#' Impute z-score from a missing test based on z-scores for other tests, and the correlation matrix between z-scores
#'
#' @param z vector of observed z-scores
#' @param Sigma matrix of correlation between z-scores
#' @param i index of z-scores to impute
#' @param lambda value used to shrink correlation matrix
#' @param useginv if \code{TRUE} use pseudoinverse
#'
#' @details Implements method by Pasaniuc, et al. (2014).
#'
#' @return \code{data.frame} storing:
#' \describe{
#'   \item{ID}{variant identifier}
#'   \item{z.stat}{imputed z-statistic}
#'   \item{se}{standard error of imputed z-statistic}
#'   \item{r2.pred}{metric of accuracy of the imputed z-statistic based on its variance}
#'   \item{lambda}{shrinkage parameter}
#' }
#'
#' @examples
#' library(GenomicDataStream)
#' library(mvtnorm)
#' library(dplyr)
#' 
#' # VCF file for reference
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#' 
#' # initialize data stream
#' gds = GenomicDataStream(file, "DS", initialize=TRUE)
#' 
#' # read genotype data from reference
#' dat = getNextChunk(gds)
#' 
#' # simulate z-statistics with correlation structure
#' # from the LD of the reference panel
#' C = cor(dat$X)
#' z = c(rmvnorm(1, rep(0, 10), C))
#' names(z) = colnames(dat$X)
#'
#' # Impute z-statistics for variants 2 and 3 
#' # using the other variants and observed z-statistics
#' # from the reference panel
#' imputez(z, C, 2:3)
#
#' @references
#' Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., ... & Price, A. L. (2014). Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), 2906-2914.
#'
#' @importFrom Matrix Diagonal
#' @importFrom MASS ginv
#' @importFrom methods is
#' @export
imputez <- function(z, Sigma, i, lambda = 0.1, useginv=FALSE) {

  stopifnot(nrow(Sigma) == ncol(Sigma))
  stopifnot(length(z) == ncol(Sigma))
  stopifnot(length(i) > 0)
  stopifnot(lambda >=0 & lambda <= 1)

  if (is(Sigma, "sparseMatrix")) {
    Sigma.shrink <- (1 - lambda) * Sigma + Diagonal(nrow(Sigma), lambda)
  } else if (is(Sigma, "matrix")) {
    Sigma.shrink <- (1 - lambda) * Sigma + diag(lambda, nrow(Sigma))
  } else {
    txt <- paste("Sigma must be a matrix or sparseMatrix. Class not supported:", class(Sigma)[1])
    stop(txt)
  }

  # impute the ith z-score using
  # impute z-score using Gaussian conditional distribution
  # z_i = Sigma.shrink[i,-i] %*% solve(Sigma.shrink[-i,-i], z[-i])
  # W = Sigma.shrink[i,-i] %*% solve(Sigma.shrink[-i,-i])
  # weights
  A <- Sigma.shrink[-i, -i]
  b <- Sigma.shrink[-i, i, drop = FALSE]
  if( useginv ){
    W <- ginv(A) %*% b
  }else{
    W <- solve(A, b)
  }

  # imputed z-scores
  z_i <- crossprod(W, z[-i])

  # compute standard error for each imputed z-score
  # Sigma_i_t <- Sigma.shrink[i, i] - Sigma.shrink[i,-i] %*% solve(Sigma.shrink[-i, -i], Sigma.shrink[i,-i])
  # Sigma_i_t <- Sigma[i, i] - (dcrossprod(W, Sigma[-i, -i]) %*% W)
  Sigma_i_t <- 1 - dcrossprod(W, crossprod(Sigma[-i, -i], W))

  # standard error of z-statistic
  se = sqrt(1 - Sigma_i_t)

  data.frame(
    ID = names(z)[i],
    z.stat = as.numeric(z_i),
    se = as.numeric(se),
    r2.pred = 1 - as.numeric(Sigma_i_t),
    lambda = lambda
  )
}



