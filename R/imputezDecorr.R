# Gabriel Hoffman
# Oct 4, 2023

#' Impute z-score based correlations
#'
#' Impute z-score from a missing test based on z-scores for other tests, and the genotype matrix
#'
#' @param z vector of observed z-scores
#' @param X matrix genotype data from a reference panel
#' @param i index of z-scores to impute
#' @param k rank of SVD
#' @param lambda (default: NULL) value used to shrink correlation matrix
#'
#' @details Uses implicit covariance and emprical Bayes shrinkage of the eigen-values to accelerate computations. \code{imputez()} is cubic in the number of features, p, but \code{imputezDecorr()}  is the minimum of O(n p^2) and O(n^2p).  For large number of features, this is can be a dramatic speedup.
#' @seealso \code{imputez()}, \code{decorrelate::eclairs()}, \code{decorrelate::decorrelate()}
#'
#' @return \code{data.frame} storing:
#' \describe{
#'   \item{ID}{variant identifier}
#'   \item{z.stat}{imputed z-statistic}
#'   \item{se}{standard error of imputed z-statistic}
#'   \item{r2.pred}{metric of accuracy of the imputed z-statistic based on its variance}
#'   \item{lambda}{shrinkage parameter}
#'   \item{averageCorrSq}{average correlation squared in X}
#' }
#'
#' @examples
#' library(GenomicDataStream)
#' library(mvtnorm)
#' 
#' # VCF file for reference
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#' 
#' # initialize data stream
#' gds <- GenomicDataStream(file, "DS", initialize=TRUE)
#' 
#' # read genotype data from reference
#' dat <- getNextChunk(gds)
#' 
#' # simulate z-statistics with correlation structure
#' # from the LD of the reference panel
#' C <- cor(dat$X)
#' set.seed(1)
#' z <- c(rmvnorm(1, rep(0, 10), C))
#' names(z) <- colnames(dat$X)
#'
#' # Impute z-statistics for variants 2 and 3 
#' # using the other variants and observed z-statistics
#' # from the reference panel
#' # Use dat$X directly instead of creating cor(dat$X)
#' imputezDecorr(z, dat$X, 2:3)
#
#' @references
#' Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., ... & Price, A. L. (2014). Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), 2906-2914.
#' 
#' @importFrom decorrelate eclairs decorrelate averageCorrSq
#' @importFrom Rfast standardise
#' @export
imputezDecorr <- function(z, X, i, k=min(nrow(X), ncol(X)-length(i)), lambda = NULL) {

  stopifnot(length(z) == ncol(X))
  stopifnot(length(i) > 0)

  # standardize cols to have mean zero and var 1
  X <- standardise(X) / sqrt(nrow(X)-1)

  X_exclude <- X[, -i, drop = FALSE]

  ecl <- eclairs(X_exclude, k=k, compute = "correlation", lambda = lambda)

  # if lambda is estiamted and k is full rank
  if( is.null(lambda) & k == min(nrow(X), ncol(X)-length(i)) ){
    # bound lambda 
    r <- min(dim(X))
    ecl$lambda <- max(ecl$lambda, 1/sqrt(r))
  }

  # cross correlation
  # g <- (1 - ecl$lambda) * cor(X_exclude, X[, i,drop=FALSE])
  g <- (1 - ecl$lambda) * crossprod(X_exclude, X[, i,drop=FALSE])

  W <- decorrelate(g, ecl, alpha = -1, transpose = TRUE)

  Sigma_i_t <- 1 - dcrossprod(W, decorrelate(W, ecl, transpose = TRUE, alpha = 1, lambda = 0))
  se <- sqrt(1 - Sigma_i_t)

  # imputed z-statistic 
  z_i <- crossprod(W, z[-i])

  data.frame(
    ID = names(z)[i],
    z.stat = as.numeric(z_i),
    se = as.numeric(se),
    r2.pred = 1 - as.numeric(Sigma_i_t),
    lambda = ecl$lambda,
    averageCorrSq = averageCorrSq(ecl)
  )
}



# diag(crossprod(X,Y))
dcrossprod <- function(X, Y){
  colSums(X * Y)
}




