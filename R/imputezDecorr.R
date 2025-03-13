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
#' }
#'
#' @references
#' Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., ... & Price, A. L. (2014). Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), 2906-2914.
#' 
#' @importFrom decorrelate eclairs decorrelate
#' @importFrom stats cor
#' @export
imputezDecorr <- function(z, X, i, k=ncol(X), lambda = NULL) {

  stopifnot(length(z) == ncol(X))
  stopifnot(length(i) > 0)

  X_exclude = X[, -i, drop = FALSE]

  ecl <- eclairs(X_exclude, k=k, compute = "correlation", lambda = lambda)

  g <- (1 - ecl$lambda) * cor(X_exclude, X[, i,drop=FALSE])

  W <- decorrelate(g, ecl, alpha = -1, transpose = TRUE)

  z_i <- crossprod(W, z[-i])

  Sigma_i_t <- 1 - dcrossprod(W, decorrelate(W, ecl, transpose = TRUE, alpha = 1, lambda = 0))

  se = sqrt(1 - Sigma_i_t)

  data.frame(
    ID = names(z)[i],
    z.stat = as.numeric(z_i),
    se = as.numeric(se),
    r2.pred = 1 - as.numeric(Sigma_i_t),
    lambda = ecl$lambda
  )
}

# diag(crossprod(X,Y))
dcrossprod = function(X, Y){
  colSums(X * Y)
}

