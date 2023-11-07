# Gabriel Hoffman
# Oct 4, 2023

# #' Impute z-score based correlations
# #' 
# #' Impute z-score from a missing test based on z-scores for other tests, and the genotype matrix 
# #' 
# #' @param z vector of observed z-scores
# #' @param X matrix genotype data from a reference panel
# #' @param i index of z-scores to impute
# #' @param lambda (default: NULL) value used to shrink correlation matrix
# #' 
# #' @details Uses implicit covariance and emprical Bayes shrinkage of the eigen-values to accelerate computations. \code{imputez()} is cubic in the number of features, p, but \code{imputez.eclairs()}  is the minimum of O(n p^2) and O(n^2p).  For large number of features, this is can be a dramatic speedup.
# #' @seealso \code{imputez}
# #'
# #' @return \code{data.frame} storing:
# #' \describe{
# #'   \item{ID}{variant identifier}
# #'   \item{z.stat}{imputed z-statistic}
# #'   \item{sigSq}{variance of imputed z-statistic}
# #'   \item{r2.pred}{metric of accuracy of the imputed z-statistic based on its variance}
# #'   \item{lambda}{correlation shrinkage parameter)
# #' }
# #' 
# #' @references{
# #'   \insertRef{pasaniuc2014fast}{imputez}
# #' }
# #' @importFrom decorrelate eclairs decorrelate
# #' @export
# imputez.eclairs = function(z, X, i, lambda = NULL){

# 	ecl = eclairs(X[,-i], compute="correlation", lambda = lambda)

# 	g = (1-ecl$lambda) * cor(X[,-i,drop=FALSE], X[,i])

# 	W = decorrelate(g, ecl, alpha = -1, transpose = TRUE)

# 	z_i = crossprod(W, z[-i])

# 	sigSq = 1 - crossprod(W, decorrelate(W, ecl, transpose=TRUE, alpha=1))

# 	data.frame(	ID = names(z)[i], 
# 				z.stat = as.numeric(z_i), 
# 				sigSq = as.numeric(sigSq),
# 				r2.pred = 1 - as.numeric(sigSq),
# 				lambda = ecl$lambda)
# }