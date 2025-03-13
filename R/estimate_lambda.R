

#' Estimate lambda shrinkage parameter
#' 
#' Estimate lambda shrinkage parameter with specified method
#' 
#' @param X matrix with features as columns
#' @param method method used to estimate lambda shrinkage parameter
#' 
#' @return estimated lambda value
#'
#' @examples
#' n = 100
#' p = 400
#' 
#' X = matrix(rnorm(n*p), n, p)
#' 
#' estimate_lambda(X, "Ledoit-Wolf")
#' 
#' @importFrom CovTools CovEst.2003LW CovEst.2010OAS
#' @importFrom ShrinkCovMat shrinkcovmat.equal
#' @importFrom corpcor estimate.lambda
#' @export
estimate_lambda = function(X, method = c("Ledoit-Wolf", "OAS", "Touloumis", "Schafer-Strimmer" )){

	method <- match.arg(method)

	switch(method, 
		"Ledoit-Wolf" = {
			fit <- CovEst.2003LW( X )
			fit$delta
		},
		"OAS" = {
			fit <- CovEst.2010OAS( X )
			fit$rho
		},
		"Touloumis" = {
			fit <- shrinkcovmat.equal( t(X) )
			fit$lambdahat
		},
		"Schafer-Strimmer" = {
			estimate.lambda(X, verbose=FALSE)
		})
}

