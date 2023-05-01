
#' Impute z-score based correlations
#' 
#' Impute z-score from a missing test based on z-scores for other tests, and the correlation matrix between z-scores
#' 
#' @param z vector of observed z-scores
#' @param Sigma matrix of correlation between z-scores
#' @param i index of z-scores to impute
#' @param lambda value used to shrink correlation matrix
#' 
#' @return \code{data.frame} with variant \code{ID}, imputed \code{z.stat} for test \code{i}, standard error \code{se} of the impute z-score, and metric of accuracy \code{r2.pred}.
#' 
#' @importFrom Matrix Diagonal
#' @export
impute_z = function(z, Sigma, i, lambda = 0.1){
	if( is(Sigma, "sparseMatrix") ){
		Sigma.shrink = (1-lambda) * Sigma + Diagonal(nrow(Sigma), lambda)
	}else{
		Sigma.shrink = (1-lambda) * Sigma + diag(lambda, nrow(Sigma))
	}
	# impute the ith z-score using 
	# impute z-score using Gaussian conditional distribution
	# z_i = Sigma.shrink[i,-i] %*% solve(Sigma.shrink[-i,-i], z[-i])
	# W = Sigma.shrink[i,-i] %*% solve(Sigma.shrink[-i,-i])
	# weights
	W = solve(Sigma.shrink[-i,-i], Sigma.shrink[-i,i,drop=FALSE])

	# imputed z-scores
	z_i = crossprod(W, z[-i])

	# compute standard error for each imputed z-score
	sig = Sigma.shrink[i,i] - (crossprod(W,Sigma.shrink[-i, -i]) %*% W)
	
	data.frame(ID = names(z)[i], 
						z.stat = as.numeric(z_i), 
						se = as.numeric(sig),
						r2.pred = 1 - as.numeric(sig))
}



# Instead do reverse diagonal
#' @export
pairwiseCompleteWindow = function(S, mid){
	# Reduce window until most distance SNPs have LD computed
	a = 1
	b = ncol(S)
	alterate = TRUE

	while(S[a,b] == 0){
		if( alterate )
			a <- pmin(a + 1, mid)
		else
			b <- pmax(b - 1, mid)
		alterate = ! alterate
	}

	c(a,b)
}

#' Construct LD matrix from data.table
#' 
#' Construct LD matrix from data.table
#' 
#' @param dfld \code{data.table} storing LD information
#' @param incl indeces to include 
#' 
#' @return \code{sparseMatrix} storing LD between variants
#' @importFrom Matrix sparseMatrix
#' @export
constructLD = function(dfld, incl){

	inclgd = expand.grid(incl, incl)

	dfldsub = dfld[.(inclgd$Var1, inclgd$Var2), mult = "first", nomatch = NULL]

	rng = dfldsub[,range(idx_A, idx_B)]
	N = rng[2] - rng[1] + 1
	IDs = dfldsub[,unique(SNP_A)]

	C = sparseMatrix( i = dfldsub$idx_A - min(dfldsub$idx_A) + 1,
							j = dfldsub$idx_B - min(dfldsub$idx_B) + 1,
							x = dfldsub$R,
							symmetric = TRUE, 
							dims = c(N,N), 
							dimnames=list(IDs, IDs) )
	C
}









