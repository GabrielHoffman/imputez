
#' Impute z-score based correlations
#' 
#' Impute z-score from a missing test based on z-scores for other tests, and the correlation matrix between z-scores
#' 
#' @param z vector of observed z-scores
#' @param Sigma matrix of correlation between z-scores
#' @param i index of z-scores to impute
#' @param lambda value used to shrink correlation matrix
#' 
#' @details Implements method by Pasaniuc, et al. (2014).
#'
#' @return \code{data.frame} with variant \code{ID}, imputed \code{z.stat} for test \code{i}, variance \code{sigSq} of the impute z-score, and metric of accuracy \code{r2.pred}.
#' 
#' @references{
#'   \insertRef{pasaniuc2014fast}{imputez}
#' }
#' @importFrom Matrix Diagonal
#' @importFrom methods is
#' @export
imputez = function(z, Sigma, i, lambda = 0.1){
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
	sigSq = Sigma.shrink[i,i] - (crossprod(W,Sigma.shrink[-i, -i]) %*% W)
	
	data.frame(	ID = names(z)[i], 
				z.stat = as.numeric(z_i), 
				sigSq = as.numeric(sigSq),
				r2.pred = 1 - as.numeric(sigSq))
}


# dcmp = eigen(Sigma, symmetric=TRUE)
# Sigma.shrink = with(dcmp, tcrossprod(vectors, vectors %*% diag(values)))
# Sigma.shrink[1:3, 1:3]

# Sigma.shrink = with(dcmp, tcrossprod(vectors, vectors %*% diag(pmax(0, values))))
# Sigma.shrink[1:3, 1:3]

#' Get complete subset of correlation matrix
#' 
#' Get complete subset of correlation matrix so that no entries are zero
#' 
#' @param S correlation matrix
#' @param mid index of target variant
#' 
#' @return begin and end indeces of completely filled matrixq
#' 
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
constructLD = function(dfld, ids){

	idx_A = idx_B = SNP_A = NULL

	inclgd = expand.grid(ids, ids)

	dfldsub = dfld[.(inclgd$Var1, inclgd$Var2), mult = "first", nomatch = 0]

	# get ordering of LD matrix subset
	IDs = dfldsub[,unique(SNP_A)]
	dfldsub$idx_A = dfldsub[,match(SNP_A, IDs)]
	dfldsub$idx_B = dfldsub[,match(SNP_B, IDs)]

	N = length(IDs)

	C = sparseMatrix( i = dfldsub$idx_A - min(dfldsub$idx_A) + 1,
							j = dfldsub$idx_B - min(dfldsub$idx_B) + 1,
							x = dfldsub$R,
							symmetric = TRUE, 
							dims = c(N,N), 
							dimnames=list(IDs, IDs) )
	C
}

#' Impute many z-statistics 
#' 
#' Impute many z-statistics given observed z-statistics and \code{data.table} of LD.
#' 
#' @param z vector of observed z-statistics with unobserved storing value \code{NA}
#' @param dfld \code{data.table} storing LD information
#' @param IDs variant names in \code{z} to impute
#' @param maxWindowSize max window size around target variant to keep
#' @param quiet default FALSE.  If TRUE, suppress progress bar
#' 
#' @details Implements method by Pasaniuc, et al. (2014).
#' 
#' @references{
#'   \insertRef{pasaniuc2014fast}{imputez}
#' }
#' 
#' @return \code{data.frame} storing imputed results
#' 
#' @importFrom progress progress_bar
#' @importFrom Rdpack reprompt
#' @export
run_imputez = function( z, dfld, IDs, maxWindowSize = 200, quiet=FALSE){

	if( ! quiet ){
		pb <- progress_bar$new(
			format = "  imputing [:bar] :percent eta: :eta",
			total = length(IDs), clear = FALSE, width= 60)
	}

	# precompute
	b = cumsum(!is.na(z))
	df_z = lapply(IDs, function(id){
		i = match(id, names(z))
		
		# if variant is not in LD reference
		if( dfld[.(id, id),is.na(R)]) return(NULL)

		# get window including maxWindowSize observed z-scores
		incl = get_window(z, i, id, maxWindowSize,b)
		z_local = z[incl]
		Sigma_local = as.matrix( constructLD(dfld, names(z)[incl]))

		# get only shared variants
		keep = intersect( names(z_local), rownames(Sigma_local))
		z_local = z_local[keep]
		Sigma_local = Sigma_local[keep,keep]

		# Reduce window until most distance SNPs have LD computed
		# This ensures that Sigma is positive definite
		mid = which(names(z_local) == id)
		window = pairwiseCompleteWindow( Sigma_local, mid)

		incl2 = seq(window[1], window[2])
		Sigma_local = Sigma_local[incl2, incl2]
		z_local = z_local[incl2]

		# keep only variants with observed z-score 
		# or the target variant to be imputed
		keep = unique(c(id, names(z_local)[!is.na(z_local)]))
		Sigma_local = Sigma_local[keep,keep]
		z_local = z_local[keep]

		k = which(names(z_local) == id)
		df = imputez(z_local, Sigma_local, k, lambda=0.1)

		if( ! quiet ) pb$tick()

		data.frame(df, width=window[2] - window[1])
	})

	if( ! quiet ) pb$terminate()

	do.call(rbind, df_z)
}

# replace 

# incl = seq(pmax(1, i-maxWindowSize), 
# 			pmin(length(z), i + maxWindowSize))
# get_window = function(z, i, id, maxWindowSize){
# 	a = cumsum(!is.na(z[seq(1, i)]))
# 	b = cumsum(!is.na(z[seq(i,length(z))]))

# 	id1 = which.min(abs(a -(a[id] - maxWindowSize)))
# 	id2 = which.min(abs(b - maxWindowSize))

# 	idx = which(names(z) %in% c(names(id1), names(id2)))
# 	seq(idx[1], idx[2])
# }

# b = cumsum(!is.na(z[seq(i,length(z))]))
# which.min(abs(b - maxWindowSize))

get_window = function(z, i, id, maxWindowSize, 
	b = cumsum(!is.na(z))){

	id1 = which.min(abs(b[seq(1,i)] -(b[id] - maxWindowSize)))
	j = which.min(abs(b - maxWindowSize)) + i - 1
	id2 = names(b)[j]

	idx = which(names(z) %in% c(names(id1), id2))
	seq(idx[1], idx[2])
}

#  range(get_window(z, i, id, maxWindowSize=50))
#  range(get_window2(z, i, id, maxWindowSize=50))



# system.time(replicate(100, get_window(z, i, id, maxWindowSize)))


# b = cumsum(!is.na(z))
# system.time(replicate(100, get_window2(z, i, id, maxWindowSize)))

# system.time(replicate(100, get_window2(z, i, id, maxWindowSize, b)))




