#' Impute many z-statistics
#'
#' Impute many z-statistics given observed z-statistics and reference panel
#'
#' @param targets variant ID's to impute
#' @param df \code{data.frame} with columns \code{ID}, \code{z}, \code{GWAS_A1}, \code{GWAS_A2}, \code{chrom}, \code{position} \code{A1}, \code{A2}.
#' @param gds \code{GenomicDataStream} of reference panel
#' @param window size of window in bp 
#' @param lambda shrinkage parameter used for LD matrix.  Defaults to \code{lambda = NULL} and is estimated from the data
#' @param method method used to estimate shrinkage parameter lambda.  default is \code{"decorrelate"}
#' @param quiet suppress messages
#'
#' @details Implements method by Pasaniuc, et al. (2014).
#'
#' @references
#' Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., ... & Price, A. L. (2014). Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), 2906-2914.
#'
#' @return \code{data.frame} storing imputed results:
#' \describe{
#'   \item{ID}{variant identifier}
#'   \item{z.stat}{imputed z-statistic}
#'   \item{sigSq}{variance of imputed z-statistic}
#'   \item{r2.pred}{metric of accuracy of the imputed z-statistic based on its variance}
#'   \item{lambda}{shrinkage parameter}
#'   \item{maf}{minor allele frequency in reference panel}
#'   \item{nVariants}{number of variants used in imputation}
#' }
#'
#' @seealso \code{imputez()}
#' @importFrom progress progress_bar
#' @importFrom Rdpack reprompt
#' @importFrom Rfast cora
#' @importFrom dplyr bind_rows as_tibble
#' @importFrom GenomicDataStream setRegion getNextChunk
#' @export
run_analysis = function(targets, df, gds, window = 100000, lambda = NULL, method = c("decorrelate", "Ledoit-Wolf", "OAS", "Touloumis", "Schafer-Strimmer", "Pseudoinverse" ), quiet=FALSE){

	method <- match.arg(method)

	if( ! all(targets %in% df$ID) ){
		stop("All target variants must be in df")
	}

	if (!quiet) {
		pb <- progress_bar$new(
		  format = "  imputing [:bar] :percent eta: :eta",
		  total = length(targets), clear = FALSE, width = 60)
	}

	res = lapply(targets, function(vID){

	  	if (!quiet) pb$tick()

	  	i <- match(vID, df$ID)

	  	# read genotypes in window
		region <- with(df[i,], paste0(chrom, ":", position-window, '-', position+window))
		gds <- setRegion(gds, region)
		dat <- getNextChunk(gds)

		# Subset to matching variants
		df_sub <- df[df$ID %in% dat$info$ID,]

		# keep variants where GWAS and reference have matching alleles
		keep <- with(df_sub, GWAS_A1 == A1 & GWAS_A2 == A2)
		df_sub <- df_sub[keep,]

		# if less than 2 variants in the window
		if( length(df_sub$ID) < 2){
			return(NULL)
		}

		# extract variant data
		X <- dat$X[,df_sub$ID]
		idx <- match(df$ID[i], df_sub$ID)
		z <- df_sub$z
		names(z) <- df_sub$ID

		# impute z-statistic
		if( method == "decorrelate" ){
			res <- imputezDecorr(z, X, idx, lambda = lambda)
		}else{
			lambda <- estimate_lambda(scale(X), method)
			res <- imputez(z, cora(X), idx, lambda = lambda)
		}

		# set MAF
		res$maf <- sum(X[,idx])  / nrow(X)
		res$maf <- pmin(res$maf, 1 - res$maf)
		res$nVariants = ncol(X)

		as_tibble( res )
	})
	res = bind_rows(res)

	if (!quiet) {
		while (!pb$finished) pb$tick()
		pb$terminate()
	}

	res
}
