
#' Impute many z-statistics
#'
#' Impute many z-statistics given observed z-statistics and reference panel
#'
#' @param df \code{data.frame} with columns \code{ID}, \code{z}, \code{GWAS_A1}, \code{GWAS_A2}, \code{CHROM}, \code{POS} \code{REF_A1}, \code{REF_A2}.
#' @param gds \code{GenomicDataStream} of reference panel
#' @param region genomic region to impute
#' @param flankWidth additional window added to \code{region} 
#' @param method method used to estimate shrinkage parameter lambda.  default is \code{"decorrelate"}
#' @param lambda (default: NULL) value used to shrink correlation matrix. Only used if method is \code{"decorrelate"}
#' 
#' @return \code{tibble} storing imputed results:
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
#' z = c(rmvnorm(1, rep(0, 10), cor(dat$X)))
#' 
#' # Combine z-statistics with variant ID, position, etc
#' df = dat$info %>%
#' 		mutate(z = z, GWAS_A1 = A1, GWAS_A2 = A2) %>%
#' 		rename(REF_A1 = A1, REF_A2 = A2)
#' 
#' # Given observed z-statistics and 
#' # GenomicDataStream of reference panel,
#' # Impute z-statistics from variants missing z-statistics.
#' # Here drop variant 2, and then impute its z-statistic
#' # Impute variants in the given region
#' region = "1:1000-100000"
#' impute_region(df[-2,], gds, region, 1000)
#
#' @importFrom GenomicDataStream setRegion getNextChunk
#' @importFrom Rfast standardise colsums
#' @importFrom dplyr as_tibble
#' @export
impute_region = function(df, gds, region, flankWidth, method = c("decorrelate", "Ledoit-Wolf", "OAS", "Touloumis", "Schafer-Strimmer" ), lambda = NULL){

	method <- match.arg(method)

	# expand region to include flanking width
	loc <- strsplit(region, ":|-")[[1]]
	start <- as.numeric(loc[2])
	end <- as.numeric(loc[3])
	region_expanded <- paste0(loc[1], ":", 
		max(0, start - flankWidth), 
		"-",
		end + flankWidth)

	# read genotype data from this region
	gds2 <- setRegion(gds, region_expanded)
	dat <- getNextChunk( gds2 )

	# if no variants found, return null
	if( length(dat) == 0){
		return(NULL)
	}

	# Modify names to include alleles
	df$ID <- with(df, paste(ID, GWAS_A1, GWAS_A2, sep="__"))
	dat$info$ID <- with(dat$info, paste(ID, A1, A2, sep="__"))
	colnames(dat$X) <- dat$info$ID

	# Subset to matching variants
	df_sub <- df[df$ID %in% dat$info$ID,]

	# keep variants where GWAS and reference have matching alleles
	keep <- with(df_sub, GWAS_A1 == REF_A1 & GWAS_A2 == REF_A2)
	df_sub <- df_sub[keep,]

	# mark variants within original region
	dat$info$inRegion <- (dat$info$POS >= start) & (dat$info$POS <= end)

	# create array of z-statisics with unobserved values as NA
	zstats <- rep(NA, ncol(dat$X))
	names(zstats) <- colnames(dat$X)
	zstats[df_sub$ID] <- df_sub$z

	# keep variants that are either 
	# 1) in the target region
	# 2) outside the target and observed
	keep <- dat$info$inRegion | (!dat$info$inRegion & !is.na(zstats))

	z <- zstats[keep]
	X <- dat$X[,keep,drop=FALSE]

	# Remove entries with duplicated names
	d <- duplicated(colnames(X))
	keep <- ! colnames(X) %in% unique(colnames(X)[d])
	X <- X[,keep,drop=FALSE]
	z <- z[keep]

	# get indeces of unobserved z-statistics
	idx <- which(is.na(z))

	if( (length(idx) == 0) | (ncol(X) - length(idx) < 3) ){
		return(NULL)
	}

	# impute z-statistic
	if( method == "decorrelate" ){
		res <- imputezDecorr(z, X, idx, lambda = lambda)
	}else{
		X_scaled = standardise(X)
		lambda <- estimate_lambda(X_scaled, method)
		C <- crossprod(X_scaled) / c(nrow(X)-1)
		res <- imputez(z, C, idx, lambda = lambda)
	}

	# Set alleles
	df_alleles <- do.call(rbind, strsplit(res$ID, "__"))
	res$ID <- df_alleles[,1]
	res$A1 <- df_alleles[,2]
	res$A2 <- df_alleles[,3]
	cols <- c("ID", "A1", "A2", "z.stat", "se", "r2.pred", "lambda")
	res <- res[,cols]

	# set MAF
	res$maf <- colsums(X[,idx,drop=FALSE]) / (2*nrow(X))
	res$maf <- pmin(res$maf, 1 - res$maf)
	res$nVariants <- ncol(X) - length(idx)

	as_tibble( res )
}

#' @importFrom dplyr group_by summarize `%>%`
get_analysis_windows = function(df, window){

	CHROM <- POS <- NULL 

	# get range of positions on each chromosome
	df_chrom <- df %>%
		group_by(CHROM) %>%
		summarize(start = min(POS), end = max(POS))

	# for each chrom
	regions <- vapply(seq(nrow(df_chrom)), function(i){

		# get width in bp and number of regions
		# to divide the chrom into to get close 
		# to the target window
		width <- df_chrom$end[i] - df_chrom$start[i]
		nregions <- max(2, floor(width / window))

		# get positions
		start <- seq(df_chrom$start[i], df_chrom$end[i], length.out=nregions)
		start <- floor(start)

		# create coordinate windows
		paste0(df_chrom$CHROM[i], ":", start[-length(start)], "-", start[-1])
	}, FUN.VALUE="character"))	

	c(regions)
}

#' Impute many z-statistics
#'
#' Impute many z-statistics given observed z-statistics and reference panel
#'
#' @param df \code{data.frame} with columns \code{ID}, \code{z}, \code{GWAS_A1}, \code{GWAS_A2}, \code{CHROM}, \code{POS} \code{A1}, \code{A2}.
#' @param gds \code{GenomicDataStream} of reference panel
#' @param window size of window in bp 
#' @param flankWidth additional window added to \code{region} 
#' @param method method used to estimate shrinkage parameter lambda.  default is \code{"decorrelate"}
#' @param lambda (default: NULL) value used to shrink correlation matrix. Only used if method is \code{"decorrelate"}
#' @param quiet suppress messages
#'
#' @details Implements method by Pasaniuc, et al. (2014).
#'
#' @references
#' Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., ... & Price, A. L. (2014). Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), 2906-2914.
#'
#' @return \code{tibble} storing imputed results:
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
#' z = c(rmvnorm(1, rep(0, 10), cor(dat$X)))
#' 
#' # Combine z-statistics with variant ID, position, etc
#' df = dat$info %>%
#' 		mutate(z = z, GWAS_A1 = A1, GWAS_A2 = A2) %>%
#' 		rename(REF_A1 = A1, REF_A2 = A2)
#' 
#' # Given observed z-statistics and 
#' # GenomicDataStream of reference panel,
#' # Impute z-statistics from variants missing z-statistics.
#' # Here drop variant 2, and then impute its z-statistic
#' res = run_imputez(df[-2,], gds, 10000, 1000)
#' 
#' # Results of imputed z-statistics
#' res
#
#' @seealso \code{imputez()}
#' @importFrom progress progress_bar
#' @importFrom GenomicDataStream setChunkSize
#' @importFrom dplyr bind_rows
#' @export
run_imputez = function(df, gds, window, flankWidth, method = c("decorrelate", "Ledoit-Wolf", "OAS", "Touloumis", "Schafer-Strimmer"), lambda = NULL, quiet=FALSE){

	method <- match.arg(method)
	gds <- setChunkSize(gds, 1e9)

	cols = c("ID", "z", "GWAS_A1", "GWAS_A2", "CHROM", "POS", "REF_A1", "REF_A2")
	if( ! any(cols %in% colnames(df)) ){
		stop("df must have colnames: ID, z, GWAS_A1, GWAS_A1, CHROM, POS, REF_A1, REF_A2")
	}

	regions <- get_analysis_windows( df, window)

	if (!quiet) {
		pb <- progress_bar$new(
		  format = "  imputing [:bar] :percent eta: :eta",
		  total = length(regions), clear = FALSE, width = 60)
	}

	# impute in each region
	res <- lapply(regions, function(region){
	  	if (!quiet) pb$tick()
		impute_region(df, gds, region, flankWidth, method, lambda)
	})
	res <- bind_rows(res)

	if (!quiet) {
		while (!pb$finished) pb$tick()
		pb$terminate()
	}

	res
}




