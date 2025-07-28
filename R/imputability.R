

#' @importFrom GenomicDataStream setRegion getNextChunk
#' @importFrom stats cor
imputability_region <- function(gds, region, flankWidth, ids, method = c("decorrelate", "Ledoit-Wolf", "OAS", "Touloumis", "Schafer-Strimmer" ), lambda = NULL,...){

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

  X <- dat$X
  z = rep(0, ncol(X))
  names(z) = colnames(X)

  ids = ids[ids %in% colnames(X)]
  include = match(ids, colnames(X))

  if( length(include) < 1){
    return( NULL)
  }

  if( ncol(X) - length(include) < 1){
    return( NULL)
  }

  # impute z-statistic
  if( method == "decorrelate" ){
    res <- imputezDecorr(z, X, include, lambda = lambda,...)
  }else{
    if( is.null(lambda) ){
      lambda <- estimate_lambda(standardise(X), method)
    }
    res <- imputez(z, cor(X), include, lambda = lambda,...)
  }

  # set MAF
  res$maf <- colsums(X[,include,drop=FALSE]) / (2*nrow(X))
  res$maf <- pmin(res$maf, 1 - res$maf)
  res$nVariants <- ncol(X) - length(include)

  with(res, 
    data.frame(ID, r2.pred, lambda, maf, nVariants))
}


#' Evaluate imputation quality
#'
#' Evaluate 
#'
#' @param df \code{data.frame} with columns \code{ID},  \code{CHROM}, \code{POS}
#' @param gds \code{GenomicDataStream} of reference panel
#' @param window size of window in bp 
#' @param flankWidth additional window added to \code{region} 
#' @param ids variant IDs to evaluate imputation r2 score
#' @param method method used to estimate shrinkage parameter lambda.  default is \code{"decorrelate"}
#' @param lambda (default: NULL) value used to shrink correlation matrix. Only used if method is \code{"decorrelate"}
#' @param quiet suppress messages
#' @param ... additional arguments passed to \code{impute_region()}
#'
#' @details Implements method by Pasaniuc, et al. (2014).
#'
#' @references
#' Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., ... & Price, A. L. (2014). Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), 2906-2914.
#'
#' @return \code{tibble} storing imputed results:
#' \describe{
#'   \item{ID}{variant identifier}
#'   \item{r2.pred}{metric of accuracy of the imputed z-statistic based on its variance}
#'   \item{lambda}{shrinkage parameter}
#'   \item{maf}{minor allele frequency in reference panel}
#'   \item{nVariants}{number of variants used in imputation}
#' }
#'
#' @examples
#' library(GenomicDataStream)
#' library(tidyverse)
#' 
#' # VCF file for reference
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#' 
#' # initialize data stream
#' gds <- GenomicDataStream(file, "DS", initialize=TRUE)
#' 
#' # read file of variant locations
#' file <- system.file("extdata", "test.map", package = "GenomicDataStream")
#' df = read_tsv(file, show_col_types=FALSE)
#' 
#' # evaluate imputation r2 for these variants
#' ids = df$ID[c(1,9)]
#' 
#' # Given GenomicDataStream of reference panel,
#' # compute imputation r2 for variants in ids
#' run_imputability(df, gds, 10000, 1000, ids)
#
#' @seealso \code{run_imputez()}
#' @importFrom progress progress_bar
#' @importFrom GenomicDataStream setChunkSize
#' @importFrom dplyr bind_rows group_by slice_head
#' @export
run_imputability <- function( df, gds, window, flankWidth, ids, method = c("decorrelate", "Ledoit-Wolf", "OAS", "Touloumis", "Schafer-Strimmer"), lambda = NULL, quiet=FALSE,...){

  method <- match.arg(method)
  gds <- setChunkSize(gds, 1e9)

  regions <- get_analysis_windows( df, window)

  if (!quiet) {
    pb <- progress_bar$new(
      format = "  imputing [:bar] :percent eta: :eta",
      total = length(regions), clear = FALSE, width = 60)
  }

  ID <- NULL

  # impute in each region
  res <- lapply(regions, function(region){
    if (!quiet) pb$tick()
    imputability_region(gds, region, flankWidth, ids, method, lambda, ...)
  })
  res <- bind_rows(res) %>%
    group_by(ID) %>% 
    slice_head

  if (!quiet) {
    while (!pb$finished) pb$tick()
    pb$terminate()
  }

  res
}






