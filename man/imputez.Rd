% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputez.R
\name{imputez}
\alias{imputez}
\title{Impute z-score based correlations}
\usage{
imputez(z, Sigma, i, lambda = 0.1)
}
\arguments{
\item{z}{vector of observed z-scores}

\item{Sigma}{matrix of correlation between z-scores}

\item{i}{index of z-scores to impute}

\item{lambda}{value used to shrink correlation matrix}
}
\value{
\code{data.frame} storing:
\describe{
  \item{ID}{variant identifier}
  \item{z.stat}{imputed z-statistic}
  \item{se}{standard error of imputed z-statistic}
  \item{r2.pred}{metric of accuracy of the imputed z-statistic based on its variance}
  \item{lambda}{shrinkage parameter}
}
}
\description{
Impute z-score from a missing test based on z-scores for other tests, and the correlation matrix between z-scores
}
\details{
Implements method by Pasaniuc, et al. (2014).
}
\references{
Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., ... & Price, A. L. (2014). Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), 2906-2914.
}
