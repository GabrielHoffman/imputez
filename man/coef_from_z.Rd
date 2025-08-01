% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef_from_z.R
\name{coef_from_z}
\alias{coef_from_z}
\title{Compute coefficient estimate from z-statistic}
\usage{
coef_from_z(z, n, sd_x, sd_y, phi)
}
\arguments{
\item{z}{z-statistic}

\item{n}{sample size}

\item{sd_x}{standard deviation of covariate}

\item{sd_y}{standard deviation of response}

\item{phi}{case ratio in the response for logistic regression.  If specified, \code{sd_y} is set to \code{phi} and the coefficient is transformed to the logistic scale}
}
\value{
\code{data.frame} storing \code{coef}, \code{se}, and \code{method}
}
\description{
Compute coefficient estimate from regression model and its standard error give z-statistic, sample size and standard deviation of response and covariate
}
\details{
Given a z-statistic, we want to obtain the coefficient value from a linear regression.  Adapting the approach from Zhu, et al (2016, Methods eqn 6), we estimate the coefficient as \deqn{\beta_{linear} = z * sd_y / (sd_x*sqrt(n + z^2)),} where \eqn{sd_x} is the standard deviation of the covariate and \eqn{sd_y} is the standard error of the response.  For a model with no covariates, this transformation gives the exact coefficient estimate.  With covariates, it is approximate.

The coeffient estimate from linear regression can be converted to the logistic scale using the approach of Pirinen, et al. (2013) according to \deqn{\beta_{logistic} = \beta_{linear} / (\phi(1-\phi)),} where \eqn{\phi} is the case ratio in the logistic regression.  This approximates the coefficient as if the model had been fit with logistic regression.
}
\examples{
# Linear regression
#------------------

# simulate data
set.seed(1)
n = 100
x = rnorm(n, 0, 7)
y = x*3 + rnorm(n)
data = data.frame(x, y)

# fit regression model
fit <- lm(y ~ x, data=data)

# get z-statistic
z = coef(summary(fit))[2,'t value']

# coef and se from regression model
coef(summary(fit))[2,-4]

# coef and se from summary statistics
coef_from_z(z, n, sd(x), sd(y))

# Logistic regression
#--------------------

# simulate data
n = 1000
p = .3
x = rbinom(n, 2, p)
eta = x*.1 
y = rbinom(n, size=1, prob=plogis(eta))
data = data.frame(x, y)

# fit regression model
fit <- glm(y  ~ x, data=data, family=binomial)

# get z-statistic
z = coef(summary(fit))[2,3]

# get case ratio 
phi = sum(y) / length(y)

# coef and se from regression model
coef(summary(fit))[2,]

# when phi is given, coef is transformed to logistic scale
coef_from_z(z, n, sd(x), phi=phi)

}
\references{
Zhu, et al. (2016). Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nature Genetics. 48:481–487 \doi{10.1038/ng.3538}

Pirinen, Donnelly and Spencer. (2013). Efficient computation with a linear mixed model on large-scale data sets with applications to genetic studies. Ann. Appl. Stat. 7(1): 369-390 \doi{10.1214/12-AOAS586}
}
