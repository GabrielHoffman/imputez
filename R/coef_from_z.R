
#' Compute coefficient estimate from z-statistic
#' 
#' Compute coefficient estimate from regression model and its standard error give z-statistic, sample size and standard deviation of response and covariate
#' 
#' @param z z-statistic
#' @param n sample size
#' @param sd_x standard deviation of covariate
#' @param sd_y standard deviation of response
#' @param phi case ratio in the response for logistic regression.  If specified, \code{sd_y} is set to \code{phi} and the coefficient is transformed to the logistic scale
#' @param p allele frequency of baseline allele
#' @param useMethod numeric.  1: first order method from Pirinen, et al. 2: second order method from Pirinen, et al. 3: method from Lloyd-Jones, et al.
#' 
#' @return \code{data.frame} storing \code{coef}, \code{se}, and \code{method}
#' 
#' @examples
#' # Linear regression
#' #------------------
#' 
#' # simulate data
#' set.seed(1)
#' n = 100
#' x = rnorm(n, 0, 7)
#' y = x*3 + rnorm(n)
#' data = data.frame(x, y)
#' 
#' # fit regression model
#' fit <- lm(y ~ x, data=data)
#' 
#' # get z-statistic
#' z = coef(summary(fit))[2,'t value']
#' 
#' # coef and se from regression model
#' coef(summary(fit))[2,-4]
#' 
#' # coef and se from summary statistics
#' coef_from_z(z, n, sd(x), sd(y))
#' 
#' # Logistic regression
#' #--------------------
#' 
#' # simulate data
#' n = 1000
#' p = .3
#' x = rbinom(n, 2, p)
#' eta = x*.1 
#' y = rbinom(n, size=1, prob=plogis(eta))
#' data = data.frame(x, y)
#' 
#' # fit regression model
#' fit <- glm(y  ~ x, data=data, family=binomial)
#' 
#' # get z-statistic
#' z = coef(summary(fit))[2,3]
#' 
#' # get case ratio 
#' phi = sum(y) / length(y)
#' 
#' # coef and se from regression model
#' coef(summary(fit))[2,]
#' 
#' # when phi is given, coef is transformed to logistic scale
#' coef_from_z(z, n, sd(x), phi=phi, p=p)
#' 
#' @details
#' Given a z-statistic, we want to obtain the coefficient value from a linear regression.  Adapting the approach from Zhu, et al (2016, Methods eqn 6), we estimate the coefficient as \deqn{\beta_{linear} = z * sd_y / (sd_x*sqrt(n + z^2)),} where \eqn{sd_x} is the standard deviation of the covariate and \eqn{sd_y} is the standard error of the response.  For a model with no covariates, this transformation gives the exact coefficient estimate.  With covariates, it is approximate.
#' 
#' The coeffient estimate from linear regression can be converted to the logistic scale using the approach of Pirinen, et al. (2013) according to \deqn{\beta_{logistic} = \beta_{linear} / (\phi(1-\phi)),} where \eqn{\phi} is the case ratio in the logistic regression.  This approximates the coefficient as if the model had been fit with logistic regression.  An alternative approximation is given by Lloyd-Jones, et al. (2018).
#' 
#' @references
#' Zhu, et al. (2016). Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nature Genetics. 48:481–487 \doi{10.1038/ng.3538}
#' 
#' Pirinen, Donnelly and Spencer. (2013). Efficient computation with a linear mixed model on large-scale data sets with applications to genetic studies. Ann. Appl. Stat. 7(1): 369-390 \doi{10.1214/12-AOAS586}
#'
#' Lloyd-Jones, Robinson, Yang and Visscher. (2018). Transformation of summary statistics from linear mixed model association on all-or-none traits to odds ratio. Genetics, 208(4), pp.1397-1408. \doi{10.1534/genetics.117.300360}
#' 
#' @export
coef_from_z <- function(z, n, sd_x, sd_y, phi, p, useMethod = 2){

	if( missing(sd_y) ){
		if(! missing(phi) ){
			sd_y <- phi
		}else{
			stop("Must specify sd_y")
		}
	}

	# convert z-statistic into coef from linear regression
	se <- sd_y / (sd_x*sqrt(n + z^2))
	beta <- z * se

	method = "linear"

	# if phi is given
	if( ! missing(phi) ){

		# convert from linear to logistic scale
		if( useMethod == 1 ){
			# first order from Pirinen, et al.
			beta <- beta / (phi*(1-phi))
		}else if( useMethod == 2){
			# second order from Pirinen, et al. 
			beta <- beta/((phi*(1-phi)) + 0.5*(1-2*phi) * (1 - 2*p)*beta)
		}else{
			# Lloyd-Jones, et al.
		  beta <- log(approx_OR(beta, p, phi))
		}

		# set se to give observed z-statistic	
		se <- beta / z

		method <- "logistic"
	}

	data.frame(coef = beta, se = se, method = method)
}




# Adapted from https://github.com/lukelloydjones/ORShiny/blob/master/shiny_lmor_func.R
approx_OR <- function(beta, p, k){
  # Odds ratio 3. A remnant.
  p_1 <- p  + (beta * p * (1 - p)) / k
  p_0 <- (p - k * p_1) / (1 - k)
  or.map.1  <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
  # -----------------------
  # Solve for the quadratic
  # -----------------------
  pp <- p
  c1 <- (1 - k) * (3 * beta - 2 * k * beta) + (beta * ((1 - k) ^ 2) * 
        (1 - 2 * k)) / k
  c2 <- (1 - k) * (1 - 4 * beta * pp) - (2 * beta * pp * (1 - k) * 
        (1 - 2 * k)) / k
  c3 <- (pp ^ 2) * beta * (1 - 2 * k) / k + pp * (k - 1 + beta)
  f0 <- c1 + c2 + c3
  # ------------------------------------------------------------------------
  # Calculate the upper and lower bounds for returning sensible ORs
  # ------------------------------------------------------------------------
  lb <- (p * k - p * (k ^ 2) - k * (1 - k)) / 
        ((1 - k) - 2 * p * (1 - k) + p ^ 2 - 2 * (p ^ 2) * k + p * k)
  ub <- (p * k - p * (k ^ 2)) / (p ^ 2 - 2 * (p ^ 2) * k + p * k)
  # Set up an empty array
  or.map.gld <- array(0, length(p))
  if (length(which(is.na(c3))) != 0 & length(which(is.na(f0))) != 0)
  {
    # --------------------
    # Odds ratio quadratic
    # --------------------  
    p.0.1 <- (-c2 + sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    p.0.2 <- (-c2 - sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    # There should be only one solution in (0, 1). Take the set that 
    # has all elements in (0, 1)
    if (all(p.0.1 > 0 & p.0.1 < 1))
    {
      p_0 <- p.0.1
    } else if (all(p.0.2 > 0 & p.0.2 < 1))
    {
      p_0 <- p.0.2
    }
    p_1 <- (p - (1 - k) *  p_0 ) / k  
    or.map.gld  <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
    or.map.gld[which(is.na(c3))] <- NA
    or.map.gld[which(is.na(f0))]  <- NA   
  } else if (!all(beta > lb & beta < ub))  
  {
    c1[which(beta < lb | beta> ub)] <- NA
    c2[which(beta < lb | beta> ub)] <- NA
    c3[which(beta < lb | beta> ub)] <- NA
    p.0.1 <- (-c2 + sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    p.0.2 <- (-c2 - sqrt(c2 ^ 2 - 4 * c1 * c3)) / (2 * c1)
    # There should be only one solution in (0, 1). 
    # Take the set that has all elements in (0, 1)
    if (!all(is.na(p.0.1)) & !all(is.na(p.0.1)))
    {
      if (all(p.0.1[-which(is.na(p.0.1))] > 0 & 
              p.0.1[-which(is.na(p.0.1))] < 1))
      {
        p_0 <- p.0.1
      } else if (all(p.0.2[-which(is.na(p.0.2))] > 0 & 
                     p.0.2[-which(is.na(p.0.2))] < 1))
      {
        p_0 <- p.0.2
      }
    } else
    {
      p_0 <- NA
    }
    p_1 <- (p - (1 - k) *  p_0 ) / k
    or.map.gld <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
  } else
  {
    # If there are no concerns with NAs and ORs on the boundary
    # calculate the full set of results
    or.map.gld <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
  }
  

  if( is.nan(or.map.gld) | is.na(or.map.gld) ){
		p_1 <- p + (beta * p * (1 - p)) / k
		p_0 <- (p - k * p_1) / (1 - k)
		or.map.1  <- (p_1 * (1 - p_0)) / (p_0 * (1 - p_1))
		or.map.gld <- log(or.map.1)
	}

  or.map.gld
}
   






