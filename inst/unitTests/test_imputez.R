

test_imputez = function(){

	library(GenomicDataStream)
	library(mvtnorm)
	library(imputez)
	library(RUnit)

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat = getNextChunk(gds)

	C = cor(dat$X)
	z = c(rmvnorm(1, rep(0, 10), C))
	names(z) = colnames(dat$X)

	res = imputez(z, cor(dat$X), 2:3, lambda=0.1)
	res2 = imputezDecorr(z, dat$X, 2:3, lambda=0.1)

	checkEqualsNumeric(res[,-1], res2[,2:ncol(res)])
}

test_coef_from_z = function(){

	# simulate data
	set.seed(1)
	n = 1000
	x = rnorm(n, 0, 7)
	y = x*3 + rnorm(n)
	data = data.frame(x, y)

	# fit regression model
	fit <- lm(y ~ x, data=data)

	# get z-statistic
	z = coef(summary(fit))[2,'t value']

	# coef and se from regression model
	res1 = coef(summary(fit))[2,-4]

	# coef and se from summary statistics
	res2 = coef_from_z(z, n, sd(x), sd(y))

	checkEqualsNumeric(res1[1:2], res2[1:2], tol=1e-5)
	checkEqualsNumeric(z, with(res2, coef / se))
}