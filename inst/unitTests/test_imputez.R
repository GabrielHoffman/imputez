

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