

# confirm that logistic approximation in coef_from_z works
res = lapply(seq(1000), function(i){
	# simulate data
	n = 1000
	p = .3
	x = rbinom(n, 2, p)
	eta = x*0  + 2
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

	res1 = coef_from_z(z, n, sd(x), phi=phi, p=p)
	res2 = coef_from_z(z, 4*n*phi*(1-phi), sd(x), sd(y), p=p)

	data.frame(coef_true = coef(fit)[2],
		coef1 = res1$coef, 
		coef2 = res2$coef)
}) %>%
	bind_rows


coef(lm(coef_true ~ 0 + coef1, res))

par(mfrow=c(1,2))
with(res, plot(coef_true, coef1))	
abline(0,1, col="red")

with(res, plot(coef_true, coef2))	
abline(0,1, col="red")



