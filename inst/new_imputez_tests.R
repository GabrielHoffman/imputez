


ml plink2/b3.43 tabix
plink --bfile /sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V4/inter/230615_V2/plinkStep1/1kg_chr22 -recode vcf  --out /sc/arion/scratch/hoffmg01/test/1kg_chr22
bgzip 1kg_chr22.vcf
tabix -p vcf 1kg_chr22.vcf.gz


library(imputez)
library(tidyverse)
library(GenomicDataStream)

# read variant positions
file = "/sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V4/inter/230615_V2/plinkStep1/1kg_chr22.bim"
df_map = read_delim(file, col_names=FALSE)
colnames(df_map) = c("chrom", "ID", "map", "position", "REF_A1" ,"REF_A2")

# read z-statistics
file = "/sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V2/exampleForGabriel/GWAS_Zscores_corrected/Bellenguez_AD_chr22.tsv"
df_z_obs = read_tsv(file) %>%
			select(-LD_A1, -LD_A2)
colnames(df_z_obs)[1:2] = c("ID", "z")

# join
df = inner_join(df_z_obs, df_map, by="ID")

# if alleles are flipped, take negative z-stat
idx = with(df, GWAS_A1 == REF_A2 & GWAS_A2 == REF_A1)
tmp = df$GWAS_A1[idx]
df$GWAS_A1[idx] = df$GWAS_A2[idx]
df$GWAS_A2[idx] = tmp
df$z[idx] = df$z[idx]


# Run analysis with multiple methods
#-----------------------------------

file = "/sc/arion/scratch/hoffmg01/test/1kg_chr22.vcf.gz"
gds = GenomicDataStream(file, field="GT", region='22', MAF=0.05)

res = run_imputez(df, gds, 1000000, 250000)	

# test exclusion
idx = seq(1, nrow(df), by=10)
res = run_imputez(df[-idx,], gds, 1000000, 250000)	
res$method = "decorrelate"

# join imputed and observed t-statistics
df_res = df[idx,] %>%
			inner_join( res, by="ID")



# res = run_imputez(df[-idx,], gds, 1000000, 250000, lambda = .7)	





# TODO
# create GenomicDataStream::chunkSize() so not dependent on user
# allow region to be "" in GenomicDataStream



idx = seq(1, nrow(df), by=10)
methods = c("decorrelate", "Ledoit-Wolf", "OAS", "Touloumis", "Schafer-Strimmer" )

df_time = list()

res = lapply(methods, function(method){
	message(method)
	tm = system.time({
		res = run_imputez(df[-idx,], gds, 1000000, 250000, method=method)%>%
		mutate(method = method)
		})
	df_time[[method]] <<- tm
	res 
})
res = bind_rows(res)

df_res = res %>%
			inner_join( df[idx,], by="ID")



# Evaluate performance
#---------------------

rmse = function(x) sqrt(mean(x^2))


# rMSE
df_res %>%
	filter(maf > 0.05) %>%
	group_by(method) %>%
	summarize(rMSE = rmse(z.stat - z), 
			rMSE.mod = rmse(z.stat/se - z))

# imputez vs observed z-statistics
lim = range(c(df_res$z, df_res$z.stat))
fig = df_res %>%
	filter(maf > 0.05) %>%
		arrange(-r2.pred) %>%
		ggplot(aes(z, z.stat, color=r2.pred)) +
			geom_point() +
			theme_classic() +
			xlab("Observed z-statistic") +
			ylab("Imputed z-statistic") +
			geom_abline(slope=1, intercept=0) +
			scale_color_gradient(low="grey30", high="red", limits=c(0,1)) +
			coord_fixed(ratio = 1) +
			ylim(lim) + 
			xlim(lim) +
			facet_wrap( ~ method)
ggsave(fig, file="~/www/test.png")


# histogram of lambda values
fig = df_res %>%
		ggplot(aes(lambda)) +
			geom_histogram() +
			theme_classic() +
			theme(aspect.ratio=1) +
			xlab(bquote(lambda)) +
			facet_wrap( ~ method)
ggsave(fig, file="~/www/test.png")






res_eclairs = run_analysis(df$ID[idx], df, gds)
res_fixed = run_analysis(df$ID[idx], df, gds, lambda = 0.01)


png("~/www/test.png")
par(mfrow=c(2,2))
hist(res_eclairs$nVariants)
hist(res_eclairs$lambda)
dev.off()







plot_concordance = function(df_res){

	# Quantify concordance
	rsq = df_res  %>%
		summarize(Rsq = cor(z, z.stat)^2) %>%
		pull %>%
		format(digits=4)

	rmse = df_res  %>%
		summarize(Rsq = sqrt(mean((z - z.stat)^2))) %>%
		pull %>%
		format(digits=4)

	lim = range(c(df_res$z, df_res$z.stat))

	df_res %>%
		arrange(-r2.pred) %>%
		ggplot(aes(z, z.stat, color=r2.pred)) +
			geom_point() +
			theme_classic() +
			xlab("Observed z-statistic") +
			ylab("Imputed z-statistic") +
			geom_abline(slope=1, intercept=0) +
			scale_color_gradient(low="grey30", high="red", limits=c(.5,1)) +
			coord_fixed(ratio = 1) +
			ylim(lim) + xlim(lim) +
			annotate("text", x=lim[1]*.6, y=lim[2]*.9, label=bquote(R^2==.(rsq))) +
			annotate("text", x=lim[1]*.6, y=lim[2]*.7, label=bquote(rMSE==.(rmse)))
}




# joint imputed and observed t-statistics
df_res1 = res_eclairs %>%
					left_join(df, by="ID") %>%
					# mutate(z.stat = z.stat / sqrt(sigSq)) %>%
					filter(maf > 0.05)
fig1 = plot_concordance(df_res1)


# joint imputed and observed t-statistics
df_res2 = res_fixed %>%
					left_join(df, by="ID") %>%
					# mutate(z.stat = z.stat / sqrt(sigSq)) %>%
					filter(maf > 0.05)
fig2 = plot_concordance(df_res2)

fig = cowplot::plot_grid(fig1, fig2)
ggsave(fig, file="~/www/test.png", width=9)






png("~/www/test.png")
par(mfrow=c(2,2))
hist(df_res1$r2.pred, main="eclairs")
hist(df_res2$r2.pred, main="fixed")
dev.off()






with(df_res1, cor(r2.pred, (z.stat-z)^2, method="sp"))
with(df_res2, cor(r2.pred, (z.stat-z)^2, method="sp"))


with(df_res1, cor(r2.pred, (z.stat-z)^2))
with(df_res2, cor(r2.pred, (z.stat-z)^2))



fig1 = ggplot(df_res1, aes(r2.pred, (z.stat-z)^2)) +
	geom_point() +
	theme_classic() +
	theme(aspect.ratio=1) +
	# scale_y_log10() +
	ggtitle("eclairs")

fig2 = ggplot(df_res2, aes(r2.pred, (z.stat-z)^2)) +
	geom_point() +
	theme_classic() +
	theme(aspect.ratio=1)+
	# scale_y_log10() +
	ggtitle("fixed")

fig = cowplot::plot_grid(fig1, fig2, ncol=1)
ggsave(fig, file="~/www/test.png", width=4)





# df_res = df_res %>%
# 	filter(r2.pred > 0.8) 



library(GenomicDataStream)


file = "~/Downloads/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38.vcf.gz"
gds = GenomicDataStream(file, field="GT", initialize=TRUE, region='chr6', chunkSize=10000)


obj = getNextChunk(gds)

library(vcfppR)

file = "/Users/gabrielhoffman/Downloads/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38.vcf.gz"
res = vcftable(file, region="chr6")





df_res[with(df_res, order(abs(z-z.stat))),] %>% 
	tail %>%
	data.frame


i = which(df$ID == 'rs9612330')



imputez(z, cor(X), idx, lambda=0.1)

Rprof()
a = replicate(200, imputez.eclairs(z, X, idx))
res = summaryRprof()


microbenchmark(imputez.eclairs(z, X, idx, lambda=0.1), imputez(z, cor(X), idx, lambda=0.1))








library(decorrelate)
library(imputez)
library(mvtnorm)

rmse = function(x) sqrt(mean(x^2))

p = 500
Sigma = autocorr.mat(p, .4)
mu = seq(-sqrt(p), sqrt(p), length.out=p)
mu[] = 0
z = c(rmvnorm(1, mu, Sigma))
names(z) = 1:p

i = floor(p/2)
mu[i]
imputez(z, Sigma, i, lambda=0.0)
z[i]

df = lapply(seq(p), function(i){
	imputez(z, Sigma, i, lambda=0.0)
	})
df = do.call(rbind, df)

plot(z, df$z.stat.normalized)
points(z, df$z.stat, col="green3", pch=20)
abline(0,1, col="red")

rmse(z - df$z.stat)
rmse(z - df$z.stat.normalized)



stat1 = with(df, (z - df$z.stat)^2)
stat2 = with(df, (z - df$z.stat)^2 / (se^2))

mean(stat1)
mean(stat2)







