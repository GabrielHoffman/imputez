



# 1000 Genomes and C4 Reference Panel
# Convert multi-allelic variants to new variants
#-----------------------------------------------
ml bcftools tabix parallel
OUT=/sc/arion/projects/roussp01a/gabriel/ref_panels/1kg
for CHR in $(seq 1 22)
do
FILE=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/VCF/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools norm -m - $FILE | bcftools view -O b -o $OUT/1kg_chr${CHR}_norm.bcf &
done

ls $OUT/*.bcf | parallel bcftools index

# Write MAP file of variant locations
cd /sc/arion/projects/roussp01a/gabriel/ref_panels/1kg/
for CHR in $(seq 1 22)
do
	echo -e "CHROM\tPOS\tID\tA1\tA2" > 1kg_chr${CHR}_norm.map
	bcftools view -s HG00096 1kg_chr${CHR}_norm.bcf | grep -v "#" | cut -f1-5 >> 1kg_chr${CHR}_norm.map
done




# C4
#-----
wget https://personal.broadinstitute.org/giulio/panels/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38.vcf.gz
bcftools norm -m -  MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38.vcf.gz | bgzip > MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38_norm.vcf.gz
tabix -p vcf MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38_norm.vcf.gz

FILE=MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38_norm.vcf.gz
echo -e "CHROM\tPOS\tID\tA1\tA2" > $(basename $FILE .vcf.gz).map
bcftools view $FILE | grep -v "#" | cut -f1-5 >> $(basename $FILE .vcf.gz).map



library(GenomicDataStream)
file = "MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38_norm.vcf.gz"
gds = GenomicDataStream(file, field="GT", region='chr6', init=TRUE)

dat = getNextChunk(gds)
i = which(nchar(dat$info$A2) > 1)

# rename C4 alleles
dat$info$ID[i] = paste0(dat$info$ID[i], "_", dat$info$A2[i])
colnames(dat$X)[i] = dat$info$ID[i]

dat$info[i,]
head(dat$X[,i])


FILE=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/VCF/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools view -r chr22:37636754-37636756 $FILE


library(imputez)
library(tidyverse)
library(GenomicDataStream)

# read variant positions
# file = "/sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V4/inter/230615_V2/plinkStep1/1kg_chr22.bim"
file = "/sc/arion/projects/roussp01a/gabriel/ref_panels/1kg/1kg_chr22_norm.map"
df_map = read_delim(file)
colnames(df_map)[4:5] = c("REF_A1", "REF_A2")

# read z-statistics
# file = "/sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V2/exampleForGabriel/GWAS_Zscores_corrected/Bellenguez_AD_chr22.tsv"
# df_z_obs = read_tsv(file) %>%
# 			select(-LD_A1, -LD_A2)
# colnames(df_z_obs)[1:2] = c("ID", "z")

# GRCh37
file = "/sc/arion/projects/roussp01a/gabriel/ref_panels/data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
df_z_obs = read_tsv(file, comment="##") %>%
		arrange(CHROM, POS) %>%
		mutate(z = BETA / SE)  %>%
		rename(GWAS_A1 = A1, GWAS_A2 = A2)


# join
df = df_z_obs %>%
	select(-CHROM, -POS) %>% 
	inner_join(df_map, by="ID")

# if alleles are flipped, take negative z-stat
idx = with(df, GWAS_A1 == REF_A2 & GWAS_A2 == REF_A1)
tmp = df$GWAS_A1[idx]
df$GWAS_A1[idx] = df$GWAS_A2[idx]
df$GWAS_A2[idx] = tmp
df$z[idx] = df$z[idx]


# file = "/sc/arion/scratch/hoffmg01/test/1kg_chr22.vcf.gz"
# region = "22:37636754-37636756"
# gds2 = GenomicDataStream(file, field="GT", init=TRUE, region=region)
# dat2 = getNextChunk(gds2)

# file = "/sc/arion/projects/roussp01a/gabriel/ref_panels/1kg/1kg_chr22_norm.bcf"
# gds = GenomicDataStream(file, field="GT", init=TRUE, samples=sampleIDs, region=region)
# dat = getNextChunk(gds)


# Run analysis with multiple methods
#-----------------------------------

# FAM file
file = "/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/PLINK/chr1.fam"
df_fam = read.table(file)
colnames(df_fam) = c("FID", "IID", "PID", "MID", "Sex", "Pop")

# Population assignments
file = "/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/release_2013502/integrated_call_samples.20130502.ALL.ped"
df_pop = read.table(file, header=TRUE, sep="\t")
df_pop = df_pop[df_pop$Relationship %in% c("father", "unrel","mother"),]
df_pop = df_pop[df_pop$Individual.ID %in% df_fam$IID,]

pops = c('GBR','IBS', 'TSI','CEU')#, "FIN")
sampleIDs = df_pop$Individual.ID[df_pop$Population %in% pops]
length(sampleIDs)



# sampleIDs = df_pop2$IID[df_pop2$Pop == "EUR"]
# sampleIDs = read.table("/sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V4/inter/230615_V2/plinkStep1/1kg_chr22.fam")$V1

file = "/sc/arion/projects/roussp01a/gabriel/ref_panels/1kg/1kg_chr22_norm.bcf"
gds = GenomicDataStream(file, field="GT", MAF=0.05, init=TRUE, samples=sampleIDs)

res = run_imputez(df, gds, 1000000, 250000)	

# test exclusion
idx = seq(1, nrow(df), by=10)
res = run_imputez(df[-idx,], gds, 1000000, 250000)	
res$method = "decorrelate"

# join imputed and observed t-statistics
df_res = df[idx,] %>%
			inner_join( res, by="ID")


# rs7284165 is 37636755 
region = "22:36208043-38253564"
a = impute_region(df[-idx,], gds, region, 250000, "Sch")

b = impute_region(df[-idx,], gds, region, 250000, "dec")


a %>%
	filter(ID == "rs7284165") %>% 
	data.frame

b %>%
	filter(ID == "rs7284165") %>% 
	data.frame






idx = seq(1, nrow(df), by=10)
methods = c("decorrelate")#, "Schafer-Strimmer" ) #"Ledoit-Wolf", "OAS", "Touloumis", 

df_time = list()

res = lapply(methods, function(method){
	message(method)
	tm = system.time({
		res = run_imputez(df[-idx,], gds, 1000000, 250000, method=method) %>%
		mutate(method = method)
		})
	df_time[[method]] <<- tm
	res 
})
res = bind_rows(res)

df_res = res %>%
			inner_join( df[idx,], by="ID")

df_res %>% 
	filter(ID == "rs7284165") %>% 
	data.frame


df_res %>%
	arrange(-abs(z-z.stat)) %>%
	head %>%
	data.frame

# Evaluate performance
#---------------------

colMethods = c("decorrelate" = "red",
              "GIW-EB (k=50)" = "#9f1214",
              "lambda = 0" = "#0000cd", 
              "Ledoit-Wolf" = "#FF7F00",    
              "OAS"= "#cccc29", 
              "Touloumis" = "#A65628", 
              "Schafer-Strimmer" = "#f569b3", 
              "Pseudoinverse" = "green3",
              "Baseline" = "grey50",
              "Oracle" = "black")


rmse = function(x) sqrt(mean(x^2))


# Time
fig = do.call(rbind, df_time) %>%
	as.data.frame %>%
	rownames_to_column('Method') %>%
	ggplot(aes(Method, elapsed, fill=Method)) +
		geom_bar(stat="identity") +
		theme_classic() +
		theme(aspect.ratio=1, legend.position="none") +
		scale_fill_manual(values = colMethods) +
		scale_y_continuous(limits=c(0, NA), expand=c(0,0)) +
		ylab("Wall time (seconds)") +
		coord_flip()
ggsave(fig, file="~/www/test.png")





# rMSE
fig = df_res %>%
	filter(maf > 0.05) %>%
	group_by(method) %>%
	summarize(rMSE = rmse(z.stat - z), 
			rMSE.mod = rmse(z.stat/se - z)) %>%
	ggplot(aes(method, rMSE, fill=method)) +
		geom_bar(stat="identity") +
		theme_classic() +
		theme(aspect.ratio=1, legend.position="none") +
		scale_fill_manual(values = colMethods) +
		scale_y_continuous(limits=c(0, NA), expand=c(0,0)) +
		ylab("Root mean squared error")
ggsave(fig, file="~/www/test.png")



# rMSE - stratified
fig = df_res %>%
	filter(maf > 0.05) %>%
	group_by(method, nVariants) %>%
	summarize(rMSE = rmse(z.stat - z), 
			rMSE.mod = rmse(z.stat/se - z)) %>%
	ggplot(aes(nVariants, rMSE, color=method)) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1, legend.position="none") +
		scale_color_manual(values = colMethods) +
		scale_y_continuous(limits=c(0, NA), expand=c(0,0)) +
		ylab("Root mean squared error")
ggsave(fig, file="~/www/test.png")


# rMSE - stratified by r2.pred
fig = df_res %>%
	# filter(maf > 0.05) %>%
	ggplot(aes(r2.pred, (z.stat - z)^2, color=method)) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1, legend.position="none") +
		scale_color_manual(values = colMethods) +
		scale_y_continuous(limits=c(0, NA), expand=c(0,0)) +
		ylab("Squared error")
ggsave(fig, file="~/www/test.png")



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
ggsave(fig, file="~/www/test.png", height=5, width=5)


fig = df_res %>%
	filter(maf > 0.05) %>%
	arrange(-r2.pred) %>%
	ggplot(aes(z.stat, se, color=r2.pred)) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1) +
		scale_color_gradient(low="grey30", high="red", limits=c(.8,1)) +
		facet_wrap( ~ method)
ggsave(fig, file="~/www/test.png")



# histogram of lambda values
fig = df_res %>%
		ggplot(aes(lambda)) +
			geom_histogram() +
			theme_classic() +
			theme(aspect.ratio=1) +
			xlab(bquote(lambda)) +
			facet_wrap( ~ method) +
			xlim(0, 1)
ggsave(fig, file="~/www/test.png")



# histogram of r2.pred values
fig = df_res %>%
		ggplot(aes(r2.pred)) +
			geom_histogram() +
			theme_classic() +
			theme(aspect.ratio=1) +
			xlab(bquote(r2.pred)) +
			facet_wrap( ~ method) +
			xlim(NA, 1)
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







