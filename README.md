
# Installation
```r
devtools::install_github("GabrielHoffman/imputez")
```


# Test imputation
NOTE: z-statistics must be input in order along chromosome

Read in z-statistics and LD information
```r
# devtools::install_github("GabrielHoffman/imputez")
library(imputez)

library(data.table)
library(Matrix)

# read z-statistics
file = "/sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V2/exampleForGabriel/GWAS_Zscores_corrected/Bellenguez_AD_chr22.tsv"
df_z_obs = fread(file)

# Specify locations of bim and ld files for each chromosome
# This example uses a data frame with 1 row
# for real analysis use rows for 1:22
df = data.frame(chrom = 22)
df$BIM = "/sc/arion/projects/roussp01a/sanan/230322_GWASimp/imputeZPipeline_V2/inter/230424/plinkStep2/1kg_chr22.bim"
df$LD = "/sc/arion/scratch/hoffmg01/chr22.ld.gz"

# read LD info
LDm = readLDMatrix( df )

# LDm = readRDS("/sc/arion/scratch/hoffmg01/LDm_chr22.RDS")

# store as RDS to reduce loading time in the future
# saveRDS(LDm, file="/sc/arion/scratch/hoffmg01/LDm_chr22.RDS")
```

Perform imputation on observed z-statistics
```r
# Create vector with z-statistics and variant names
# the z-statistics must be sorted by chromosome location
z = df_z_obs$Z
names(z) = df_z_obs$SNP

# Run z-statistic imputation on one chromosome
df_z = run_imputez(z, LDm[['22']]$dfld, 1:1000)
```

Plot results
```r
# Plot comparing observed and imputed values
library(ggplot2)
library(tidyverse)

df = data.frame(df_z, z.orig = z[df_z$ID]) 

r = with(df, cor(z.orig, z.stat))

ggplot(df, aes(z.orig, z.stat, color=r2.pred)) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) +
		geom_abline(color="red") + 
		xlab("Observed z") + 
		ylab("Imputed z") + 
		ggtitle("Compare observed and imputed") +
		scale_color_gradient(low="lightblue", high="black", limits=c(0,1)) +
		annotate(geom="text", x=-2, y=2, label=paste0("R=",format(r, digits=4)))
```

Additional plots
```r
plot(df_z$width, abs(z[df_z$ID]-df_z$z.stat))

plot(df_z$r2.pred, abs(z[df_z$ID]-df_z$z.stat))
```

Impute unobserved z-statistics
```r
# Impute unobserved z-statistics
idx = 50:255
z[idx] = NA
df_z = run_imputez(z, LDm[['22']]$dfld, idx)
```

