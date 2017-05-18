sessionInfo()

#################################################################################################################################################################################################
#  																																									                            #
#R version 3.2.3 (2015-12-10)                                                                                                                                                                   #
#Platform: x86_64-w64-mingw32/x64 (64-bit)                                                                                                                                                      #
#Running under: Windows 10 x64 (build 10240)                                                                                                                                                    #
#                                                                                                                                                                                               #
#locale:                                                                                                                                                                                        #
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C      LC_TIME=English_United States.1252                  #
#                                                                                                                                                                                               #
#attached base packages:                                                                                                                                                                        #
#[1] grid      stats     graphics  grDevices utils     datasets  methods   base                                                                                                                 #
#                                                                                                                                                                                               #
#other attached packages:                                                                                                                                                                       #
# [1] phylotools_0.1.2 fields_8.3-6     spam_1.3-0       spaa_0.2.1       picante_1.6-2    nlme_3.1-124     vegan_2.3-3      lattice_0.20-33  permute_0.9-0    seqRFLP_1.0.1    Cairo_1.5-9     #
#[12] phytools_0.5-10  maps_3.0.2       ape_3.4                                                                                                                                                 #
#                                                                                                                                                                                               #
#loaded via a namespace (and not attached):                                                                                                                                                     #
# [1] igraph_1.0.1            cluster_2.0.3           magrittr_1.5       splines_3.2.3      nnls_1.4          MASS_7.3-45             mnormt_1.5-3    scatterplot3d_0.3-36                      #
# [9] quadprog_1.5-5          tools_3.2.3             parallel_3.2.3     msm_1.6            mgcv_1.8-11       clusterGeneration_1.3.4 phangorn_2.0.2  plotrix_3.6-1                             #
#[17] survival_2.38-3         numDeriv_2014.2-1       Matrix_1.2-3       animation_2.4      expm_0.999-0      mvtnorm_1.0-5                                                                     #
#                                                                                                                                                                                               #
#################################################################################################################################################################################################

library(ape)
library(phytools)
library(phylotools)
library(geiger)

tcaf <- read.nexus("caffeine.nxs")
tabl <- as.matrix(read.table("caffeine.table", header=TRUE, sep="\t", row.names=1))

new.tip.labels <- read.table("tip.translation.table", header = FALSE, sep="\t")
caftree <- sub.tip.label(caftree, new.tip.labels)

caf <- tabl[,1]
name.check(tcaf,caf)

BM <- fitContinuous(tcaf, caf, model="BM")
OU <- fitContinuous(tcaf, caf, model="OU")
LA <- fitContinuous(tcaf, caf, model="lambda")
KA <- fitContinuous(tcaf, caf, model="kappa")
DE <- fitContinuous(tcaf, caf, model="delta")
EB <- fitContinuous(tcaf, caf, model="EB")
WH <- fitContinuous(tcaf, caf, model="white")
TR <- fitContinuous(tcaf, caf, model="trend")

## WHICH EVOLUTIONARY MODEL IS BETTER? Compare Likelihoods using LRT and compare AIC

df <- as.matrix(data.frame(mods=c("BM", "OU", "LA", "KA", "DE", "EB", "WH", "TR"), lnL="NA", aic="NA"))
mdf <- df
mdf[1,2] <- BM$opt$lnL
mdf[2,2] <- OU$opt$lnL
mdf[3,2] <- LA$opt$lnL
mdf[4,2] <- KA$opt$lnL
mdf[5,2] <- DE$opt$lnL
mdf[6,2] <- EB$opt$lnL
mdf[7,2] <- WH$opt$lnL
mdf[8,2] <- TR$opt$lnL
mdf[1,3] <- BM$opt$aic
mdf[2,3] <- OU$opt$aic
mdf[3,3] <- LA$opt$aic
mdf[4,3] <- KA$opt$aic
mdf[5,3] <- DE$opt$aic
mdf[6,3] <- EB$opt$aic
mdf[7,3] <- WH$opt$aic
mdf[8,3] <- TR$opt$aic

mdf <- apply(mdf[,2:3],2,as.numeric)
row.names(mdf)=df[,1]

> mdf
          lnL        aic
BM   3.381097 -2.7621933
OU -16.789586 39.5791718
LA   3.381097 -0.7621933
KA   3.381097 -0.7621933
DE   3.509125 -1.0182497
EB   3.381097 -0.7621933
WH   3.364309 -2.7286170
TR   3.403729 -0.8074585

