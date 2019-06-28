#*****************************************************************************
#****** Detecting multilocus adaptation using Redundancy Analysis (RDA) ******
#*****************************************************************************

### vignette by Brenna Forester
### https://popgen.nescent.org/2018-03-27_RDA_GEA.html

setwd("/Users/eboulanger-admin/Documents/project_SEACONNECT/seaConnect--radFishComp/03-RDA/01-Diplodus/")

# load packages
library(adegenet)
library(psych)
library(vegan)

##### ALL SNPs #####

### load and prepare data ----

# read in genetic data and extract minor allele frequency matrix
diplodus_all <- read.PLINK("a-all/diplodus.raw")
dip_m <-as.matrix(diplodus_all)

diplodus_adpt <- read.PLINK("c-adaptive/diplodus_adaptive.raw")
dip_m <-as.matrix(diplodus_adpt)


# explore data
  # jpeg(filename = "glPlot_dip_all.jpeg", res = 200, width = 1200, height = 1200)
#glPlot(diplodus_all, posi="topleft") # allows to visualize missing data
  # dev.off()

# RDA requires complete data frames
sum(is.na(dip_m)) # 63 835 NA's in the matrix
# impute missing values using the most common genotype at each SNP across all individuals
dip.imp <- apply(dip_m, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(dip.imp)) # No NAs

# load environmental data
env <- read.table("../../00-envdata-EB/Diplodus_envData.txt", sep = "\t", header = T, stringsAsFactors = F)
str(env)
env$Ecoregion <- as.factor(env$Ecoregion)
env$Ecoregion_adj <- as.factor(env$Ecoregion_adj)

# Confirm that genotypes and environmental data are in the same order
identical(rownames(dip.imp), env$ID) # FALSE
temp <- cbind(rownames(dip.imp), env$ID) # order is indeed different
rownames(env) <- env$ID
env <- env[rownames(dip.imp),]
identical(rownames(dip.imp), rownames(env)) #TRUE

# check correlation between predictors
  # jpeg(filename = "dip_RDA_pred_corr4.jpeg", res = 200, width = 1200, height = 1200)
pairs.panels(env[,5:20], scale=T)
  # dev.off()
  # SST and SSS are highly correlated amongst eachother, of course, but only 60% between eachother -> ok for RDA
  # keep only max of each for now
  # also keep lat and lon -> test for autocorrelation / geographic distance?
  # PAR highly correlated with latitude and temperature
pred <- subset(env, select=c(lat, lon, SST_yearly_max, SSS_yearly_max, SST_april_mean, SSS_april_mean, mean_bat))
names(pred) <- c("lat", "lon","SST_yearly", "SSS_yearly","SST_april", "SSS_april", "bathy")
pairs.panels(pred, scale=T)


### run the RDA ----

dip.rda <- rda(dip.imp ~ ., data=pred, scale=T)
dip.rda

# calculate adjusted Rsquared
RsquareAdj(dip.rda)
 # $adj.r.squared
 # only SSS + SST: 0.004287249 -> very very low? 0,4% (wolf vignette: 5%)
 # include lat + lon: 0.01090912
 # include april + bathy: 0.01381814
# eigenvalues constrained axes
summary(eigenvals(dip.rda, model = "constrained"))
screeplot(dip.rda) # first axis explains most, rest all small portions

# check RDA model for significance
# null-hypothesis: no linear relationship exists between SNP data and environmental predictors
Sys.time()
signif.full <- anova.cca(dip.rda, parallel=getOption("mc.cores")) # default is permutation=999
Sys.time()
signif.full
  # SSS + SST :

  #  Permutation test for rda under reduced model
  #  Permutation: free
  #  Number of permutations: 999
  #  
  #  Model: rda(formula = dip.imp ~ SST_yearly_max + SSS_yearly_max, data = pred, scale = T)
  #  Df Variance      F Pr(>F)    
  #  Model      2      146 1.6372  0.001 ***
  #    Residual 294    13111 

  # include lat + lon : 

  #  Permutation test for rda under reduced model
  #  Permutation: free
  #  Number of permutations: 999
  #  
  #  Model: rda(formula = dip.imp ~ lat + lon + SST + SSS, data = pred, scale = T)
  #  Df Variance      F Pr(>F)    
  #  Model      4    321.8 1.8162  0.001 ***
  #    Residual 292  12935.2     

# check each constrained axis for significance
# run during lunch because can take a long time 
Sys.time() # SSS + SST: "2019-06-19 12:31:59 CEST"
           # + lat + lon: 
#signif.axis <- anova.cca(dip.rda, by="axis", parallel=getOption("mc.cores")) 
Sys.time() # SSS + SST:  "2019-06-19 12:53:06 CEST"
           # + lat + lon: 
signif.axis

  #  Permutation test for rda under reduced model
  #  Forward tests for axes
  #  Permutation: free
  #  Number of permutations: 999
  #  
  #  Model: rda(formula = dip.imp ~ SST_yearly_max + SSS_yearly_max, data = pred, scale = T)
  #               Df Variance      F Pr(>F)    
  #    RDA1       1     98.4 2.2062  0.001 ***
  #    RDA2       1     47.6 1.0683  0.013 *  
  #    Residual 294  13111.0       



### plot the RDA ----

# base 
plot(dip.rda, scaling=3)          # default is axes 1 and 2
plot(dip.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3

# more informative plot
#levels(env$Ecoregion_adj) <- c("Alboran Sea", "Western Mediterranean", "Central Mediterranean","Adriatic Sea","Aegean Sea","Levantine Sea")
#env$Ecoregion_adj <- factor(env$Ecoregion_adj, levels=c("Alboran Sea", "Western Mediterranean", "Central Mediterranean","Adriatic Sea","Aegean Sea","Levantine Sea")) # set levels order
eco <- env$Ecoregion_adj
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecoregions_adj
#bg8 <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a') # 8 nice colors for our ecoregions

# axes 1 & 2
  # jpeg(filename = "RDAplot_dip_all_3.jpeg", res = 200, width = 1200, height = 1200)
plot(dip.rda, type="n", scaling=3, xlim = c(-5, 4))
points(dip.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(dip.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the individuals
text(dip.rda, scaling=3, display="bp", col="black", cex=1)    #  col="#0868ac"                       # the predictors
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  # dev.off()

### Identify candidate SNPs ----
### involved in local adaptation 

load.rda <- scores(dip.rda, choices=c(1:2), display="species")  # Species scores for the first two constrained axes

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 231 candidate SNPs SSS + SST, 325 SNPs + lat + lon
cand2 <- outliers(load.rda[,2],3) # 49 candidate SNPs / 48

ncand <- length(cand1) + length(cand2)
ncand # 280 / 373 / 387

# create dataframe with results
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with the eight environmental predictors:
foo <- matrix(nrow=(ncand), ncol=7)  # 7 columns for 7 predictors
colnames(foo) <- c("lat", "lon","SST_yearly", "SSS_yearly","SST_april", "SSS_april", "bathy")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- dip.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

# look for duplicate candidates
length(cand$snp[duplicated(cand$snp)])  # 0 duplicate detections

# see which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,11] <- names(which.max(abs(bar[4:10]))) # gives the variable
  cand[i,12] <- max(abs(bar[4:7]))              # gives the correlation
}

colnames(cand)[11] <- "predictor"
colnames(cand)[12] <- "correlation"
table(cand$predictor) 

### plot candidate SNPs ----

# plot these highly correlated SNPs in the ordination space
sel <- cand$snp
envp <- cand$predictor
envp[envp=="lat"] <- '#1f78b4'
envp[envp=="lon"] <- '#a6cee3'
envp[envp=="SST_yearly"] <- '#6a3d9a'
envp[envp=="SSS_yearly"] <- '#e31a1c'
envp[envp=="SST_april"] <- '#33a02c'
envp[envp=="SSS_april"] <- '#ffff33'
envp[envp=="bathy"] <- '#fb9a99'

# color by predictor:
col.pred <- rownames(dip.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- envp[i]
}

col.pred[grep("scaffold",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a')

# plot the SNPs
# axes 1 & 2
  # jpeg(filename = "RDAplot_dip_SNP_3.jpeg", res = 200, width = 1200, height = 1200)
plot(dip.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(dip.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(dip.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(dip.rda, scaling=3, display="bp", col="black", cex=1)
legend("bottomleft", legend=c("lat", "lon","SST_yearly", "SSS_yearly","SST_april", "SSS_april", "bathy"), 
       bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  # dev.off()

