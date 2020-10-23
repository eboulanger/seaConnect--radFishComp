# Run RDA using previously compiled explanatory vars containing

# - environmental data: mean SST
# - habitat data: 4 habitat PC's
# - MEM: selected through OrdiStep
# - AEM: selected through Ordistep

# load libraries
library(adegenet)
library(psych)
library(reshape)
library(ape)
library(dplyr)
library(vegan)
library(car)

##### load data ####
# load genetic data
genepop.dip <- read.genepop("data/dip_neutral_7655.gen", ncode = 3L)
genepop.mul <- read.genepop("data/mul_neutral_2462.gen", ncode = 3L)

# load explanatory vars
expl.dip <- read.csv("output/rda_expl_vars_dip_neutral.csv",row.names= 1)
expl.mul <- read.csv("output/rda_expl_vars_mul_neutral.csv",row.names= 1)
# add SSS (wanneer tijd: proper in voorbereidende script doen)
SSS_dip <- read.csv("../00-Data/20-envdata-EB/dip_individual_envData.csv", row.names = 1, stringsAsFactors = F) %>% 
   select(c("Individual", "SSS_yearly_mean"))
SSS_mul <- read.csv("../00-Data/20-envdata-EB/mul_individual_envData.csv", row.names = 1, stringsAsFactors = F) %>% 
   select(c("Individual","SSS_yearly_mean"))
expl.dip <- left_join(expl.dip, SSS_dip)
expl.mul <- left_join(expl.mul, SSS_mul)

# make sure response and explanatory data are in same order
rownames(genepop.dip@tab) == expl.dip$Individual
rownames(genepop.mul@tab) == expl.mul$Individual

##### prepare response variables #####
## Calculate Euclidean distances from genepop object 

distgenEUCL.dip <- dist(genepop.dip, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
distgenEUCL.mul <- dist(genepop.mul, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

# Perform a Principal Coordinates Analysis (PCoA) on the Euclidean distances 
# Genetic distances are then in adequate multivariate space format for running the db-RDA.

pcoa.dip <- pcoa(distgenEUCL.dip)
pcoa.dip
pcoa.mul <- pcoa(distgenEUCL.mul)
pcoa.mul

# Extract Pcoa principal components, which will be the response variable in the db-RDA.

X.dip <- pcoa.dip$vectors
X.mul <- pcoa.mul$vectors

# Check out explanatory variables
Y.dip <- expl.dip[,6:58]
Y.mul <- expl.mul[,6:52]

pairs.panels(Y.dip) # too many variables to see anything

cor(Y.dip) %>% melt() %>% mutate(abs.value = abs(value)) %>% filter(abs.value > 0.6 & abs.value < 1) # 4 sets of variables are highly correlated: MEM4-SST_yearly_mean, MEM2-habitatPC4
cor(Y.dip) %>% melt() %>% filter(X1 == "SST_yearly_mean")

cor(Y.mul) %>% melt() %>% mutate(abs.value = abs(value)) %>% filter(abs.value > 0.6 & abs.value < 1) # 3 set of variables is highly correlated
cor(Y.mul) %>% melt() %>% filter(X1 == "SST_yearly_mean")

#### run preliminary rda ####
# perform a db-RDA global model including all explanatory variables.

rda1.dip <- rda(X.dip, Y.dip)
rda1.mul <- rda(X.mul, Y.mul)

#Looking at the Variance Inflation Factor (VIF) for multicolinearity within the model. 
# When VIF is greater than around 10, this is problematic.

vif(rda1.dip) # SST_yearly_mean is problematic
vif(rda1.mul) # many vars > 10

# for now: don't exclude vars.

# look at explained variance
RsquareAdj(rda1.dip) # $r.squared = 0.1071277, $adj.r.squared = 0.0101491
RsquareAdj(rda1.mul) # $r.squared = 0.09897339, $adj.r.squared = 0.01007209

# check out significance of model
anova(rda1.dip, perm = 999)
anova(rda1.mul, perm = 999)

#### iter 1: pre-select 4 variables per category ####
# take the first four dbMEMs and AEMs. order of variables stems for preparatory ordistep procedure

# for dip
colnames(Y.dip)
Y.dip.pre <- Y.dip %>% select("SST_yearly_mean", "SSS_yearly_mean",
                              "habitatPC1","habitatPC2","habitatPC3","habitatPC4",
                              "MEM6","MEM17","MEM7","MEM9",
                              "AEM23","AEM1","AEM35","AEM55")

# for mul: MEM3 and AEM1 (corr with SST) and AEM5 (least important from the two, based on the order of ordistep results)
colnames(Y.mul)
Y.mul.pre <- Y.mul %>% select("SST_yearly_mean", "SSS_yearly_mean",
                              "habitatPC1","habitatPC2","habitatPC3","habitatPC4",
                              "MEM5","MEM1","MEM13","MEM2",
                              "AEM7","AEM30","AEM17","AEM9")

# run rda
rda.dip.pre <- rda(X.dip, Y.dip.pre)
rda.mul.pre <- rda(X.mul, Y.mul.pre)

vif(rda.dip.pre) 
vif(rda.mul.pre) 

# look at explained variance
RsquareAdj(rda.dip.pre) 
RsquareAdj(rda.mul.pre) 

# check model significance
anova(rda.dip.pre, perm = 999) # P = 0.001
anova(rda.mul.pre, perm = 999) # P = 0.001

## Ordi2Step to select most informative variables
# start working form an empty model
rda.dip.0 <- rda(X.dip ~ 1, Y.dip.pre)
rda.mul.0 <- rda(X.mul ~ 1, Y.mul.pre)

# OrdiR2step will move towards the global model with all explanatory variables
rda.dip.G<- rda(X.dip ~ ., Y.dip.pre)
rda.mul.G<- rda(X.mul ~ ., Y.mul.pre)

# run selection of variables
Sel.dip <- ordiR2step(rda.dip.0, scope = formula(rda.dip.G), direction="both") 
Sel.mul <- ordiR2step(rda.mul.0, scope = formula(rda.mul.G), direction="both") 

# check out selected variables 
Sel.dip$anova
Sel.mul$anova

# check out spatial scale selected MEM
s.value(expl.dip[,4:5], expl.dip$MEM6)
s.value(expl.dip[,4:5], expl.dip$AEM23)

s.value(expl.mul[,4:5], expl.mul$MEM13)
s.value(expl.mul[,4:5], expl.mul$AEM9)


#### build final model 
# Build a model with the selected variables and visualize the results
Sel.dip.vars <-  rownames(Sel.dip$anova)[-nrow(Sel.dip$anova)] %>% sub("..", "",.)
Sel.mul.vars <-  rownames(Sel.mul$anova)[-nrow(Sel.mul$anova)] %>% sub("..", "",.)

Y.dip.sel <- Y.dip %>% select(all_of(Sel.dip.vars)) # fill in vars selected by ordistep
Y.mul.sel <- Y.mul %>% select(all_of(Sel.mul.vars)) 

# run final rda
rda.dip.S<- rda(X.dip~.,Y.dip.sel)
rda.mul.S<- rda(X.mul~.,Y.mul.sel)

summary(rda.dip.S)   
summary(rda.mul.S)   

RsquareAdj(rda.dip.S) 
RsquareAdj(rda.mul.S)

vif(rda.dip.S)
vif(rda.mul.S)

anova.dip.S         <- anova(rda.dip.S,permutations = 9999)
anova.dip.S.terms   <- anova(rda.dip.S,permutations = 9999, by = "terms")
anova.dip.S.axis    <- anova(rda.dip.S,permutations = 9999, by = "axis")
anova.dip.S.margin  <- anova(rda.dip.S,permutations = 9999, by = "margin")

anova.mul.S         <- anova(rda.mul.S,permutations = 9999)
anova.mul.S.terms   <- anova(rda.mul.S,permutations = 9999, by = "terms")
anova.mul.S.axis    <- anova(rda.mul.S,permutations = 999, by = "axis")
anova.mul.S.margin  <- anova(rda.mul.S,permutations = 9999, by = "margin")

# variance explained by axes
summary(rda.dip.S)$cont$importance[,1:5]
summary(rda.mul.S)$cont$importance[,1:5]


#### plotting ####
plot(rda.dip.S)
plot(rda.mul.S)

# add ecoregion colours
eco.dip <- expl.dip$Ecoregion_adj
eco.mul <- expl.mul$Ecoregion_adj
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecoregions_adj

plot(rda.dip.S, type="n", scaling=3, main = "Diplodus sargus | 7 570 neutral SNPs")
points(rda.dip.S, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda.dip.S, scaling=3, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.dip]) # the individuals
text(rda.dip.S,display="bp", col="black", cex=1, scaling=3)  # the predictors
text(x = -6, y = 3, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.dip.S)$r.squared, 5)))
text(x = -6, y = 2.5, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.dip.S)$adj.r.squared, 5)))
   # save as pdf
   #dev.print(pdf, 'figures/rda_dip_ntrl_1.pdf')


plot(rda.mul.S, type="n", main = "Mullus surmuletus | 13 548 neutral SNPs")
points(rda.mul.S, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.mul.S, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.mul]) # the individuals
text(rda.mul.S,display="bp", col="black", cex=1, scaling=3)                              # the predictors
text(x = -21, y = 9, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.mul.S)$r.squared, 5)))
text(x = -21, y = 8, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.mul.S)$adj.r.squared, 5)))
   # save as pdf
   #dev.print(pdf, 'figures/rda_mul_ntrl_1.pdf')

# add dapc colours
dapc.dip <- read.table("../02-DAPC/01-Diplodus/b-neutral/DAPCgrp_dip_neutral_K2.txt", header = T)
rownames(rda.dip.S$Ybar) == dapc.dip$INDIVIDUALS # check if individuals in same order
dapc.dip <- factor(dapc.dip$STRATA)

dapc.mul <- read.table("../02-DAPC/02-Mullus/b-neural/DAPCgrp_mul_neutral_K2.txt", header = T)
rownames(rda.mul.S$Ybar) == dapc.mul$INDIVIDUALS # check if individuals in same order
dapc.mul <- factor(dapc.mul$STRATA)

plot(rda.dip.S, type="n", scaling=1, main = "Diplodus sargus | 494 adaptive SNPs")
points(rda.dip.S, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda.dip.S, scaling=1, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[dapc.dip]) # the individuals
text(rda.dip.S,display="bp", col="black", cex=1, scaling=3)  # the predictors
text(x = -11, y = 2.4, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.dip.S)$r.squared, 5)))
text(x = -11, y = 2.1, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.dip.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(dapc.dip), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_dip_adpt_1.pdf')

plot(rda.mul.S, type="n", main = "Mullus surmuletus | x adaptive SNPs")
points(rda.mul.S, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.mul.S, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[dapc.mul]) # the individuals
text(rda.mul.S,display="bp", col="black", cex=1, scaling=3)                              # the predictors
text(x = -17, y = 7, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.mul.S)$r.squared, 5)))
text(x = -17, y = 6.5, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.mul.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(dapc.mul), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_mul_adpt_1.pdf')

#### iter 1b: partial RDA with vars iter1 ####
rda.part.dip <- rda(X.dip ~ SSS_yearly_mean +  SST_yearly_mean +
                       Condition(MEM6+MEM17+AEM23+MEM9+AEM1+AEM55+MEM7),
                    data = Y.dip.pre)
anova.part.dip <- anova(rda.part.dip,permutations = 9999)
anova.part.dip.terms <- anova(rda.part.dip,permutations = 9999, by = "terms")
anova.part.dip.axis <- anova(rda.part.dip,permutations = 9999, by = "axis")
anova(rda.part.dip,permutations = 9999)
RsquareAdj(rda.part.dip)
plot(rda.part.dip)

rda.part.mul <- rda(X.mul ~ SSS_yearly_mean + SST_yearly_mean +
                       Condition(MEM13+MEM5+MEM1+AEM9+AEM30+AEM17+MEM2+AEM7),
                    data = Y.mul.pre)
anova.part.mul <- anova(rda.part.dip,permutations = 9999)
anova.part.mul.terms <- anova(rda.part.mul,permutations = 9999, by = "terms")
anova.part.mul.axis <- anova(rda.part.mul,permutations = 999, by = "axis")

RsquareAdj(rda.part.mul)
plot(rda.part.mul)

#### pure ordi2step results #### 
vif(Sel.dip)
vif(Sel.mul)

eco.dip <- expl.dip$Ecoregion_adj
eco.mul <- expl.mul$Ecoregion_adj
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecoregions_adj

plot(Sel.dip, type="n", scaling=3, main = "Diplodus sargus | 7 570 neutral SNPs")
points(Sel.dip, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(Sel.dip, scaling=3, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.dip]) # the individuals
text(Sel.dip,display="bp", col="black", cex=1, scaling=3)  # the predictors
text(x = 9, y = 3, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(Sel.dip)$r.squared, 5)))
text(x = 9, y = 2.5, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(Sel.dip)$adj.r.squared, 5)))
legend("bottomright", legend=levels(eco.dip), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_dip_ntrl_2.pdf')


plot(Sel.mul, type="n", main = "Mullus surmuletus | 13 548 neutral SNPs")
points(Sel.mul, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(Sel.mul, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.mul]) # the individuals
text(Sel.mul,display="bp", col="black", cex=1, scaling=3)                              # the predictors
text(x = -19, y = 9, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(Sel.mul)$r.squared, 5)))
text(x = -19, y = 8, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(Sel.mul)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.mul), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_mul_ntrl_2.pdf')
#### iter 2: VIF on MEM & AEM ####

# run VIF fonction on MEM and AEM variables, then add SSS and habitat again
# select spatial variables
Y.dip.spatial <- select(Y.dip, starts_with(c("MEM","AEM")))
Y.mul.spatial <- select(Y.mul, starts_with(c("MEM","AEM")))
# calculate correlation matrix
cor.dip.spatial <- cor(Y.dip.spatial)
cor.mul.spatial <- cor(Y.mul.spatial)
# apply VIF function to stepwise remove vars
vif.dip.spatial <- vif_func(cor.dip.spatial,thresh=5,trace=T)
vif.mul.spatial <- vif_func(cor.mul.spatial, thresh = 5, trace = T)
# see which vars were removed
rm.dip.spatial <- setdiff(colnames(Y.dip.spatial), vif.dip.spatial)
rm.mul.spatial <- setdiff(colnames(Y.mul.spatial), vif.mul.spatial)
# remove these from expl var dataset
Y.dip.spVif <- select(Y.dip, -rm.dip.spatial)
Y.mul.spVif <- select(Y.mul, -all_of(rm.mul.spatial))

# run RDA and check out VIF (global model to use for Ordi2Step)
rda.spVif.dip.G <- rda(X.dip ~ ., Y.dip.spVif)
rda.spVif.mul.G <- rda(X.mul ~ ., Y.mul.spVif)

vif(rda.spVif.dip.G) # looks good!
vif(rda.spVif.mul.G) # also good!

RsquareAdj(rda.spVif.dip.G)
RsquareAdj(rda.spVif.mul.G)

# run Ordi2Step to select most informative variables
# start working form an empty model
rda.spVif.dip.0 <- rda(X.dip ~ 1, Y.dip.spVif)
rda.spVif.mul.0 <- rda(X.mul ~ 1, Y.mul.spVif)

# run selection of variables
Sel.spVif.dip <- ordiR2step(rda.spVif.dip.0, scope = formula(rda.spVif.dip.G), direction="both") 
Sel.spVif.mul <- ordiR2step(rda.spVif.mul.0, scope = formula(rda.spVif.mul.G), direction="both") 

# check out selected variables 
Sel.spVif.dip$anova
Sel.spVif.mul$anova

# check out spatial scale selected MEM
s.value(expl.dip[,4:5], expl.dip$MEM17)
s.value(expl.dip[,4:5], expl.dip$AEM3)
s.value(expl.dip[,4:5], expl.dip$MEM17)

s.value(expl.mul[,4:5], expl.mul$AEM9)
s.value(expl.mul[,4:5], expl.mul$MEM5)

#### build final model 
# Build a model with the selected variables and visualize the results
Sel.spVif.dip.vars <-  rownames(Sel.spVif.dip$anova)[-nrow(Sel.spVif.dip$anova)] %>% sub("..", "",.)
Sel.spVif.mul.vars <-  rownames(Sel.spVif.mul$anova)[-nrow(Sel.spVif.mul$anova)] %>% sub("..", "",.)

Y.dip.spVif.sel <- Y.dip %>% select(all_of(Sel.spVif.dip.vars)) # fill in vars selected by ordistep
Y.mul.spVif.sel <- Y.mul %>% select(all_of(Sel.spVif.mul.vars)) 

# run final rda
rda.spVif.dip.S<- capscale(X.dip~.,Y.dip.spVif.sel)
rda.spVif.mul.S<- capscale(X.mul~.,Y.mul.spVif.sel)

summary(rda.spVif.dip.S)   
summary(rda.spVif.mul.S)   

RsquareAdj(rda.spVif.dip.S)
RsquareAdj(rda.spVif.mul.S)

vif(rda.spVif.dip.S)
vif(rda.spVif.mul.S)

anova(rda.spVif.dip.S, permutations = 1000)
anova(rda.spVif.dip.S, permutations = 1000, by = "margin")

anova(rda.spVif.mul.S, permutations = 1000)
anova(rda.spVif.mul.S, permutations = 1000, by = "margin")
anova(rda.spVif.mul.S, permutations = 1000, by = "terms")

#anova(rda.spVif.dip.S, by = "axis")
#anova(rda.spVif.mul.S, by = "axis")


### plotting
# add ecoregion colours
eco.dip <- factor(expl.dip$Ecoregion_adj)
eco.mul <- factor(expl.mul$Ecoregion_adj)
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecoregions_adj

plot(rda.spVif.dip.S)
plot(rda.spVif.dip.S, type="n", main = "Diplodus sargus | 7 570 neutral SNPs")
points(rda.spVif.dip.S, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.spVif.dip.S, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.dip]) # the individuals
text(rda.spVif.dip.S,display="bp", col="black", cex=1)  # the predictors
text(x = -25, y = 4, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.spVif.dip.S)$r.squared, 5)))
text(x = -25, y = 2.9, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.spVif.dip.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.dip), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_dip_ntrl_spatialVIF_ordi2step.pdf')

plot(rda.spVif.mul.S)
plot(rda.spVif.mul.S, type="n", main = "Mullus surmuletus | 13 548 neutral SNPs")
points(rda.spVif.mul.S, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.spVif.mul.S, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.mul]) # the individuals
text(rda.spVif.mul.S,display="bp", col="black", cex=1)                              # the predictors
text(x = -19, y = 9, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.spVif.mul.S)$r.squared, 5)))
text(x = -19, y = 8.5, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.spVif.mul.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.mul), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
#text(rda.spVif.mul.S, display = "sites")
# save as pdf
#dev.print(pdf, 'figures/rda_mul_ntrl_spatialVIF_ordi2step.pdf')

# correlation between vars
pairs.panels(Y.dip.spVif.sel[,c(1,7)])
pairs.panels(Y.mul.spVif.sel[,c(1,21)])

#### iter2b: partial RDA with significant vars iter2 ####

# extract the spatial variables retained after vif & ordi2step
Y.dip.spVif.sel %>% select(starts_with(c("MEM", "AEM"))) %>% colnames() %>% paste(., collapse="+")
Y.mul.spVif.sel %>% select(starts_with(c("MEM", "AEM"))) %>% colnames() %>% paste(., collapse="+")
# use these as conditional variables for partial rda
rda.part.dip <- capscale(X.dip ~ SSS_yearly_mean + SST_yearly_mean + habitatPC4 + 
                           Condition(MEM17+MEM16+MEM9+MEM14+MEM10+AEM18+AEM3+AEM4+AEM10+AEM12+AEM16),
                         data = Y.dip)
rda.part.mul <- capscale(X.mul ~ SSS_yearly_mean + habitatPC3 + SST_yearly_mean + 
                           Condition(MEM5+MEM18+MEM7+MEM14+MEM2+MEM16+MEM15+MEM6+MEM13+MEM17+AEM8+AEM9+AEM2+AEM3+AEM12+AEM14+AEM17+AEM10+AEM4+AEM11+AEM13),
                         data = Y.mul)

plot(rda.part.dip)
plot(rda.part.mul)

RsquareAdj(rda.part.dip)
RsquareAdj(rda.part.mul)

anova(rda.part.dip, permutations = 1000)
anova(rda.part.mul, permutations = 1000)

anova(rda.part.dip, by = "margin", permutations = 1000)
anova(rda.part.mul, by = "margin", permutations = 1000)
anova(rda.part.mul, by = "terms", permutations = 1000)


### plotting
col.SSS <- colorRampPalette(c("red", "yellow"))
col.SST <- colorRampPalette(c())

plot(rda.part.dip)
plot(rda.part.dip, type="n", main = "Diplodus sargus | neutral SNPs")
points(rda.part.dip, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
# colour by salinity
points(rda.part.dip, display="sites", pch=21, cex=1.3, col="gray32", bg=col.SSS(length(Y.dip))) # the individuals
text(rda.part.dip,display="bp", col="black", cex=1)  # the predictors
text(x = -18, y = 2, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.part.dip)$r.squared, 5)))
text(x = -18, y = 1.4, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.part.dip)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.dip), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_dip_ntrl_spatialVIF_ordi2step_partial.pdf')

plot(rda.part.mul)
plot(rda.part.mul, type="n", main = "Mullus surmuletus | 1 558 adaptive SNPs")
points(rda.part.mul, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.part.mul, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.mul]) # the individuals
text(rda.part.mul,display="bp", col="black", cex=1)                              # the predictors
text(x = -16, y = 9, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.part.mul)$r.squared, 5)))
text(x = -16, y = 8.2, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.part.mul)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.mul), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_mul_adpt_spatialVIF_ordi2step_partial.pdf')


