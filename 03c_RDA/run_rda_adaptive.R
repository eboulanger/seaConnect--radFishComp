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

select <- dplyr::select

##### load data ####
# load genetic data
genepop.dip <- read.genepop("data/dip_adaptive_413.gen", ncode = 3L)
genepop.mul <- read.genepop("data/mul_adaptive_291.gen", ncode = 3L)

# load explanatory vars
expl.dip <- read.csv("output/rda_expl_vars_dip_adaptive.csv",row.names= 1, stringsAsFactors = F)
expl.mul <- read.csv("output/rda_expl_vars_mul_adaptive.csv",row.names= 1, stringsAsFactors = F)
# add SSS (wanneer tijd: proper in voorbereidende script doen)
# SSS_dip <- read.csv("../00-Data/20-envdata-EB/dip_individual_envData.csv", row.names = 1, stringsAsFactors = F) %>% 
#   select(c("Individual", "SSS_yearly_mean"))
# SSS_mul <- read.csv("../00-Data/20-envdata-EB/mul_individual_envData.csv", row.names = 1, stringsAsFactors = F) %>% 
#   select(c("Individual","SSS_yearly_mean"))

# add seasonal temp & chlorophyll a
env.dip <- read.csv("../00-Data/20-envdata-EB/dip_individual_envData.csv", stringsAsFactors = F)
env.mul <- read.csv("../00-Data/20-envdata-EB/mul_individual_envData.csv", stringsAsFactors = F)

# add chlorophyll a

expl.dip <- left_join(expl.dip, select(env.dip, Individual, SSS_yearly_mean, SST_winter_mean, SST_summer_mean, chloSurfaceMean))
expl.mul <- left_join(expl.mul, select(env.mul, Individual, SSS_yearly_mean, SST_winter_mean, SST_summer_mean, chloSurfaceMean))
# monthly
expl.dip.monthly <- left_join(expl.dip, sst_det_dip)
expl.mul.monthly <- left_join(expl.mul, sst_det_mul)

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
Y.dip <- expl.dip[,6:46]
Y.mul <- expl.mul[,6:56]

pairs.panels(Y.dip) # too many variables to see anything

cor(Y.dip) %>% melt() %>% mutate(abs.value = abs(value)) %>% filter(abs.value > 0.6 & abs.value < 1) # 4 sets of variables are highly correlated: MEM4-SST_yearly_mean, MEM2-habitatPC4
cor(Y.dip) %>% melt() %>% filter(X1 == "MEM1")

#                          value   abs.value
#  MEM4 SST_yearly_mean  0.7271697 0.7271697
#  AEM1 SST_yearly_mean -0.6724576 0.6724576
#  MEM2      habitatPC4  0.6111314 0.6111314
# AEM18            MEM7 -0.7162332 0.7162332

cor(Y.mul) %>% melt() %>% mutate(abs.value = abs(value)) %>% filter(abs.value > 0.5 & abs.value < 1) # 3 set of variables is highly correlated
cor(Y.mul) %>% melt() %>% filter(X1 == "SSS_yearly_mean")

#                         value   abs.value
# MEM3 SST_yearly_mean -0.7546721 0.7546721
# AEM1 SST_yearly_mean -0.7124675 0.7124675
# AEM7            MEM5  0.7026138 0.7026138

#### run preliminary rda #####
# perform a db-RDA global model including all explanatory variables.

rda1.dip <- rda(X.dip, Y.dip)
rda1.mul <- rda(X.mul, Y.mul)

#Looking at the Variance Inflation Factor (VIF) for multicolinearity within the model. 
# When VIF is greater than around 10, this is problematic.

vif(rda1.dip) # SST_yearly_mean is problematic, AEM1 limit
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
Y.dip.pre <- Y.dip %>% select("SSS_yearly_mean",
                              "SST_winter_mean", "SST_summer_mean",
                              "chloSurfaceMean",
                            "habitatPC1","habitatPC2","habitatPC3","habitatPC4",
                            "MEM6","MEM16","MEM7","MEM9",
                            "AEM23","AEM35","AEM51","AEM52")
                            
# for mul: MEM3 and AEM1 (corr with SST) and AEM5 (least important from the two, based on the order of ordistep results)
colnames(Y.mul)
Y.mul.pre <- Y.mul %>% select("SSS_yearly_mean",
                              "SST_winter_mean", "SST_summer_mean",
                              "chloSurfaceMean",
                              "habitatPC1","habitatPC2","habitatPC3","habitatPC4",
                            "MEM5","MEM1","MEM12","MEM17",
                            "AEM41","AEM40","AEM17","AEM3")

# run rda
rda.dip.pre <- rda(X.dip, Y.dip.pre)
rda.mul.pre <- rda(X.mul, Y.mul.pre)

vif(rda.dip.pre) 
vif(rda.mul.pre) 

# look at explained variance
RsquareAdj(rda.dip.pre) # r.squared = 0.1071277,  adj.r.squared = 0.0101491 | r.squared = 0.1965732,  adj.r.squared = 0.1566867
RsquareAdj(rda.mul.pre) # r.squared = 0.09897339, adj.r.squared = 0.01007209 |r.squared = 0.1280608,  adj.r.squared = 0.1010538

# check model significance
anova(rda.dip.pre, perm = 9999) # P = 0.001
anova(rda.mul.pre, perm = 9999) # P = 0.001

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
s.value(expl.dip[,4:5], expl.dip$habitatPC2)

s.value(expl.mul[,4:5], expl.mul$MEM5)
s.value(expl.mul[,4:5], expl.mul$AEM40)
s.value(expl.mul[,4:5], expl.mul$habitatPC2)
s.value(expl.mul[,4:5], expl.mul$habitatPC4)

#### build final model 
# Build a model with the selected variables and visualize the results
Sel.dip.vars <-  rownames(Sel.dip$anova)[-nrow(Sel.dip$anova)] %>% sub("..", "",.)
          # [1] "MEM6"            "AEM23"           "MEM9"           
          # [4] "AEM52"           "AEM51"           "AEM35"          
          # [7] "habitatPC2"      "MEM7"            "SST_winter_mean"
          # [10] "SSS_yearly_mean"
Sel.mul.vars <-  rownames(Sel.mul$anova)[-nrow(Sel.mul$anova)] %>% sub("..", "",.)
          # [1] "SSS_yearly_mean" "MEM5"            "AEM40"          
          # [4] "MEM12"           "MEM1"            "chloSurfaceMean"
          # [7] "AEM3"            "MEM17"           "AEM41"          
          # [10] "habitatPC2"      "habitatPC4"      "habitatPC3"     
          # [13] "SST_summer_mean"

Y.dip.sel <- Y.dip %>% select(all_of(Sel.dip.vars)) # fill in vars selected by ordistep
Y.mul.sel <- Y.mul %>% select(all_of(Sel.mul.vars)) 

# run final rda
rda.dip.S<- rda(X.dip~.,Y.dip.sel)
rda.mul.S<- rda(X.mul~.,Y.mul.sel)

    # export rda
    #saveRDS(rda.dip.S, file = "output/rda.4sel.dip.adpt.rds")
    #saveRDS(rda.mul.S, file = "output/rda.4sel.mul.adpt.rds")

summary(rda.dip.S)   
summary(rda.mul.S)   

RsquareAdj(rda.dip.S) # 0.1551034
RsquareAdj(rda.mul.S)

vif(rda.dip.S)
vif(rda.mul.S)

anova.dip.S         <- anova(rda.dip.S,permutations = 9999)
anova.dip.S.terms   <- anova(rda.dip.S,permutations = 9999, by = "terms")
anova.dip.S.axis    <- anova(rda.dip.S,permutations = 9999, by = "axis")
anova.dip.S.margin  <- anova(rda.dip.S,permutations = 9999, by = "margin")

anova.mul.S         <- anova(rda.mul.S,permutations = 9999)
anova.mul.S.terms   <- anova(rda.mul.S,permutations = 9999, by = "terms")
anova.mul.S.axis    <- anova(rda.mul.S,permutations = 9999, by = "axis")
anova.mul.S.margin  <- anova(rda.mul.S,permutations = 9999, by = "margin")

# variance explained by axes
summary(rda.dip.S)$cont$importance[,1:5]
summary(rda.mul.S)$cont$importance[,1:5]

## variation partitioning ####
anova.dip.S.margin
Y.dip.env <- Y.dip.pre[,c("SST_winter_mean", "SSS_yearly_mean")]
Y.dip.hab <- Y.dip.pre[,c("habitatPC2")]
Y.dip.MEM <- Y.dip.pre[,c("MEM6","MEM9","MEM7")]
Y.dip.AEM <- Y.dip.pre[,c("AEM23", "AEM52","AEM51", "AEM35")]
Y.dip.spatial <- cbind(Y.dip.MEM, Y.dip.AEM)

vp.dip <- varpart(X.dip, Y.dip.env, Y.dip.hab, Y.dip.spatial)
vp.dip

anova.mul.S.margin
Y.mul.env <- Y.mul.pre[,c("SSS_yearly_mean", "SST_summer_mean", "chloSurfaceMean")]
Y.mul.hab <- Y.mul.pre[,c("habitatPC2", "habitatPC4", "habitatPC3")]
Y.mul.MEM <- Y.mul.pre[,c("MEM5","MEM12","MEM1","MEM17")]
Y.mul.AEM <- Y.mul.pre[,c("AEM40", "AEM3", "AEM41")]
Y.mul.spatial <- cbind(Y.mul.MEM, Y.mul.AEM)

vp.mul <- varpart(X.mul, Y.mul.env, Y.mul.hab, Y.mul.spatial)
vp.mul

### Represent the partitioning using a Venn diagram
     pdf("figures/varpart_dip_venn.pdf", width = 5, height = 5)
plot(vp.dip, digits= 2, bg = 2:4, Xnames = c("ENV", "HAB", "SPA"))
     dev.off()
     pdf("figures/varpart_mul_venn.pdf", width = 5, height = 5)
plot(vp.mul, digits = 2, bg = 2:5, Xnames = c("ENV", "HAB", "SPA"))
     dev.off()

### custom Venn diagram with proportional size
library(VennDiagram)
library(eulerr)

fit.mul.c   <- euler(c("ENV"       = round(vp.mul$part$indfract$Adj.R.square[1],3),
                       "HAB"         = round(vp.mul$part$indfract$Adj.R.square[2],3),
                       "SPA"         = round(vp.mul$part$indfract$Adj.R.square[3],3),
                       "ENV&HAB"     = round(vp.mul$part$indfract$Adj.R.square[4],3),
                       "ENV&SPA"     = round(vp.mul$part$indfract$Adj.R.square[5],3),
                       "HAB&SPA"     = round(vp.mul$part$indfract$Adj.R.square[6],3),
                       "ENV&HAB&SPA" = round(vp.mul$part$indfract$Adj.R.square[7],3)),
                   shape = "circle")
fit.mul.e   <- euler(c("ENV"         = round(vp.mul$part$indfract$Adj.R.square[1],3),
                       "HAB"         = round(vp.mul$part$indfract$Adj.R.square[2],3),
                       "SPA"         = round(vp.mul$part$indfract$Adj.R.square[3],3),
                       "ENV&HAB"     = round(vp.mul$part$indfract$Adj.R.square[4],3),
                       "ENV&SPA"     = round(vp.mul$part$indfract$Adj.R.square[5],3),
                       "HAB&SPA"     = round(vp.mul$part$indfract$Adj.R.square[6],3),
                       "ENV&HAB&SPA" = round(vp.mul$part$indfract$Adj.R.square[7],3)),
                     shape = "ellipse")
     pdf("figures/varpart_mul_euler_ellipse.pdf", width = 5, height = 5)
plot(fit.mul.e, key = TRUE, counts = TRUE, 
     quantities = list(type = c("counts"), font=3, round=2, cex=0.8), 
     fills =list(fill=c(RColorBrewer::brewer.pal(3, "Greys"))),alpha = 0.7,
     factor_names = TRUE, labels=list(font=2, cex=1), legend = FALSE)
     dev.off()
     
     pdf("figures/varpart_mul_euler_ellipse_percent.pdf", width = 5, height = 5)
     plot(fit.mul.e, key = TRUE, counts = TRUE, 
          quantities = list(type = c("percent"), font=3, round=2, cex=0.8), 
          fills =list(fill=c(RColorBrewer::brewer.pal(3, "Greys"))),alpha = 0.7,
          factor_names = TRUE, labels=list(font=2, cex=1), legend = FALSE)
     dev.off()

fit.dip.c   <- euler(c("ENV"         = round(vp.dip$part$indfract$Adj.R.square[1],3),
                       "HAB"         = round(vp.dip$part$indfract$Adj.R.square[2],3),
                       "SPA"         = round(vp.dip$part$indfract$Adj.R.square[3],3), 
                       "ENV&HAB"     = round(vp.dip$part$indfract$Adj.R.square[4],3), 
                       "ENV&SPA"     = round(vp.dip$part$indfract$Adj.R.square[5],3), 
                       "HAB&SPA"     = round(vp.dip$part$indfract$Adj.R.square[6],3), 
                       "ENV&HAB&SPA" = 0),
                   shape = "circle")
fit.dip.e <- euler(c("ENV"         = round(vp.dip$part$indfract$Adj.R.square[1],3),
                     "HAB"         = round(vp.dip$part$indfract$Adj.R.square[2],3),
                     "SPA"         = round(vp.dip$part$indfract$Adj.R.square[3],3), 
                     "ENV&HAB"     = round(vp.dip$part$indfract$Adj.R.square[4],3), 
                     "ENV&SPA"     = round(vp.dip$part$indfract$Adj.R.square[5],3), 
                     "HAB&SPA"     = round(vp.dip$part$indfract$Adj.R.square[6],3), 
                     "ENV&HAB&SPA" = 0), 
                   shape = "ellipse")
     pdf("figures/varpart_dip_euler_ellipse.pdf", width = 5, height = 5)
plot(fit.dip.e, key = TRUE, counts = TRUE, 
     quantities = list(type = c("counts"), font=3, round=2, cex=0.8), 
     fills =list(fill=c(RColorBrewer::brewer.pal(3, "Greys"))),alpha = 0.7,
     factor_names = TRUE, labels=list(font=2, cex=1), legend = FALSE)
     dev.off()
     
     pdf("figures/varpart_dip_euler_ellipse_percent.pdf", width = 5, height = 5)
     plot(fit.dip.e, key = TRUE, counts = TRUE, 
          quantities = list(type = c("percent"), font=3, round=2, cex=0.8), 
          fills =list(fill=c(RColorBrewer::brewer.pal(3, "Greys"))),alpha = 0.7,
          factor_names = TRUE, labels=list(font=2, cex=1), legend = FALSE)
     dev.off()
     


grid.newpage()
draw.triple.venn(area1 = vp.mul$part$fract$Adj.R.square[1],
                 area2 = vp.mul$part$fract$Adj.R.square[2],
                 area3 = vp.mul$part$fract$Adj.R.square[3],
                 n12 = 0.01401542, 
                 n13 = 0.02198+0.01289, 
                 n23 = 0.01318+0.01289,
                 n123 = 0.01289, scaled = T)


#### plotting ####
plot(rda.dip.S)
plot(rda.mul.S)

# add ecoregion colours
eco.dip <- factor(expl.dip$Ecoregion_adj)
eco.mul <- factor(expl.mul$Ecoregion_adj)
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecoregions_adj

plot(rda.dip.S, type="n", scaling=1, main = "Diplodus sargus | 494 adaptive SNPs")
points(rda.dip.S, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda.dip.S, scaling=1, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.dip]) # the individuals
text(rda.dip.S,display="bp", col="black", cex=1, scaling=3)  # the predictors
text(x = -11, y = 2.4, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.dip.S)$r.squared, 5)))
text(x = -11, y = 2.1, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.dip.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.dip), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
   # save as pdf
   #dev.print(pdf, 'figures/rda_dip_adpt_1.pdf')

plot(rda.mul.S, type="n", main = "Mullus surmuletus | 1 558 adaptive SNPs")
points(rda.mul.S, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.mul.S, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.mul]) # the individuals
text(rda.mul.S,display="bp", col="black", cex=1, scaling=3)                              # the predictors
text(x = -17, y = 7, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.mul.S)$r.squared, 5)))
text(x = -17, y = 6.5, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.mul.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.mul), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
   #dev.print(pdf, 'figures/rda_mul_adpt_1.pdf')

# add dapc colours
dapc.dip <- read.table("../02-DAPC/01-Diplodus/c-adaptive/DAPCgrp_dip_adaptive_K3.txt", header = T)
rownames(rda.dip.S$Ybar) == dapc.dip$INDIVIDUALS # check if individuals in same order
dapc.dip <- factor(dapc.dip$STRATA)

dapc.mul <- read.table("../02-DAPC/02-Mullus/c-adaptive/DAPCgrp_mul_adaptive_K3.txt", header = T)
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
# Diplodus non spatial terms
rda.part.dip.hab2 <- rda(X.dip ~ habitatPC2 +
                          Condition(MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+MEM7+SST_winter_mean+SSS_yearly_mean),
                        data = Y.dip.pre)
rda.part.dip.sss <- rda(X.dip ~ SSS_yearly_mean +
                          Condition(MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+habitatPC2+MEM7+SST_winter_mean),
                        data = Y.dip.pre)
rda.part.dip.sst <- rda(X.dip ~ SST_winter_mean +
                          Condition(MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+habitatPC2+MEM7+SSS_yearly_mean),
                        data = Y.dip.pre)

apd.hab2 <- anova(rda.part.dip.hab2, permutations = 9999)  # 0.0191 *
apd.sst  <- anova(rda.part.dip.sst , permutations = 9999)  # 8e-04 ***
apd.sss  <- anova(rda.part.dip.sss , permutations = 9999)  # 0.019 *

RsquareAdj(rda.part.dip.sst)  # 0.003189305
RsquareAdj(rda.part.dip.sss)  # 0.001906393
RsquareAdj(rda.part.dip.hab2) # 0.001555188

# Diplodus variation partitions
rda.part.dip.hab <- rda(X.dip ~ habitatPC2 +
                           Condition(MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+MEM7+SST_winter_mean+SSS_yearly_mean),
                         data = Y.dip.pre)
rda.part.dip.env <- rda(X.dip ~ SSS_yearly_mean +SST_winter_mean+
                          Condition(MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+habitatPC2+MEM7),
                        data = Y.dip.pre)
rda.part.dip.spa <- rda(X.dip ~ MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+MEM7+
                          Condition(habitatPC2+SSS_yearly_mean +SST_winter_mean),
                        data = Y.dip.pre)

apd.hab <- anova(rda.part.dip.hab, permutations = 9999) # 0.0183 *
apd.env <- anova(rda.part.dip.env, permutations = 9999) # 0.0078 **
apd.spa <- anova(rda.part.dip.spa, permutations = 9999) # 1e-04 ***

RsquareAdj(rda.part.dip.hab) # 0.001555188
RsquareAdj(rda.part.dip.env) # 0.003242096
RsquareAdj(rda.part.dip.spa) # 0.09121215

# Mullus non spatial terms
rda.part.mul.chla <- rda(X.mul ~ chloSurfaceMean+
                      Condition(SSS_yearly_mean +MEM5+AEM40+MEM12+MEM1+AEM3+MEM17+AEM41+habitatPC2+habitatPC4+habitatPC3+SST_summer_mean),
                    data = Y.mul.pre)
rda.part.mul.sst <- rda(X.mul ~ SST_summer_mean +
                      Condition(SSS_yearly_mean +MEM5+AEM40+MEM12+MEM1+chloSurfaceMean+AEM3+MEM17+AEM41+habitatPC2+habitatPC4+habitatPC3),
                    data = Y.mul.pre)
rda.part.mul.sss <- rda(X.mul ~ SSS_yearly_mean +
                      Condition(MEM5+AEM40+MEM12+MEM1+chloSurfaceMean+AEM3+MEM17+AEM41+habitatPC2+habitatPC4+habitatPC3+SST_summer_mean),
                    data = Y.mul.pre)
rda.part.mul.hab2 <- rda(X.mul ~ habitatPC2+
                          Condition(SSS_yearly_mean+MEM5+AEM40+MEM12+MEM1+chloSurfaceMean+AEM3+MEM17+AEM41+habitatPC4+habitatPC3+SST_summer_mean),
                        data = Y.mul.pre)
rda.part.mul.hab4 <- rda(X.mul ~ habitatPC4+
                          Condition(SSS_yearly_mean+MEM5+AEM40+MEM12+MEM1+chloSurfaceMean+AEM3+MEM17+AEM41+habitatPC2+habitatPC3+SST_summer_mean),
                        data = Y.mul.pre)
rda.part.mul.hab3 <- rda(X.mul ~ habitatPC3+
                          Condition(SSS_yearly_mean+MEM5+AEM40+MEM12+MEM1+chloSurfaceMean+AEM3+MEM17+AEM41+habitatPC2+habitatPC4+SST_summer_mean),
                        data = Y.mul.pre)

apm.sss  <- anova(rda.part.mul.sss , permutations = 9999)   # 0.0012 **
apm.chla <- anova(rda.part.mul.chla, permutations = 9999)   # 1e-04 ***
apm.sst  <- anova(rda.part.mul.sst , permutations = 9999)   # 0.0541 .
apm.hab2 <- anova(rda.part.mul.hab2, permutations = 9999)   # 0.022 *
apm.hab4 <- anova(rda.part.mul.hab4, permutations = 9999)   # 0.0031 **
apm.hab3 <- anova(rda.part.mul.hab3, permutations = 9999)   # 0.0529 .


RsquareAdj(rda.part.mul.chla)  # 0.003407195
RsquareAdj(rda.part.mul.sst)   # 0.0006874253
RsquareAdj(rda.part.mul.sss)   # 0.003035332
RsquareAdj(rda.part.mul.hab2)  # 0.001048933
RsquareAdj(rda.part.mul.hab4)  # 0.00193056
RsquareAdj(rda.part.mul.hab3)  # 0.0007496847

# Mullus variation partitions
rda.part.mul.hab <- rda(X.mul ~ habitatPC2+habitatPC4+habitatPC3+
                           Condition(SSS_yearly_mean+MEM5+AEM40+MEM12+MEM1+chloSurfaceMean+AEM3+MEM17+AEM41+SST_summer_mean),
                         data = Y.mul.pre)
rda.part.mul.env <- rda(X.mul ~ SSS_yearly_mean+chloSurfaceMean+SST_summer_mean +
                          Condition(MEM5+AEM40+MEM12+MEM1+AEM3+MEM17+AEM41+habitatPC2+habitatPC4+habitatPC3),
                        data = Y.mul.pre)
rda.part.mul.spa <- rda(X.mul ~ MEM5+AEM40+MEM12+MEM1+AEM3+MEM17+AEM41+
                          Condition(SSS_yearly_mean+chloSurfaceMean+habitatPC2+habitatPC4+habitatPC3+SST_summer_mean),
                        data = Y.mul.pre)
rda.mul.pc2 <- rda(X.mul ~ habitatPC2, data = Y.mul.pre)
anova(rda.mul.pc2)
RsquareAdj(rda.mul.pc2)

apm.hab <- anova(rda.part.mul.hab, permutations = 9999) # 9e-04 ***
apm.env <- anova(rda.part.mul.env, permutations = 9999) # 1e-04 ***
apm.spa <- anova(rda.part.mul.spa, permutations = 9999) # 1e-04 ***

RsquareAdj(rda.part.mul.hab) # 0.003663831
RsquareAdj(rda.part.mul.env) # 0.01144185
RsquareAdj(rda.part.mul.spa) # 0.03977473

#### iter 1c: partial RDA with only habitat ####
rda.part.dip.hab <- rda(X.dip ~ habitatPC2 +
                          Condition(MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+MEM7+ SSS_yearly_mean + SST_yearly_mean),
                        data = Y.dip.pre)
anova(rda.part.dip.hab, permutations = 999)
RsquareAdj(rda.part.dip.hab)
plot(rda.part.dip.hab)

rda.part.mul.hab <- rda(X.mul ~ habitatPC2 + habitatPC4 +
                      Condition(MEM5+AEM40+MEM12+MEM1+AEM3+MEM17+AEM41+SST_yearly_mean+SSS_yearly_mean),
                    data = Y.mul.pre)
anova(rda.part.mul.hab, permutations = 999)
RsquareAdj(rda.part.mul.hab)
plot(rda.part.mul.hab)

anova(rda.part.mul.hab, by = "terms", permutations = 999)
###
rda.part.dip.SS <- rda(X.dip ~ SSS_yearly_mean + SST_yearly_mean +
                          Condition(MEM6+AEM23+MEM9+AEM52+AEM51+AEM35+MEM7+habitatPC2 ),
                        data = Y.dip.pre)
anova(rda.part.dip.SS, permutations = 999)
RsquareAdj(rda.part.dip.SS)
plot(rda.part.dip.SS)

rda.part.mul.sss <- rda(X.mul ~ SSS_yearly_mean +
                          Condition(MEM5+AEM40+MEM12+MEM1+AEM3+MEM17+AEM41+SST_yearly_mean+habitatPC2 + habitatPC4),
                        data = Y.mul.pre)
anova(rda.part.mul.sss, permutations = 999)
RsquareAdj(rda.part.mul.sss)
plot(rda.part.mul.sss)
rda.part.mul.sst <- rda(X.mul ~ SST_yearly_mean +
                          Condition(MEM5+AEM40+MEM12+MEM1+AEM3+MEM17+AEM41+SSS_yearly_mean+habitatPC2 + habitatPC4),
                        data = Y.mul.pre)
anova(rda.part.mul.sst, permutations = 999)
RsquareAdj(rda.part.mul.sst)
plot(rda.part.mul.sst)

rda.part.mul.sss.sst <- rda(X.mul ~ SST_yearly_mean+SSS_yearly_mean +
                          Condition(MEM5+AEM40+MEM12+MEM1+AEM3+MEM17+AEM41+habitatPC2 + habitatPC4),
                        data = Y.mul.pre)
anova(rda.part.mul.sss.sst, permutations = 999)
RsquareAdj(rda.part.mul.sss.sst)
plot(rda.part.mul.sss.sst)

anova(rda.part.mul.hab, by = "terms", permutations = 999)

#### iter 2: only env variables ####
# take the first four dbMEMs and AEMs. order of variables stems for preparatory ordistep procedure

# for dip
colnames(Y.dip)
Y.dip.pre2 <- Y.dip %>% select("SSS_yearly_mean",
                              "SST_winter_mean", "SST_summer_mean",
                              "chloSurfaceMean")
                              "habitatPC1","habitatPC2","habitatPC3","habitatPC4")

# for mul: MEM3 and AEM1 (corr with SST) and AEM5 (least important from the two, based on the order of ordistep results)
colnames(Y.mul)
Y.mul.pre2 <- Y.mul %>% select("SSS_yearly_mean",
                              "SST_winter_mean", "SST_summer_mean",
                              "chloSurfaceMean")
                              #"habitatPC1","habitatPC2","habitatPC3","habitatPC4")

# run rda
rda.dip.pre2 <- rda(X.dip, Y.dip.pre2)
rda.mul.pre2 <- rda(X.mul, Y.mul.pre2)

vif(rda.dip.pre) 
vif(rda.mul.pre) 

# look at explained variance
RsquareAdj(rda.dip.pre2) # adj.r.squared = 0.1578281
RsquareAdj(rda.mul.pre2) # adj.r.squared = 0.1049795

# check model significance
anova(rda.dip.pre2, perm = 9999) # P = 0.001
anova(rda.mul.pre2, perm = 9999) # P = 0.001

## Ordi2Step to select most informative variables
# start working form an empty model
rda2.dip.0 <- rda(X.dip ~ 1, Y.dip.pre2)
rda2.mul.0 <- rda(X.mul ~ 1, Y.mul.pre2)

# OrdiR2step will move towards the global model with all explanatory variables
rda2.dip.G<- rda(X.dip ~ ., Y.dip.pre2)
rda2.mul.G<- rda(X.mul ~ ., Y.mul.pre2)

# run selection of variables
Sel.dip2 <- ordiR2step(rda2.dip.0, scope = formula(rda2.dip.G), direction="both") 
Sel.mul2 <- ordiR2step(rda2.mul.0, scope = formula(rda2.mul.G), direction="both") 

# check out selected variables 
Sel.dip2$anova
Sel.mul2$anova

#### build final model 
# Build a model with the selected variables and visualize the results
Sel.dip.vars2 <-  rownames(Sel.dip2$anova)[-nrow(Sel.dip2$anova)] %>% sub("..", "",.)
# [1] "SST_summer_mean" "SST_winter_mean" "SSS_yearly_mean"
# [4] "habitatPC1"      "chloSurfaceMean"
Sel.mul.vars2 <-  rownames(Sel.mul2$anova)[-nrow(Sel.mul2$anova)] %>% sub("..", "",.)
# [1] "SSS_yearly_mean" "SST_summer_mean" "habitatPC4"     
# [4] "habitatPC2"      "chloSurfaceMean" "habitatPC3" 

Y.dip.sel2 <- Y.dip %>% select(all_of(Sel.dip.vars2)) # fill in vars selected by ordistep
Y.mul.sel2 <- Y.mul %>% select(all_of(Sel.mul.vars2)) 

# run final rda
rda.dip.S2<- rda(X.dip~.,Y.dip.sel2)
rda.mul.S2<- rda(X.mul~.,Y.mul.sel2)

# export rda
#saveRDS(rda.dip.S, file = "output/rda.4sel.dip.adpt.rds")
#saveRDS(rda.mul.S, file = "output/rda.4sel.mul.adpt.rds")

summary(rda.dip.S2)   
summary(rda.mul.S2)   

RsquareAdj(rda.dip.S2) # 0.1551034
RsquareAdj(rda.mul.S2)

vif(rda.dip.S2)
vif(rda.mul.S2)

anova.dip.S2         <- anova(rda.dip.S2,permutations = 9999)
anova.dip.S2.terms   <- anova(rda.dip.S2,permutations = 9999, by = "terms")
anova.dip.S2.axis    <- anova(rda.dip.S2,permutations = 9999, by = "axis")
anova.dip.S2.margin  <- anova(rda.dip.S2,permutations = 9999, by = "margin")

anova.mul.S2         <- anova(rda.mul.S2,permutations = 9999)
anova.mul.S2.terms   <- anova(rda.mul.S2,permutations = 9999, by = "terms")
anova.mul.S2.axis    <- anova(rda.mul.S2,permutations = 9999, by = "axis")
anova.mul.S2.margin  <- anova(rda.mul.S2,permutations = 9999, by = "margin")

# variance explained by axes
summary(rda.dip.S)$cont$importance[,1:5]
summary(rda.mul.S)$cont$importance[,1:5]

#### iter 3: env + hab variables ####
# take the first four dbMEMs and AEMs. order of variables stems for preparatory ordistep procedure

# for dip
colnames(Y.dip)
Y.dip.pre3 <- Y.dip %>% select("SSS_yearly_mean",
                               "SST_winter_mean", "SST_summer_mean",
                               "chloSurfaceMean",
                               "habitatPC1","habitatPC2","habitatPC3","habitatPC4")

# for mul: MEM3 and AEM1 (corr with SST) and AEM5 (least important from the two, based on the order of ordistep results)
colnames(Y.mul)
Y.mul.pre3 <- Y.mul %>% select("SSS_yearly_mean",
                               "SST_winter_mean", "SST_summer_mean",
                               "chloSurfaceMean",
                               "habitatPC1","habitatPC2","habitatPC3","habitatPC4")

## Ordi3Step to select most informative variables
# start working form an empty model
rda3.dip.0 <- rda(X.dip ~ 1, Y.dip.pre3)
rda3.mul.0 <- rda(X.mul ~ 1, Y.mul.pre3)

# OrdiR3step will move towards the global model with all explanatory variables
rda3.dip.G<- rda(X.dip ~ ., Y.dip.pre3)
rda3.mul.G<- rda(X.mul ~ ., Y.mul.pre3)

# run selection of variables
Sel.dip3 <- ordiR2step(rda3.dip.0, scope = formula(rda3.dip.G), direction="both") 
Sel.mul3 <- ordiR2step(rda3.mul.0, scope = formula(rda3.mul.G), direction="both") 

# check out selected variables 
Sel.dip3$anova
Sel.mul3$anova

#### build final model 
# Build a model with the selected variables and visualize the results
Sel.dip.vars3 <-  rownames(Sel.dip3$anova)[-nrow(Sel.dip3$anova)] %>% sub("..", "",.)
# [1] "SST_summer_mean" "SST_winter_mean" "SSS_yearly_mean"
# [4] "habitatPC1"      "chloSurfaceMean"
Sel.mul.vars3 <-  rownames(Sel.mul3$anova)[-nrow(Sel.mul3$anova)] %>% sub("..", "",.)
# [1] "SSS_yearly_mean" "SST_summer_mean" "habitatPC4"     
# [4] "habitatPC2"      "chloSurfaceMean" "habitatPC3" 

Y.dip.sel3 <- Y.dip %>% select(all_of(Sel.dip.vars3)) # fill in vars selected by ordistep
Y.mul.sel3 <- Y.mul %>% select(all_of(Sel.mul.vars3)) 

# run final rda
rda.dip.S3<- rda(X.dip~.,Y.dip.sel3)
rda.mul.S3<- rda(X.mul~.,Y.mul.sel3)

# export rda
#saveRDS(rda.dip.S, file = "output/rda.4sel.dip.adpt.rds")
#saveRDS(rda.mul.S, file = "output/rda.4sel.mul.adpt.rds")

summary(rda.dip.S3)   
summary(rda.mul.S3)   

RsquareAdj(rda.dip.S3) # 0.1551034
RsquareAdj(rda.mul.S3)

vif(rda.dip.S3)
vif(rda.mul.S3)

anova.dip.S3         <- anova(rda.dip.S3,permutations = 9999)
anova.dip.S3.terms   <- anova(rda.dip.S3,permutations = 9999, by = "terms")
anova.dip.S3.axis    <- anova(rda.dip.S3,permutations = 9999, by = "axis")
anova.dip.S3.margin  <- anova(rda.dip.S3,permutations = 9999, by = "margin")

anova.mul.S3         <- anova(rda.mul.S3,permutations = 9999)
anova.mul.S3.terms   <- anova(rda.mul.S3,permutations = 9999, by = "terms")
anova.mul.S3.axis    <- anova(rda.mul.S3,permutations = 9999, by = "axis")
anova.mul.S3.margin  <- anova(rda.mul.S3,permutations = 9999, by = "margin")

# variance explained by axes
summary(rda.dip.S3$cont$importance[,1:5])
summary(rda.mul.S3$cont$importance[,1:5])


#### iter 4: 4sel + env, no hab variables ####
# take the first four dbMEMs and AEMs. order of variables stems for preparatory ordistep procedure
# for dip
Y.dip.pre4 <- Y.dip %>% select("SSS_yearly_mean",
                               "SST_winter_mean", "SST_summer_mean",
                               "chloSurfaceMean",
                               "MEM6","MEM16","MEM7","MEM9",
                               "AEM23","AEM35","AEM51","AEM52")
Y.mul.pre4 <- Y.mul %>% select("SSS_yearly_mean",
                               "SST_winter_mean", "SST_summer_mean",
                               "chloSurfaceMean",
                               "MEM5","MEM1","MEM12","MEM17",
                               "AEM41","AEM40","AEM17","AEM3")
## Ordi2Step to select most informative variables
# start working form an empty model
rda4.dip.0 <- rda(X.dip ~ 1, Y.dip.pre4)
rda4.mul.0 <- rda(X.mul ~ 1, Y.mul.pre4)

# OrdiR4step will move towards the global model with all explanatory variables
rda4.dip.G<- rda(X.dip ~ ., Y.dip.pre4)
rda4.mul.G<- rda(X.mul ~ ., Y.mul.pre4)

# run selection of variables
Sel.dip4 <- ordiR2step(rda4.dip.0, scope = formula(rda4.dip.G), direction="both") 
Sel.mul4 <- ordiR2step(rda4.mul.0, scope = formula(rda4.mul.G), direction="both") 

# check out selected variables 
Sel.dip4$anova
Sel.mul4$anova

#### build final model 
# Build a model with the selected variables and visualize the results
Sel.dip.vars4 <-  rownames(Sel.dip4$anova)[-nrow(Sel.dip4$anova)] %>% sub("..", "",.)
# [1] "MEM6"            "AEM23"           "MEM9"           
# [4] "AEM52"           "AEM51"           "AEM35"          
# [7] "MEM7"            "SSS_yearly_mean" "SST_winter_mean"
Sel.mul.vars4 <-  rownames(Sel.mul4$anova)[-nrow(Sel.mul4$anova)] %>% sub("..", "",.)
#  [1] "SSS_yearly_mean" "MEM5"            "AEM40"          
# [4] "MEM12"           "MEM1"            "chloSurfaceMean"
# [7] "AEM3"            "MEM17"           "AEM41"          
# [10] "SST_summer_mean"

Y.dip.sel4 <- Y.dip %>% select(all_of(Sel.dip.vars4)) # fill in vars selected by ordistep
Y.mul.sel4 <- Y.mul %>% select(all_of(Sel.mul.vars4)) 

# run final rda
rda.dip.S4<- rda(X.dip~.,Y.dip.sel4)
rda.mul.S4<- rda(X.mul~.,Y.mul.sel4)

# export rda
#saveRDS(rda.dip.S, file = "output/rda.4sel.dip.adpt.rds")
#saveRDS(rda.mul.S, file = "output/rda.4sel.mul.adpt.rds")

RsquareAdj(rda.dip.S4) 
RsquareAdj(rda.mul.S4)

vif(rda.dip.S4)
vif(rda.mul.S4)

anova.dip.S4         <- anova(rda.dip.S4,permutations = 9999)
anova.dip.S4.terms   <- anova(rda.dip.S4,permutations = 9999, by = "terms")
anova.dip.S4.axis    <- anova(rda.dip.S4,permutations = 9999, by = "axis")
anova.dip.S4.margin  <- anova(rda.dip.S4,permutations = 9999, by = "margin")

anova.mul.S4         <- anova(rda.mul.S4,permutations = 9999)
anova.mul.S4.terms   <- anova(rda.mul.S4,permutations = 9999, by = "terms")
anova.mul.S4.axis    <- anova(rda.mul.S4,permutations = 9999, by = "axis")
anova.mul.S4.margin  <- anova(rda.mul.S4,permutations = 9999, by = "margin")

# variance explained by axes
summary(rda.dip.S4)$cont$importance[,1:5]
summary(rda.mul.S4)$cont$importance[,1:5]
                
## variation partitioning ####
anova.dip.S4.margin
Y.dip.sss <- Y.dip.pre[,c("SSS_yearly_mean")]
Y.dip.sst <- Y.dip.pre[,c("SST_winter_mean")]
Y.dip.hab <- Y.dip.pre[,c("habitatPC2", "habitatPC1", "habitatPC4")]
Y.dip.MEM <- Y.dip.pre[,c("MEM6","MEM9","MEM7")]
Y.dip.AEM <- Y.dip.pre[,c("AEM23", "AEM52","AEM51", "AEM35")]
Y.dip.spatial <- cbind(Y.dip.MEM, Y.dip.AEM)

vp.dip4 <- varpart(X.dip, Y.dip.sss, Y.dip.sst, Y.dip.spatial, Y.dip.hab)
vp.dip4

anova.mul.S.margin
Y.mul.sss <- Y.mul.pre[,c("SSS_yearly_mean")]
Y.mul.sst <- Y.mul.pre[,c("SST_summer_mean")]
Y.mul.chl <- Y.mul.pre[,c("chloSurfaceMean")]
Y.mul.env <- cbind(Y.mul.sss, Y.mul.sst, Y.mul.chl)
Y.mul.hab <- Y.mul.pre[,c("habitatPC2", "habitatPC4", "habitatPC3")]
Y.mul.MEM <- Y.mul.pre[,c("MEM5","MEM12","MEM1","MEM17")]
Y.mul.AEM <- Y.mul.pre[,c("AEM40", "AEM3", "AEM41")]
Y.mul.spatial <- cbind(Y.mul.MEM, Y.mul.AEM)

vp.mul4 <- varpart(X.mul, Y.mul.env, Y.mul.spatial, Y.mul.hab)
vp.mul4

### Represent the partitioning using a Venn diagram
plot(vp.dip4, digits= 2, bg = 2:4, Xnames = c("SSS", "SST", "SPA", "HAB"))


plot(vp.mul4, digits = 2, bg = 2:5, Xnames = c("ENV", "SPA","HAB"), main= "M.surmuletus - outlier SNPs")

                
#### iter 6: VIF on MEM & AEM ####

# run VIF fonction on MEM and AEM variables, then add SSS and habitat again
# select spatial variables
Y.dip.spatial <- select(Y.dip, starts_with(c("MEM","AEM")))
Y.mul.spatial <- select(Y.mul, starts_with(c("MEM","AEM")))
# calculate correlation matrix
cor.dip.spatial <- cor(Y.dip.spatial)
cor.mul.spatial <- cor(Y.mul.spatial)
# apply VIF function to stepwise remove vars
source("vif_func.R")
vif.dip.spatial <- vif_func(cor.dip.spatial,thresh=10,trace=T)
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

# check out vif
vif(Sel.spVif.dip)
vif(Sel.spVif.mul)

# check out spatial scale selected MEM
s.value(expl.dip[,4:5], expl.dip$MEM16)
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
anova(rda.spVif.dip.S, permutations = 1000, by = "terms")

anova(rda.spVif.mul.S, permutations = 1000)
anova(rda.spVif.mul.S, permutations = 1000, by = "terms")
anova(rda.spVif.mul.S, permutations = 1000, by = "axis")

#anova(rda.spVif.dip.S, by = "axis")
#anova(rda.spVif.mul.S, by = "axis")

plot(rda.spVif.dip.S)
plot(rda.spVif.mul.S)

### plotting
# add ecoregion colours
eco.dip <- factor(expl.dip$Ecoregion_adj)
eco.mul <- factor(expl.mul$Ecoregion_adj)
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecoregions_adj
# dapc clusters
dapc.dip <- read.csv("../02-DAPC/01-Diplodus/c-adaptive/DAPCassign_dip_adaptive_K3.csv", row.names = 1)
rownames(dapc.dip) <- dapc.dip$IND
dapc.dip <- dapc.dip[rownames(rda.spVif.dip.S$Ybar),]
dapc.dip <- factor(dapc.dip$cluster)

dapc.mul <- read.table("../02-DAPC/02-Mullus/c-adaptive/DAPCassign_mul_adaptive_K7.txt", header = T)
rownames(dapc.mul) <- dapc.mul$INDIVIDUALS
dapc.mul <- dapc.mul[rownames(rda.spVif.mul.S$Ybar),]
dapc.mul <- factor(dapc.mul$STRATA)

standard_12=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")

plot(rda.spVif.dip.S)
plot(rda.spVif.dip.S, type="n", main = "Diplodus sargus | xx adaptive SNPs")
points(rda.spVif.dip.S, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.spVif.dip.S, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.dip]) # the individuals
points(rda.spVif.dip.S, display="sites", pch=21, cex=1.3, col="gray32", bg=standard_12[dapc.dip]) # the individuals
text(rda.spVif.dip.S,display="bp", col="black", cex=1)  # the predictors
text(x = -10, y = 6, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.spVif.dip.S)$r.squared, 5)))
text(x = -10, y = 5.5, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.spVif.dip.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.dip), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_dip_adpt_spatialVIF_ordi2step.pdf')

plot(rda.spVif.mul.S)
plot(rda.spVif.mul.S, type="n", main = "Mullus surmuletus | xx adaptive SNPs")
points(rda.spVif.mul.S, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.spVif.mul.S, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.mul]) # the individuals
points(rda.spVif.mul.S, display="sites", pch=21, cex=1.3, col="gray32", bg=standard_12[dapc.mul]) # the individuals
text(rda.spVif.mul.S,display="bp", col="black", cex=1)                              # the predictors
text(x = -10, y = 4, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.spVif.mul.S)$r.squared, 5)))
text(x = -10, y = 3.5, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.spVif.mul.S)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.mul), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_mul_adpt_spatialVIF_ordi2step.pdf')

# correlation between vars
pairs.panels(Y.mul.spVif.sel)

# get individual scores to find individuals clustered separately
scores(rda.spVif.mul.S)$sites %>%
  as.data.frame() %>% 
  mutate(ind = rownames(.)) %>% 
  filter(CAP1 < -2 & CAP2 < -2) 


#### iter 6b: partial RDA with significant vars iter2 ####

# extract the spatial variables retained after vif & ordi2step
Y.dip.spVif.sel %>% select(starts_with(c("MEM", "AEM"))) %>% colnames() %>% paste(., collapse="+")
Y.mul.spVif.sel %>% select(starts_with(c("MEM", "AEM"))) %>% colnames() %>% paste(., collapse="+")
# use these as conditional variables for partial rda
rda.part.dip <- capscale(X.dip ~ SSS_yearly_mean + habitatPC4 + 
                           Condition(MEM16+MEM9+MEM17+AEM3+AEM4+AEM16+AEM10+AEM18+AEM7),
                         data = Y.dip)
rda.part.mul <- capscale(X.mul ~ SSS_yearly_mean + habitatPC1 + habitatPC2 + 
                           Condition(MEM5+MEM18+MEM7+MEM14+MEM17+MEM16+MEM2+MEM13+MEM15+MEM3+AEM8+AEM9+AEM2+AEM3+AEM13+AEM4+AEM14+AEM10),
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
plot(rda.part.dip)
plot(rda.part.dip, type="n", main = "Diplodus sargus | 494 adaptive SNPs")
points(rda.part.dip, display="species", pch=20, cex=0.7, col="gray32")           # the SNPs
points(rda.part.dip, display="sites", pch=21, cex=1.3, col="gray32", bg=bg[eco.dip]) # the individuals
text(rda.part.dip,display="bp", col="black", cex=1)  # the predictors
text(x = -18, y = 2, adj = c(0,NA), 
     labels = paste0("Rsquare = ",round(RsquareAdj(rda.part.dip)$r.squared, 5)))
text(x = -18, y = 1.4, adj = c(0,NA), 
     labels = paste0("RsquareAdj = ",round(RsquareAdj(rda.part.dip)$adj.r.squared, 5)))
legend("bottomleft", legend=levels(eco.dip), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# save as pdf
#dev.print(pdf, 'figures/rda_dip_adpt_spatialVIF_ordi2step_partial.pdf')

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



#### iter 7: no VIF ####

# run RDA and check out VIF (global model to use for Ordi2Step)
rda1.dip <- rda(X.dip ~ ., Y.dip)
rda1.mul <- rda(X.mul ~ ., Y.mul)

vif(rda1.dip) # many quite high
vif(rda1.mul) # 

RsquareAdj(rda1.dip) # $adj.r.squared = 0.2083966
RsquareAdj(rda1.mul) # $adj.r.squared = 0.1342474

# run Ordi2Step to select most informative variables
# start working form an empty model
rda0.dip <- rda(X.dip ~ 1, Y.dip)
rda0.mul <- rda(X.mul ~ 1, Y.mul)

# run selection of variables
Sel.dip <- ordiR2step(rda0.dip, scope = formula(rda1.dip), direction="both") 
Sel.mul <- ordiR2step(rda0.mul, scope = formula(rda1.mul), direction="both") 

# check out selected variables 
Sel.dip$anova
Sel.mul$anova

# check out vif
vif(Sel.dip)
vif(Sel.mul)

#### iter 8: (suitable) habitat vars ####
# get detailed habitat vars pre-PCA
substrate    <- read.csv("../03-RDA/00-envdata-EB/seaconnect_substrate_res200m.csv",   stringsAsFactor = FALSE, sep = ";")
head(substrate)
rocky <- substrate$Coarse.and.mixed.sediment + substrate$Rock.or.other.hard.substrata
sandy <- substrate$Coarse.and.mixed.sediment + substrate$Sand + substrate$Sandy.mud + substrate$Sandy.mud.to.muddy.sand

habitat.dip <- as.data.frame(cbind(substrate$SamplingCell, rocky, sandy), stringsAsFactors = F)
colnames(habitat.dip) <-c("SamplingCell", "rocky","sandy")
rownames(habitat.dip) <- habitat.dip$SamplingCell
habitat.dip <- habitat.dip[expl.dip$SamplingCell,]
habitat.dip <- cbind(habitat.dip, expl.dip)
habitat.dip$rocky <- as.numeric(habitat.dip$rocky)
habitat.dip$sandy <- as.numeric(habitat.dip$sandy)
Y.dip.iter3 <- habitat.dip[,c(2:3,14:38)]
head(Y.dip.iter3)

habitat.mul <- as.data.frame(cbind(substrate$SamplingCell, rocky, sandy), stringsAsFactors = F)
colnames(habitat.mul) <-c("SamplingCell", "rocky","sandy")
rownames(habitat.mul) <- habitat.mul$SamplingCell
habitat.mul <- habitat.mul[expl.mul$SamplingCell,]
habitat.mul <- cbind(habitat.mul, expl.mul)
habitat.mul$rocky <- as.numeric(habitat.mul$rocky)
habitat.mul$sandy <- as.numeric(habitat.mul$sandy)
Y.mul.iter3 <- habitat.mul[,c(2:3, 9,14:45)]
head(Y.mul.iter3)

# rda
rda.dip.iter3 <- capscale(X.dip ~ ., data = Y.dip.iter3)
summary(rda.dip.iter3)
vif(rda.dip.iter3)
RsquareAdj(rda.dip.iter3)
anova(rda.dip.iter3, permutations = 999)
anova(rda.dip.iter3, by = "margin")

rda.mul.iter3 <- capscale(X.mul ~ ., data = Y.mul.iter3)
summary(rda.mul.iter3)
vif(rda.mul.iter3)
RsquareAdj(rda.mul.iter3)
anova(rda.mul.iter3, permutations = 999)
anova(rda.mul.iter3, by = "margin")

## Ordi2Step again to select most informative variables
# start working form an empty model
rda.dip.iter3.0 <- rda(X.dip ~ 1, Y.dip.iter3)
rda.mul.iter3.0 <- rda(X.mul ~ 1, Y.mul.iter3)

# OrdiR2step will move towards the global model with all explanatory variables
rda.dip.iter3.G<- rda(X.dip ~ ., Y.dip.iter3)
rda.mul.iter3.G<- rda(X.mul ~ ., Y.mul.iter3)

# run selection of variables
Sel.dip.iter3 <- ordiR2step(rda.dip.iter3.0, scope = formula(rda.dip.iter3.G), direction="both") 
Sel.mul.iter3 <- ordiR2step(rda.mul.iter3.0, scope = formula(rda.mul.iter3.G), direction="both") 

Sel.dip.iter3$anova
Sel.mul.iter3$anova

RsquareAdj(Sel.dip.iter3)
RsquareAdj(Sel.mul.iter3)

vif(Sel.dip.iter3)
vif(Sel.mul.iter3)

## plot