#*****************************************************************************
#******             Detecting population clusters using DAPC            ******
#*****************************************************************************

### tutorial by adegenet
### http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

setwd("/Users/eboulanger-admin/Documents/project_SEACONNECT/seaConnect--radFishComp/02-DAPC/01-Diplodus")

# load packages
library(adegenet)
library(dplyr)
library(parallel)
library(png)
library(stringr)
# for bar chart
library(reshape)
library(ggplot2)
# for pie map
library(maps)
library(mapdata)
library(plotrix)

# set variables
sp <- "dip"
subset <- "neutral"

# import stamps 
data_pic <- list.files(path = "../../00-Misc/stamps/", pattern="*.png",recursive = FALSE)
source_pic <- paste0("../../00-Misc/stamps/", data_pic)
pic <- lapply(source_pic,readPNG)
names(pic) <- str_sub(data_pic,1, -5)

# load data
sampling <-  read.table("../../00-Data/coord_seaconnect_tous.txt", sep = "\t", header = T, stringsAsFactors = F)
if (sp == "mul") {
  all      <-"../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_all_filtered.raw" 
  adaptive <-"../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_adaptive_1558.raw" 
  neutral  <-"../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_neutral_13548.raw"
} else {
  all      <-"../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/dip_all_filtered.raw" 
  adaptive <-"../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/dip_adaptive_494.raw" 
  neutral  <-"../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/dip_neutral_7570.raw"
}
# read plink file
dat <- read.PLINK(get(subset))

### Discriminant Analysis of Principal Components -----

## DAPC with sample sites as a priori cluster
# depracated
# pop(dip_adpt) <- gsub("_.*","", pop(dip_adpt))
# dapc.dip.pop <- dapc(dip_adpt)
# 50
# 2
# scatter(dapc.dip.pop)

## DAPC with ecoregions as a priori clusters ----
sites <- gsub("i.*","", pop(dat)) %>% 
  as.data.frame() %>% 
  `colnames<-` ("SamplingCell")
ecor <- sampling %>% 
  .[c("SamplingCell", "Ecoregion_adj")] %>% 
  right_join(sites)
pop(dat) <- ecor$Ecoregion_adj
dapc.ecor <- dapc(dat, n.pca = round(nrow(dat)/3))
# check the alpha-score
a.score(dapc.ecor) # adaptive: -0.02606066
                   #           which is very low, and indicates that the result is very unstable
                   # neutral: 0.1054604
# find the optimal number of PC's to retain based on the a-score
dapc.ecor.optimscore <- optim.a.score(dapc.ecor) # adaptive: $25
                                                 # neutral : $43
# run the dapc with the optimal number of PCs
dapc.ecor.best <- dapc(dat, n.pca = dapc.ecor.optimscore$best, n.da = 2)
dapc.ecor.best.score <- a.score(dapc.ecor.best, n.sim = 1000) # adaptive: $mean = 0.2095513
                                                              # neutral : $mean =       

 #jpeg(filename = "b-neutral/DAPCplot_dip_neutral_ecor_optim.jpeg", res = 200, width = 1200, height = 1200)
scatter(dapc.ecor.best)
title("DAPC D. sargus adaptive | ecoregions a priori\n")
 dev.off()
 dev.set(dev.prev())

## DAPC without a priori populations ----

# search for clusters
grp <- find.clusters(dat, max.n.clust=20, n.pca = round(nrow(dat)/3)) # n.pca = 1/3 n.ind
# run DAPC with new a priori clusters
dapc.clust <- dapc(dat, grp$grp, n.pca = nrow(dat)/3, n.da = 2)
scatter(dapc.clust)
# find the optimal number of PC's to retain based on the a-score
dapc.clust.optimscore <- optim.a.score(dapc.clust) # neutral : $1 
                                                   # adaptive: $1 for K = 2, $12 for K=3, 15 for K = 4, 13 for K = 5
# run the dapc with the optimal number of PCs
dapc.clust.best <- dapc(dat, grp$grp, n.pca = dapc.clust.optimscore$best, n.da=2)
scatter(dapc.clust.best)
dapc.clust.best.score <- a.score(dapc.clust.best, n.sim = 1000) # neutral : $mean = 0.4999182 for K = 2, n.pc = 1
                                                                # adaptive: $mean = 0.4996529 for K = 2, n.pc = 1
                                                                #           $mean = 0.6446182 for K = 3
                                                                #           $mean = 0.6668141 for K = 4
                                                                #           $mean = 0.6384799 for K = 5
 
# plot DAPC
# use pophelper standard colours
standard_12=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")
 
#jpeg(filename = "b-neutral/DAPCplot_dip_neutral_K2.jpeg", res = 200, width = 1200, height = 1200)
scatter(dapc.clust.best, posi.da="bottomright",  bg="white",
        cstar=1, col=standard_12, scree.pca=TRUE, clab=0,
        posi.pca="topright",
        leg=TRUE, txt.leg=paste("Cluster",levels(grp$grp)), posi.leg = "topleft")
# dev.off()
# dev.set(dev.prev())
 
# barplot
# own barplot, so I can plot by longitude
ind_posterior <- dapc.clust.best$posterior %>% 
  as.data.frame() %>% 
  `colnames<-`(paste0("Cluster",colnames(.))) %>% 
  mutate(IND = rownames(.)) %>% 
  mutate(SamplingCell = gsub("i.*","", IND)) %>% 
  left_join(sampling, by = "SamplingCell") %>% 
  select(-Sampler, -ComplCell, -Area_km2, -Ecoregion, -Ecoregion_adj, -Temp_mean) %>% 
  arrange(Longitude)
# set correct levels order, so that ggplot plots in the order we want
ind_posterior$IND <- factor(ind_posterior$IND, levels = ind_posterior$IND[order(ind_posterior$Longitude)])
# long format for gglot
ggdata_posterior <- melt(ind_posterior, id.vars = c("IND", "SamplingCell", "Longitude", "Latitude"))
 
#stacked barplot
 # pdf(file="b-neutral/DAPCbar_dip_neutral_K2.pdf", width = 16, height = 9)
ggplot(ggdata_posterior) +
  geom_bar(aes(x = IND, y = value, fill = variable), stat = "identity") +
  scale_fill_manual(values = standard_12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) 
 # dev.off()
 # dev.set(dev.prev())
 
# pie map
cell_posterior <- summarise(group_by(ind_posterior, SamplingCell),
                            Cluster1 = mean(Cluster1),
                            Cluster2 = mean(Cluster2),
                            #Cluster3 = mean(Cluster3),
                            #Cluster4 = mean(Cluster4),
                            #Cluster5 = mean(Cluster5),
                            Longitude = first(Longitude),
                            Latitude = first(Latitude))
# plot map 
pie_cell <- cell_posterior[, 2:(ncol(ind_posterior)-3)] %>% 
  data.matrix(rownames.force = NA)
 
lon_cell <- as.numeric(as.vector(cell_posterior$Longitude))
lat_cell <- as.numeric(as.vector(cell_posterior$Latitude))
 
 #pdf(file="b-neutral/DAPCmap_dip_neutral_K2.pdf", width = 16, height = 9)
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47),col = "gray80", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
 #points(membership2_pop$lon, membership2_pop$lat, pch=19, col="red", cex=0.5) 
for(i in 1:nrow(cell_posterior)) {
  floating.pie(lon_cell[i], lat_cell[i], pie_cell[i, ], radius = 0.4, 
              col = standard_12)
}
legend("bottomleft", 
      legend = colnames(pie_cell), 
      title =  "Posterior Membership Probabilities",
      pch=19, cex = 0.8, ncol = 2,
      col= standard_12,
      bty ="o", bg ="gray90",box.col = "gray90")
map.scale(3, 31, ratio=FALSE, relwidth=0.15, cex=1)
map.axes(cex.axis=0.8)
rasterImage(pic$diplodus_sargus_adaptive, 
            xleft = 33, xright = 37, 
            ybottom = 45.1, ytop = 47)
dev.off()  
dev.set(dev.prev())
 
 
 x
### Interpreting variable contributions ----

# which alleles highlight most the separation of the two population clusters?
set.seed(4)
contrib <- loadingplot(dapc1$var.contr, thres = .01, lab.jitter = 1)

# look at allele frequencies of these SNPs
# first convert genlight to genind



