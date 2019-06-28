#*****************************************************************************
#******             Detecting population clusters using DAPC            ******
#*****************************************************************************

### tutorial by adegenet
### http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

setwd("/Users/eboulanger-admin/Documents/project_SEACONNECT/seaConnect--radFishComp/02-DAPC/01-Diplodus/c-adaptive/")

# load packages
library(adegenet)
library(dplyr)

# load data
all <-"/Users/eboulanger-admin/Documents/project_SEACONNECT/seaConnect--radFishComp/00-Data/01-Diplodus/a-all/diplodus.raw" 
adaptive <-"/Users/eboulanger-admin/Documents/project_SEACONNECT/seaConnect--radFishComp/00-Data/01-Diplodus/c-adaptive/sar_adaptive.raw" 
dip_adpt <- read.PLINK(adaptive)

sampling <-  read.table("../../../00-Data/coord_seaconnect_tous.txt", sep = "\t", header = T, stringsAsFactors = F)

### Discriminant Analysis of Principal Components -----

## DAPC with sample sites as a priori clusters
pop(dip_adpt) <- gsub("_.*","", pop(dip_adpt))
dapc.dip.pop <- dapc(dip_adpt)
50
2
scatter(dapc.dip.pop)

## DAPC with ecoregions as a priori clusters
sites <- gsub("_.*","", pop(dip_adpt)) %>% 
  as.data.frame() %>% 
  `colnames<-` ("SamplingCell")
ecor <- sampling %>% 
  .[c("SamplingCell", "Ecoregion_adj")] %>% 
  right_join(sites)
pop(dip_adpt) <- ecor$Ecoregion_adj
dapc.dip.pop <- dapc(dip_adpt)
50
2
     jpeg(filename = "DAPCplot_ecor.jpeg", res = 200, width = 1200, height = 1200)
scatter(dapc.dip.pop, posi.da = "topright")
     dev.off()

## DAPC without a priori populations

# search for clusters
grp <- find.clusters(dip_adpt, max.n.clust=20)
350
2
# run DAPC with new a priori clusters
dapc1 <- dapc(dip_adpt, grp$grp)
50
1

# plot DAPC
# use pophelper standard colours
standard_12=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")

# for one discriminant function
     jpeg(filename = "DAPCplot_K2.jpeg", res = 200, width = 1200, height = 1200)
scatter(dapc1,1,1,  bg="white",
        cstar=0, col=standard_12, scree.pca=TRUE, clab=0,
        posi.pca="none",
        leg=TRUE, txt.leg=paste("Cluster",1:2), posi.leg = "topright")
     dev.off()
# more than one discriminant function
scatter(dapc1, posi.da="bottomright",  bg="white",
        cstar=0, col=standard_12, scree.pca=TRUE, clab=0,
        posi.pca="topright",
        leg=TRUE, txt.leg=paste("Cluster",levels(grp$grp)), posi.leg = "topleft")

# barplot
compoplot(dapc1, col=standard_12[1:2], border = NA, show.lab = F)

# pie map

library(maps)
library(mapdata)
library(plotrix)

posterior2 <- dapc1$posterior
colnames(posterior2) <- c("Cluster1", "Cluster2")
ID <- rownames(posterior2)
SamplingCell <- gsub("_.*","", ID)
membership2 <- cbind("ID" = ID, "SamplingCell" = SamplingCell) %>% 
  as.data.frame()
membership2$Cluster1 = as.numeric(posterior2[,1])
membership2$Cluster2 = as.numeric(posterior2[,2])
str(membership2)

membership2_cell <- summarise(group_by(membership2, SamplingCell),
                              Cluster1 = mean(Cluster1),
                              Cluster2 = mean(Cluster2)) %>% 
  left_join(sampling[,c("SamplingCell", "Longitude","Latitude","Ecoregion_adj")], by = "SamplingCell") %>% 
  as.data.frame()
  
# plot map 
pie2_pop <- membership2_cell[, 2:3] %>% 
  data.matrix(rownames.force = NA)

lon_pop <- as.numeric(as.vector(membership2_cell$Longitude))
lat_pop <- as.numeric(as.vector(membership2_cell$Latitude))

     pdf(file="DAPCmap_dip_adpt_K2.pdf", width = 7, height = 5)
carte=map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE)
#points(membership2_pop$lon, membership2_pop$lat, pch=19, col="red", cex=0.5) 
for(i in 1:nrow(membership2_cell)) {
  floating.pie(lon_pop[i], lat_pop[i], pie2_pop[i, ], radius = 0.4, 
               col = c("#2121D9", "#9999FF"))
}
legend("bottomleft", 
       legend = c("Cluster 1", "Cluster2"), 
       title =  "Posterior Membership Probabilities",
       pch=19, cex = 0.8, ncol = 2,
       col=c("#2121D9","#9999FF"),
       bty ="o", bg ="gray90",box.col = "gray90")
map.scale(5.7, 31.1,relwidth = 0.15, metric = TRUE, ratio = FALSE, cex = 0.8)
map.axes(cex.axis=0.8)
     dev.off()

  

### Interpreting variable contributions ----

# which alleles highlight most the separation of the two population clusters?
set.seed(4)
contrib <- loadingplot(dapc1$var.contr, thres = .01, lab.jitter = 1)

# look at allele frequencies of these SNPs
# first convert genlight to genind



