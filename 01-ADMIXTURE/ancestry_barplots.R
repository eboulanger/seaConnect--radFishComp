##### cross-validation plot ADMIXTURE -----
getwd()
setwd("/Users/eboulanger-admin/Documents/project_SEACONNECT/seaConnect--radFishComp/01-ADMIXTURE/01-Diplodus/c-adaptive")

library(ggplot2)
library(reshape2)
crossval <- read.table("cross_validation.txt")
K <- c(1: nrow(crossval))
crossval_gg <- cbind(crossval, K) 
pdf("cross_validation.pdf", width = 5, height = 5)
ggplot(data = crossval_gg, aes(x = K, y = V4)) +
  geom_line() +
  geom_point() +
  ylab("Cross-validation error")
dev.off()  

##### stacked barplot ADMIXTURE -----

library(pophelper)
library(gridExtra)

## import ancestry proportions
ancestry2 <- read.table("sar_adaptive.2.Q", col.names = c("Cluster1","Cluster2"))
ancestry3 <- read.table("sar_adaptive.3.Q", col.names = c("Cluster1","Cluster2","Cluster3"))

ancestry3 <- read.table("dip_adaptive_wrongid.3.Q", col.names = c("Cluster1","Cluster2","Cluster3"))

# import necessary metadata
id <- read.table("id_admixture.txt", col.names = "ID",stringsAsFactors = F)
id$SamplingCell<- gsub("_.*","", id$ID)
#coord <- read.table("../../../00-Data/01-Diplodus/coord_sar_297Ind.txt")
sampling <- read.table("../../../00-Data/coord_seaconnect_tous.txt", sep = "\t", header = T,
                       strip.white = T, stringsAsFactors = F)
#sampling$Ecoregion[sampling$Ecoregion == "Saharan Upwelling"] <- "Alboran Sea"

ancestry2_coord <- cbind(id, ancestry2)
ancestry2_coord <- merge(ancestry2_coord, sampling[,c("SamplingCell","Ecoregion_adj","Longitude","Latitude")])
ancestry2_lon <- ancestry2_coord[order(ancestry2_coord$Longitude),]  # order by longitude

ancestry3_coord <- cbind(id, ancestry3)
ancestry3_coord <- merge(ancestry3_coord, sampling[,c("SamplingCell","Ecoregion_adj","Longitude","Latitude")])
ancestry3_lon <- ancestry3_coord[order(ancestry3_coord$Longitude),]  # order by longitude

# stacked barplot base R
barplot(t(ancestry2_lon[,c("Cluster1","Cluster2")]), col = rainbow(3),
        xlab = "Individuals by longitude", ylab = "Ancestry", 
        border = NA)

# stacked barplot using pophelper -----

# read in / transform data as qlist
# make own qlist from dataframes
qlist <- list("sar_adaptive.2.Q"=ancestry2_lon[,3:4],"sar_adaptive.3.Q"=ancestry3_lon[,3:5])
qlist <- lapply(qlist,"rownames<-",ancestry3_lon$ID) # add individual labels
attributes(qlist)

# add grouping labels
label <- as.data.frame(ancestry2_lon$Ecoregion_adj, stringsAsFactors = F)
colnames(label) <- "Ecoregion"
str(label)
label2 <- as.data.frame(gsub(" ","\n", label[,1]), stringsAsFactors = F) # labels will appear on two lines
colnames(label2) <- "Ecoregion"
str(label2)

# Q plot
# only K = x
plotQ(qlist[1], returnplot = F, exportplot = T, showindlab = F,
      grplab=label2, grplabsize=1,linesize=0.1,pointsize=1, grplabpos = 0.4, 
      splabsize = 3, divsize = 0.1, 
      showlegend = T, legendtextsize = 2, legendkeysize = 2,
      ordergrp = TRUE, subsetgrp=c("Alboran\nSea", "Western\nMediterranean","Central\nMediterranean","Adriatic\nSea","Ionian\nSea","Tunisian\nPlateau","Aegean\nSea", "Levantine\nSea"),
      grplabangle=0, imgtype = "pdf")

# both K = 2 and K = 3
plotQ(qlist[1:2], returnplot = F, exportplot = T, showindlab = F,
      imgoutput="join",
      grplab=label2, grplabsize=1,linesize=0.1,pointsize=1, grplabpos = 0.3, 
      splabsize = 3, divsize = 0.1, 
      showlegend = T,legendtextsize = 2, legendkeysize = 2,
      ordergrp = TRUE, subsetgrp=c("Alboran\nSea", "Western\nMediterranean","Central\nMediterranean","Adriatic\nSea","Ionian\nSea","Tunisian\nPlateau","Aegean\nSea", "Levantine\nSea"),
      grplabangle=0, imgtype = "pdf")

##### ancestry pie maps ----
library(maps)
library(mapdata)

# get mean ancestry by sampling cell and add geographic data
ancestry2_cell <- summarise(group_by(ancestry2_lon, SamplingCell), 
                            Cluster1 = mean(Cluster1),
                            Cluster2 = mean(Cluster2)) %>% 
  left_join(sampling[,c("SamplingCell", "Longitude","Latitude","Ecoregion")], by = "SamplingCell") %>% 
  as.data.frame()

# plot map 
pie2_pop <- ancestry2_cell[, 2:3] %>% 
  data.matrix(rownames.force = NA)

lon_pop <- as.numeric(as.vector(ancestry2_cell$Longitude))
lat_pop <- as.numeric(as.vector(ancestry2_cell$Latitude))

pdf(file="sar_adaptive.2.Q.piemap.pdf", width = 7, height = 5)
carte=map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE)
#points(membership2_pop$lon, membership2_pop$lat, pch=19, col="red", cex=0.5) 
for(i in 1:nrow(ancestry2_cell)) {
  floating.pie(lon_pop[i], lat_pop[i], pie2_pop[i, ], radius = 0.4, 
               col = c("#2121D9","#9999FF"))
}
legend("bottomleft", 
       legend = c("Cluster 1", "Cluster2"), 
       title = "Ancestry fractions",
       pch=19, cex = 0.8, ncol = 2,
       col=c("#2121D9","#9999FF"),
       bty ="o", bg ="gray90",box.col = "gray90")
map.scale(5.7, 31.1,relwidth = 0.15, metric = TRUE, ratio = FALSE, cex = 0.8)
map.axes(cex.axis=0.8)
dev.off()




x
### reserve ####

subsetgrp=c("Saharan Upwelling","Alboran Sea", "Western Mediterranean","Adriatic Sea","Ionian Sea","Tunisian Plateau","Aegean Sea", "Levantine Sea"),

p2 <- plotQ(q3, returnplot = T, exportplot = F, showindlab = F,
            grplab=label2, grplabsize=2,linesize=0.3,pointsize=3, grplabpos = 0.5, splabsize = 2,
            ordergrp = TRUE, subsetgrp=c("Saharan\nUpwelling","Alboran\nSea", "Western\nMediterranean","Adriatic\nSea","Ionian\nSea","Tunisian\nPlateau","Aegean\nSea", "Levantine\nSea"),
            grplabangle=20, grplabheight = 10)
grid.arrange(p2$plot[[1]])

pdf("test.pdf", width = 7, height = 3)
grid.arrange(p2$plot[[1]])
dev.off()

