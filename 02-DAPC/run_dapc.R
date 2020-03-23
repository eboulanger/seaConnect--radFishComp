#*****************************************************************************
#******             Detecting population clusters using DAPC            ******
#*****************************************************************************

### tutorial by adegenet
### http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

# load packages
library(adegenet)
library(dplyr)
library(maps)
library(mapdata)
library(plotrix)
library(pophelper)

# define arguments
inputFolderSp <- "01-Diplodus/"
inputFolderLc <- "c-adaptive"
inputFile <- "dip_adaptive_494"

workDir <- paste0("~/Documents/project_SEACONNECT/seaConnect--radFishComp/02-DAPC/",
                  inputFolderSp, inputFolderLc)
inputPath <- paste0("~/Documents/project_SEACONNECT/seaConnect--dataPrep/04-finalPrep/",
                    inputFolderSp, inputFile,".raw")

# set working directory
setwd(workDir)

# load data
genlight <- read.PLINK(inputPath)
sampling <-  read.table("../../../00-Data/coord_seaconnect_tous.txt", sep = "\t", header = T, stringsAsFactors = F)

### Discriminant Analysis of Principal Components -----

## DAPC with sample sites as a priori clusters
pop(genlight) <- gsub("_.*","", pop(genlight))
     # for mullus: need to alter SamplingCell / pop names, nl. add "C" in front.
     pop(genlight) <- paste0("C",pop(genlight))
dapc.pop <- dapc(genlight)
300
2
scatter(dapc.pop)

## DAPC with ecoregions as a priori clusters
sites <- gsub("_.*","", pop(genlight)) %>% 
  as.data.frame() %>% 
  `colnames<-` ("SamplingCell")
ecor <- sampling %>% 
  .[c("SamplingCell", "Ecoregion_adj")] %>% 
  right_join(sites)
pop(genlight) <- ecor$Ecoregion_adj
dapc.pop <- dapc(genlight)
300
2
     jpeg(filename = paste0("DAPCplot_ecor_",inputFile,".jpeg"), res = 200, width = 1200, height = 1200)
scatter(dapc.pop, posi.da = "bottomleft")
     dev.off()

## DAPC without a priori populations

# search for clusters
grp <- find.clusters(genlight, max.n.clust=20)
250
2
     # manually export BIC plot. name: BICplot_inputFile, dimensions: 700 x 600

# run DAPC with new a priori clusters
dapc1 <- dapc(genlight, grp$grp)
250
1

#### DAPC scatterplot ----

# use pophelper standard colours
standard_12=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")

# for two groups
     jpeg(filename = paste0("DAPCplot_K2_",inputFile,".jpeg"), res = 200, width = 1200, height = 1200)
scatter(dapc1,1,1,  bg="white",
        cstar=0, col=standard_12, scree.pca=TRUE, clab=0,
        posi.pca="none",
        leg=TRUE, txt.leg=paste("Cluster",levels(grp$grp)), posi.leg = "topright")
     dev.off()
# more than two groups
     jpeg(filename = paste0("DAPCplot_K3_",inputFile,".jpeg"), res = 200, width = 1200, height = 1200)
scatter(dapc1, posi.da="none",  bg="white",
        cstar=0, col=standard_12, scree.pca=TRUE, clab=0,
        posi.pca="none",
        leg=TRUE, txt.leg=paste("Cluster",levels(grp$grp)), posi.leg = "topleft")
     dev.off()

# barplot
compoplot(dapc1, col=standard_12[1:4], border = NA, show.lab = F)

### DAPC pie map ----
# pie map when K = 2. if K > 2: remember to adjust script

posterior2 <- dapc1$posterior
colnames(posterior2) <- c("Cluster1", "Cluster2")
ID <- rownames(posterior2)
   # for mullus only: alter individual id
   # id1 <- gsub("_.*","", ID)
   # id2 <- gsub(".*_","", ID)
   # Cid <- paste("C", id1, "_", "D", id1, "i", id2, sep ="")
   # ID <- rownames(posterior2) <- Cid
SamplingCell <- gsub("_.*","", ID)
membership2 <- cbind("ID" = ID, "SamplingCell" = SamplingCell) %>% 
  as.data.frame()
membership2$Cluster1 = as.numeric(posterior2[,1])
membership2$Cluster2 = as.numeric(posterior2[,2])
#membership2$Cluster3 = as.numeric(posterior2[,3])
str(membership2)

membership2_cell <- summarise(group_by(membership2, SamplingCell),
                              Cluster1 = mean(Cluster1),
                              Cluster2 = mean(Cluster2)) %>% 
  left_join(sampling[,c("SamplingCell", "Longitude","Latitude","Ecoregion_adj")], by = "SamplingCell") %>% 
  as.data.frame()
  
# plot map 
# create separate datasets for the function floating.pie
pie2_pop <- membership2_cell[, 2:3] %>% 
  data.matrix(rownames.force = NA)
  # slightly alter proportions to avoid wrongfull coloration when using floating.pie
  # problem occurs when the cluster values are 1 or 0
  pie2_pop[pie2_pop == 0] <- 0.0001
  pie2_pop[pie2_pop == 1] <- 0.9999

lon_pop <- as.numeric(as.vector(membership2_cell$Longitude))
lat_pop <- as.numeric(as.vector(membership2_cell$Latitude))

     pdf(file=paste0("DAPCmap_K2_",inputFile,".pdf"), width = 7, height = 5)
carte=map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE)
#points(membership2_pop$lon, membership2_pop$lat, pch=19, col="red", cex=0.5) 
for(i in 1:nrow(membership2_cell)) {
  floating.pie(lon_pop[i], lat_pop[i], pie2_pop[i, ], radius = 0.4, 
               col = c("#2121D9", "#9999FF"))
}
legend("bottomleft", 
       legend = c("Cluster 1", "Cluster 2"), 
       title =  "Posterior Membership \nProbabilities",
       pch=19, cex = 0.8, ncol = 2,
       col=c("#2121D9","#9999FF","#DF0101","#04B404"),
       bty ="o", bg ="gray90",box.col = "gray90")
map.scale(5.7, 31.1,relwidth = 0.15, metric = TRUE, ratio = FALSE, cex = 0.8)
map.axes(cex.axis=0.8)
     dev.off()

### DAPC ancestry barplot ----

ancestry2_coord <- membership2 %>% 
       select("SamplingCell", "ID", everything())
ancestry2_coord <- merge(ancestry2_coord, sampling[,c("SamplingCell","Ecoregion_adj","Longitude","Latitude")])
ancestry2_lon <- ancestry2_coord[order(ancestry2_coord$Longitude),]  # order by longitude
     
# stacked barplot base R
barplot(t(ancestry2_lon[,c("Cluster1","Cluster2")]), col = rainbow(3),
        xlab = "Individuals by longitude", ylab = "Ancestry", 
        border = NA)

# stacked barplot using pophelper 

# transform data as qlist
qlist <- list("2_clusters"=ancestry2_lon[,3:4])
qlist <- lapply(qlist,"rownames<-",ancestry2_lon$ID) # add individual labels
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
      splab = paste0(inputFile,"\nK = 2"), # label to appear on the right side of plot
      splabsize = 3, divsize = 0.1, 
      showlegend = T, legendtextsize = 2, legendkeysize = 2,
      ordergrp = TRUE, subsetgrp=c("Alboran\nSea", "Western\nMediterranean","Central\nMediterranean","Adriatic\nSea","Ionian\nSea","Tunisian\nPlateau","Aegean\nSea", "Levantine\nSea"),
      grplabangle=0, imgtype = "pdf",
      outputfilename = paste0("DAPCbarplot_",inputFile,"_K2"))

### export DAPC results ----
# create new populations
membership2$Pop <- as.numeric(membership2$Cluster1)
membership2$Pop[membership2$Pop > 0.99] <- "DAPCluster1"
membership2$Pop[membership2$Pop < 0.01] <- "DAPCluster2"
membership2$Pop[membership2$Pop == 1.65575725825131e-318] <- "DAPCluster2"

write.table(membership2, file = paste0("DAPC_newPop_", inputFile,".csv"), quote = F, sep = ";")
### Interpreting variable contributions ----

# which alleles highlight most the separation of the two population clusters?
set.seed(4)
contrib <- loadingplot(dapc1$var.contr, thres = .005, lab.jitter = 1)

# look at allele frequencies of these SNPs
# first convert genlight to genind



