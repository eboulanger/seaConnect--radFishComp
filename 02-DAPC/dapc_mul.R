#*****************************************************************************
#******             Detecting population clusters using DAPC            ******
#*****************************************************************************

### tutorial by adegenet: http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

setwd("/Users/eboulanger-admin/Documents/project_SEACONNECT/seaConnect--radFishComp/02-DAPC/02-Mullus")

# load packages
library(adegenet)
library(dplyr)
library(parallel)
library(png)
library(stringr)
library(tidyr)
library(patchwork)
 # for bar chart
library(reshape)
library(ggplot2)
 # for pie map
library(maps)
library(mapdata)
library(plotrix)


# set variables
sp <- "mul"
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
  adaptive <-"../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_adaptive_291.raw" 
  neutral  <-"../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_neutral_2462.raw"
} else {
  all      <-"../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/dip_all_filtered.raw" 
  adaptive <-"../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/dip_adaptive_413.raw" 
  neutral  <-"../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/dip_neutral_7655.raw"
}
  # read plink file
  dat <- read.PLINK(get(subset))

# alternatively, if don't know how many markers, use regex!
# 1 specify path to file
if (sp == "mul") {
  path = "../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/"
} else {
  path = "../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/"
}
# 2 get file using regex
datafile <- list.files(path = path, pattern = paste0(sp,"_",subset,"_[0-9]+.raw"))

# read plink file
dat <- read.PLINK(paste0(path, datafile))

### Discriminant Analysis of Principal Components ###

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
a.score(dapc.ecor) # adaptive: for N PCA = 400, mean a-score =  0.01149647
                   #           which is very low, and indicates that the result is very unstable
                   # neutral: 0.05590811
# find the optimal number of PC's to retain based on the a-score
dapc.ecor.optimscore <- optim.a.score(dapc.ecor) # adaptive: $1
                                                 # neutral : $19
# run the dapc with the optimal number of PCs
dapc.ecor.best <- dapc(dat, n.pca = dapc.ecor.optimscore$best, n.da = 2)
dapc.ecor.best.score <- a.score(dapc.ecor.best, n.sim = 1000) # adaptive: $mean = 0.2095513
                                                              # neutral : $mean =       

     #jpeg(filename = "c-adaptive/DAPCplot_mul_adaptive_ecor_optim.jpeg", res = 200, width = 1200, height = 1200)
scatter(dapc.ecor.best)
title("DAPC M. surmuletus adaptive\necoregions a priori")
     dev.off()
     dev.set(dev.prev())

## DAPC without a priori populations ----

# search for clusters
grp <- find.clusters(dat, max.n.clust=20, n.pca = round(nrow(dat)/3)) # n.pca = 1/3 n.ind
# run DAPC with new a priori clusters
dapc.clust <- dapc(dat, grp$grp, n.pca = nrow(dat)/3, n.da = 2)
scatter(dapc.clust)
a.score(dapc.clust, n.sim = 100)
# find the optimal number of PC's to retain based on the a-score
dapc.clust.optimscore <- optim.a.score(dapc.clust) # neutral : $66 for K = 2
                                                   # adaptive: $33 for K = 7, $1 for K = 3, $27 for K=5
# run the dapc with the optimal number of PCs
dapc.clust.best <- dapc(dat, grp$grp, n.pca = dapc.clust.optimscore$best, n.da=2)
scatter(dapc.clust.best)
dapc.clust.best.score <- a.score(dapc.clust.best, n.sim = 1000) # neutral : $mean = 0.3741986 for K = 2, n.pc = 75
                                                                # adaptive: $mean = 0.4777093 for K = 2, n.pc = 5
                                                                #           $mean = 0.609534 for K = 3
                                                                #           $mean = 0.6706575 for K = 4
                                                                #           $mean = 0.6384799 for K = 5
     # write Rdata file for neutral DAPC results, because take a long time to run
     saveRDS(dapc.clust.best, file = "b-neutral/DAPC_mul_neutral_K2.rds")
     
     # read in Rdata to not have to run dapc again
     dapc.clust.best <- readRDS("c-adaptive/DAPC_mul_adaptive_K4.rds")
     
     # create indmap for later analyses
     indmap <- as.data.frame(dapc.clust.best$grp)
     indmap$INDIVIDUALS <- rownames(indmap)
     indmap$STRATA <- paste0("Cluster", indmap$`dapc.clust.best$grp`)
     indmap <- indmap[,c("INDIVIDUALS", "STRATA")]
     head(indmap)
     unique(indmap$STRATA)
     #write.table(indmap, file = "c-adaptive/DAPCgrp_mul_adaptive_K3.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# plot DAPC
# use pophelper standard colours
standard_12=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")
# colourblind palette
   library(ggthemes)
   library(scales)
   show_col(colorblind_pal()(8))
   colblind_3 <- c("#0072B2","#009E73","#D55E00") # 3 clusters: group1 = west med, group 2 = gibraltar, group 3 = east med
   colblind_4 <- c("#0072B2","#E69F00","#D55E00", "#009E73") # 4 clusters: group1 = west med 1, group2 = west med 2, group 3 = east med, group 4 = alboran
   
colblind_5 <- c("#0072B2","#009E73","#D55E00", "#E69F00", "#F0E442")
colblind_7 <- c("#0072B2","#009E73","#D55E00", "#E69F00", "#F0E442","#CC79A7","#56B4E9")


#      #jpeg(filename = "c-adaptive/DAPCplot_mul_adaptive_K3.jpeg", res = 200, width = 1200, height = 1200)
# scatter(dapc.clust.best, posi.da="bottomright",  bg="white",
#         cstar=1, col=colblind_7, scree.pca=FALSE, clab=0,
#         posi.pca="bottomleft",
#         leg=TRUE, txt.leg=paste("Cluster",levels(dapc.clust.best$grp)), posi.leg = "topleft")
#      # dev.off()
#      # dev.set(dev.prev())

# with ggplot
ggdata.dapc <- as.data.frame(dapc.clust.best$ind.coord)
ggdata.dapc$group <- dapc.clust.best$grp
ggdata.dapc$ind <- rownames(ggdata.dapc)
centroid <- as.data.frame(dapc.clust.best$grp.coord)
colnames(centroid) <- paste0(colnames(centroid), ".ctr")
centroid$group <- rownames(centroid)
ggdata.dapc.ctr <- plyr::join(ggdata.dapc, centroid)

scatter1 <- ggplot(ggdata.dapc.ctr, aes(x = LD1, y = LD2, color = group, fill = group)) + 
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_point(size = 3, shape = 21) + 
  scale_color_manual(values=c(colblind_4)) + 
  scale_fill_manual(values=c(colblind_4)) +
  labs(x= "DPC1", y="DPC2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect()) 
# add stars and circle
dapc.scatter <- scatter1 + geom_segment(aes(x = LD1, y = LD2, xend = LD1.ctr, yend = LD2.ctr, colour = group)) +
  geom_label(data = unique(ggdata.dapc.ctr[,c("group", "LD1.ctr", "LD2.ctr")]),
             aes(x = LD1.ctr, y = LD2.ctr, label = group), fill = "white") +
  guides(colour = "none", fill = "none")

# K=3: only one DA so density plot
dapc.density <- ggplot(ggdata.dapc.ctr, aes(x = LD1, color = group, fill = group)) +
  geom_density() +
  scale_color_manual(values=c(colblind_2)) + 
  scale_fill_manual(values=c(colblind_2)) +
  labs(x= "DPC1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect()) +
  geom_label(data = unique(ggdata.dapc.ctr[,c("group", "LD1.ctr")]),
             aes(x = LD1.ctr, y = 0.1, label = group), fill = "white") +
  guides(colour = "none", fill = "none")

# barplot
  #pdf("b-neutral/DAPCcompoplot_mul_neutral_k3.pdf", width = 10, height = 6)
compoplot(dapc.clust.best, col=standard_12, border = NA, show.lab = F)
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

ind_grp <- select(ggdata.dapc, ind, group) %>% 
  mutate(SamplingCell = gsub("i.*","", ind)) %>% 
  mutate(group = paste0("group",group)) %>% 
  left_join(select(sampling, SamplingCell, Longitude, Latitude), by = "SamplingCell") %>% 
  mutate(value = rep(1, nrow(.))) %>% 
  pivot_wider(names_from = group, values_from = value) %>% 
  replace(is.na(.),0) %>% 
  arrange(Longitude)
# set correct levels order, so that ggplot plots in the order we want
ind_grp$ind <- factor(ind_grp$ind, levels = ind_grp$ind[order(ind_grp$Longitude)])
# long format
ggdata_grp <- pivot_longer(ind_grp, cols = starts_with("group"), names_to = "group")

#stacked barplot
# pdf(file="c-adaptive/DAPCbar_mul_adaptive_K3.pdf", width = 16, height = 9)
ggplot(ggdata_posterior) +
  geom_bar(aes(x = IND, y = value, fill = variable), stat = "identity") +
  scale_fill_manual(values = standard_12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) 
 # dev.off()
 # dev.set(dev.prev())

# pdf(file="c-adaptive/DAPCbar_mul_adaptive_K5_grp.pdf", width = 16, height = 9)
grp.bar <- ggplot(ggdata_grp) +
  geom_bar(aes(x = ind, y = value, fill = group), stat = "identity")  +
  scale_fill_manual(values = colblind_4) +
  guides(fill = "none") +
  labs(x = "Individuals by longitude", y = "") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
 
combi1 <- dapc.density / grp.bar + plot_layout(heights = c(2,1))
ggsave(combi1, file = "c-adaptive/DAPC_scatter_bar_K3.pdf")
combi1 <- dapc.scatter / grp.bar + plot_layout(heights = c(2,1))
ggsave(combi1, file = "c-adaptive/DAPC_scatter_bar_K4.pdf")

# pie map
cell_posterior <- dplyr::summarise(group_by(ind_posterior, SamplingCell),
                            Cluster1 = mean(Cluster1),
                            Cluster2 = mean(Cluster2),
                            Cluster3 = mean(Cluster3),
                            #Cluster4 = mean(Cluster4),
                            #Cluster5 = mean(Cluster5),
                            #Cluster6 = mean(Cluster6),
                            #Cluster7 = mean(Cluster7),
                            Longitude = first(Longitude),
                            Latitude = first(Latitude))
cell_grp <- dplyr::summarise(group_by(ind_grp, SamplingCell),
                      group1 = mean(group1),
                      group2 = mean(group2),
                      group3 = mean(group3),
                      #group4 = mean(group4),
                      #group5 = mean(group5),
                      #group6 = mean(group6),
                      #group7 = mean(group7),
                      Longitude = first(Longitude),
                      Latitude = first(Latitude))

# plot map 
pie_cell <- cell_grp[, 2:(ncol(ind_grp)-3)] %>% 
  data.matrix(rownames.force = NA)
pie_cell[pie_cell == 0] <- 0.000001 # otherwise 0 values are skipped when assigning colours

lon_cell <- as.numeric(as.vector(cell_grp$Longitude))
lat_cell <- as.numeric(as.vector(cell_grp$Latitude))

  #pdf(file="c-adaptive/DAPCgrpmap_mul_adaptive_K3.pdf", width = 16, height = 9)
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47),col = "gray80", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
#points(membership2_pop$lon, membership2_pop$lat, pch=19, col="red", cex=0.5) 
for(i in 1:nrow(cell_grp)) {
  floating.pie(lon_cell[i], lat_cell[i], pie_cell[i, ], radius = 0.4, 
               col = colblind_7) }
legend("bottomleft", 
       legend = colnames(pie_cell), 
       #title =  "Posterior Membership Probabilities",
       title = "Group membership",
       pch=19, cex = 1, ncol = 2,
       col = colblind_7,
       bty ="o", bg ="gray90",box.col = "gray90")
map.scale(3, 31, ratio=FALSE, relwidth=0.15, cex=1)
map.axes(cex.axis=0.8)
#rasterImage(pic$mullus_surmuletus_adaptive, 
#            xleft = 33, xright = 37, 
#            ybottom = 45.1, ytop = 47)
    dev.off()  
    dev.set(dev.prev())

# regular map
ind.grp.coord <- select(ggdata.dapc, ind, group) %>% 
  mutate(SamplingCell = gsub("i.*","", ind)) %>% 
  mutate(group = paste0("group",group)) %>% 
  left_join(select(sampling, SamplingCell, Longitude, Latitude), by = "SamplingCell") 
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47),col = "gray80", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
points(x=ind.grp.coord$Longitude, y=ind.grp.coord$Latitude, col = )
library(basicPlotteR)
temp <- addPoints(ind.grp.coord$Longitude, ind.grp.coord$Latitude, cex=1, col.line="red")

    
    x
    
### extract individuals ####
    # individuals from cluster 6 consistently cluster together. see if same as individuals that cluster separately in RDA
outl.ind.mul <- dapc.clust.best$posterior %>% 
      as.data.frame() %>% 
      `colnames<-`(paste0("Cluster",colnames(.))) %>% 
      mutate(ind = rownames(.)) %>% # filter drops rownames so add individual identifiers
      filter(Cluster6 > 0.05) %>% 
      select(ind) %>% 
      pull()
    #    "C47i4"  "C47i5"  "C47i6"  "C47i7"  "C47i8"  "C80i10" "C80i1"  "C80i3"  "C80i4" 
    #    "C80i5"  "C80i6"  "C80i7"  "C90i5"  "C90i6" 
    #  manually
    outl.ind.mul <- c("C47i4","C47i5","C47i6","C47i7","C47i8","C80i10","C80i1","C80i3","C80i4",
                      "C80i5","C80i6","C80i7","C90i5","C90i6")
### iter 2: remove outlier clustered individuals and run dapc ####
dat.sub <- dat[!indNames(dat) %in% outl.ind.mul]

## DAPC without a priori populations
# search for clusters
grp <- find.clusters(dat.sub, max.n.clust=20, n.pca = round(nrow(dat)/3), n.clust = 3) # n.pca = 1/3 n.ind
# run DAPC with new a priori clusters
dapc.sub.clust <- dapc(dat.sub, grp$grp, n.pca = nrow(dat)/3, n.da = 2)
scatter(dapc.sub.clust)
a.score(dapc.sub.clust, n.sim = 100)
# plot DAPC
# use pophelper standard colours
# standard_12=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")
#jpeg(filename = "b-neutral/DAPCplot_mul_neutral_sub_K3.jpeg", res = 200, width = 1200, height = 1200)
scatter(dapc.sub.clust, posi.da="topleft",  bg="white",
        cstar=1, col=standard_12, scree.pca=TRUE, clab=0,
        posi.pca="topright",
        leg=TRUE, txt.leg=paste("Cluster",levels(grp$grp)), posi.leg = "bottomright")
# dev.off()
# dev.set(dev.prev())
# own barplot, so I can plot by longitude
ind_posterior <- dapc.sub.clust$posterior %>% 
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
# pdf(file="b-neutral/DAPCbar_mul_neutral_sub_K3.pdf", width = 16, height = 9)
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
                            Cluster3 = mean(Cluster3),
                            Cluster4 = mean(Cluster4),
                            Cluster5 = mean(Cluster5),
                            Cluster6 = mean(Cluster6),
                            Longitude = first(Longitude),
                            Latitude = first(Latitude))

# plot map 
pie_cell <- cell_posterior[, 2:(ncol(ind_posterior)-3)] %>% 
  data.matrix(rownames.force = NA)

lon_cell <- as.numeric(as.vector(cell_posterior$Longitude))
lat_cell <- as.numeric(as.vector(cell_posterior$Latitude))

#pdf(file="b-neutral/DAPCmap_mul_neutral_sub_K3.pdf", width = 16, height = 9)
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
rasterImage(pic$mullus_surmuletus_adaptive, 
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



