# create maps of sampling cell with size indicating sample size

# libraries
library(png)
library(maps)
library(mapdata)
library(dplyr)
library(ggplot2)
library(scales)

library(stringr)
library(reshape)

source("scale_bar.R")
# import stamps 
data_pic <- list.files(path = "../../00-Misc/stamps/", pattern="*.png",recursive = FALSE)
source_pic <- paste0("../../00-Misc/stamps/", data_pic)
pic <- lapply(source_pic,readPNG)
names(pic) <- str_sub(data_pic,1, -5)

# coord data
coord <- read.table("data/coord_seaconnect_tous.txt", sep = "\t", head = TRUE) 
pop_dip <- read.table("data/dip_population_map_297ind.txt", sep = "\t", head = TRUE)
pop_mul <- read.table("data/mul_population_map_467ind.txt", sep = "\t", head = TRUE)
# wrangle to one dataset with coords and sample size diplodus and mullus
n_dip <-pop_dip %>% 
  group_by(STRATA) %>% 
  summarise(length(INDIVIDUALS))
colnames(n_dip)[2] <- "n_dip"
n_mul <-pop_mul %>% 
  group_by(STRATA) %>% 
  summarise(length(INDIVIDUALS))
colnames(n_mul)[2] <- "n_mul"

coord_size <- coord %>% 
  mutate(STRATA = SamplingCell) %>% 
  select(STRATA, Longitude, Latitude) %>% 
  left_join(n_dip, by = "STRATA") %>% 
  left_join(n_mul, by = "STRATA") 
summary(coord_size$n_dip)
summary(coord_size$n_mul)

# how many sites in total? remove rows with twice NA
coord_size[rowSums(is.na(coord_size)) != 2,] %>% nrow()

# separate maps ----
# map D sargus sampling
#pdf(file="sampling_map_diplodus_297ind_cellNum.pdf", width = 16, height = 9)
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47),col = "gray80", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
points(coord_size$Longitude, coord_size$Latitude, pch=19, col="#053061", cex=scales::rescale(coord_size$n_dip,to = c(1,4))) 
legend("bottomleft", 
       legend = c(1, 3, 5, 7, 10), 
       title = "# Diplodus sargus",
       pch=20,col="#053061",cex=1.2,
       pt.cex = c(1, 1.7, 2.5, 3.3, 4),
       ncol=3,
       bg = "transparent", bty = "n")
map.axes(cex.axis=1)
map.scale(3, 31, ratio=FALSE, relwidth=0.15, cex=1)
text(labels = coord_size$STRATA[!is.na(coord_size$n_dip)], 
     coord_size$Longitude[!is.na(coord_size$n_dip)], 
     coord_size$Latitude[!is.na(coord_size$n_dip)] + 0.5, cex = 0.7)


# map M surmuletus sampling
#pdf(file="sampling_map_mullus_424ind_cellNum.pdf", width = 16, height = 9)
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47),col = "gray80", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
points(coord_size$Longitude, coord_size$Latitude, pch=19, col="#67001f", cex=rescale(coord_size$n_mul,to = c(1, 4)))
legend("bottomleft", 
       legend = c(1, 3, 5, 7, 10), 
       title = "# Mullus surmuletus",
       pch=20,col="#67001f",cex=1.2,
       pt.cex = c(1, 1.7, 2.5, 3.3, 4),
       ncol=3,
       bg = "transparent", bty = "n")
map.axes(cex.axis=1)
map.scale(3, 31, ratio=FALSE, relwidth=0.15, cex=1)
text(labels = coord_size$STRATA[!is.na(coord_size$n_mul)], 
     coord_size$Longitude[!is.na(coord_size$n_mul)], 
     coord_size$Latitude[!is.na(coord_size$n_mul)] + 0.5, cex = 0.7)

# both with different symbols
pdf(file="maps/sampling_map_both_triangle.pdf", width = 16, height = 9)
map("world", xlim=c(-8,37), ylim=c(29.5,47),col = "gray80", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
points(coord_size$Longitude, coord_size$Latitude, pch=2, col="#67001f", cex=rescale(coord_size$n_mul,to = c(1, 4)))
points(coord_size$Longitude, coord_size$Latitude, pch=6, col="#053061", cex=scales::rescale(coord_size$n_dip,to = c(1,4))) 
legend("bottomleft", 
       legend = c(1, 3, 5, 7, 10), 
       title = "# individuals",
       pch=11,col="black",cex=1.2,
       pt.cex = c(1, 1.7, 2.5, 3.3, 4),
       ncol=3,
       bg = "transparent", bty = "n")
legend("topleft", 
       legend = c("Diplodus sargus", "Mullus surmuletus"),
       pch=c(6,2),col=c("#053061","#67001f"),cex=1.2,
       pt.cex = 1.2,
       bg = "transparent", bty = "n")
dev.off()


# combined map ----
# two species on one map: pie charts where both present
map_data <- coord_size[,c("STRATA", "Longitude", "Latitude")]
map_data$diplodus <- coord_size$n_dip
map_data$mullus   <- coord_size$n_mul
map_data$diplodus[is.na(map_data$diplodus)] <- 0
map_data$diplodus[map_data$diplodus>0] <- 1
map_data$mullus[is.na(map_data$mullus)] <- 0
map_data$mullus[map_data$mullus>0] <- 1
# remove empty rows
map_data <- map_data[rowSums(map_data == 0) != 2,]

##### pie maps #####
pie_cell <- map_data[,c("diplodus", "mullus")] %>% 
  data.matrix(rownames.force = NA)
# replace 0's with 0.00001 so they are not ignored
pie_cell[pie_cell == 0] <- 0.0001
lon_cell <- map_data$Longitude
lat_cell <- map_data$Latitude


#pdf(file="sampling_map_both_cellNum.pdf", width = 16, height = 9)
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47),col = "gray80", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
for(i in 1:nrow(pie_cell)) {
  floating.pie(lon_cell[i], lat_cell[i], pie_cell[i, ], radius = 0.4, col = c("#4393c3", "#d6604d"))
}
legend("bottomleft", 
       legend = c("Diplodus sargus", "Mullus surmuletus"), 
       pch=19, cex = 1.5, ncol = 1,
       col= c("#4393c3", "#d6604d"),
       bty ="o", bg ="gray90",box.col = "gray90")
map.scale(3, 31, ratio=FALSE, relwidth=0.15, cex=1)
map.axes(cex.axis=0.8)
#rasterImage(pic$diplodus_sargus, 
#            xleft = 2, xright = 3.7, 
#            ybottom = 31, ytop = 32)
text(labels = map_data$STRATA, 
     map_data$Longitude, 
     map_data$Latitude+ 0.5, cex = 0.7)
dev.off()  
dev.set(dev.prev())

##### symbol map #####

map_data$species <- rep(0, nrow(map_data))
map_data$species[map_data$diplodus == 1 & map_data$mullus == 1] <- "both"
map_data$species[map_data$diplodus == 1 & map_data$mullus == 0] <- "diplodus"
map_data$species[map_data$diplodus == 0 & map_data$mullus == 1] <- "mullus"

# map in ggplot2 ----
# tutorial: http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html

wH <- map_data("world",  xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea
# further subset dataset so don't plot whole polygons
# wH_sub <- wH[wH$long<37 & wH$long > c(-8) & wH$lat < 47 & wH$lat > 29.5,] # creates weird margins around map. rather set limits in ggplot

med_base <- ggplot() + 
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = "black") + 
  coord_fixed( xlim=c(-8,37), ylim=c(29.5,46.5), ratio = 1.3) +
  labs(x= "Longitude (°)", y = "Latitude (°)") +
  #theme_nothing() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))  # add black border on top again

sampling <- med_base +
  geom_point(data = map_data, aes(x=Longitude, y = Latitude, shape = species), cex = 3) +
  geom_text(data = map_data, aes(x=Longitude, y=Latitude + 0.4, label = STRATA), cex = 3) +
  scale_shape_discrete(labels=c("both species", "Diplodus sargus", "Mullus surmuletus")) +
  theme(legend.position = c(0.07, 0.07),
        legend.background = element_blank(),
        legend.title = element_blank())
#sampling
ggsave(sampling, filename= "maps/sampling_map_both_cellNum_shape.pdf", width = 13, height = 7)
# export map data for adding shapes to other figures
write.csv(map_data, "sampling_species_data.csv", row.names = F)

# add fish silhouettes
# install EBImage
library(ggimage)

med_base +
  geom_image(data = map_data, aes(x=Longitude, y = Latitude), 
             image="../../00-Misc/stamps/diplodus_sargus.png", size = 0.05) 


#### combined but overlapping ####
med_base +
  geom_point(data = select(coord_size, -n_mul), aes(x= Longitude, y= Latitude, size = n_dip), pch = 21,col = "blue") + #, col = "white", alpha = 0.5) +
  geom_point(data = select(coord_size, -n_dip), aes(x= Longitude, y= Latitude, size = n_mul), pch = 21,col = "red" ) + #, col = "white", alpha = 0.5) +
  scale_size_continuous(name = "# of fish", breaks = c(1, 3, 7,10)) +
  theme(legend.position = c(0.06, 0.12))

ggdata_coord_size <- pivot_longer(coord_size, cols = c(n_dip, n_mul),names_to = "species", values_to = "sample_size")

sampling_both <- med_base +
  geom_point(data = ggdata_coord_size, aes(x= Longitude, y= Latitude, size = sample_size, col = species), pch = 21) +
  scale_size_continuous(name = "# of fish", breaks = c(1, 3, 7,10)) +
  scale_colour_manual(values= c("blue", "red"), labels = c("Diplodus sargus", "Mullus surmuletus")) +
  guides(colour=guide_legend(ncol=2),
         size  =guide_legend(ncol=1)) +
  theme(legend.position = c(0.2, 0.12),legend.direction = "horizontal")
sampling_both

ggsave(sampling_both, filename ="maps/sampling_map_both_size_col.pdf", width = 13, height = 7)  


#### map EEZs ####
# source : https://www.marineregions.org/downloads.php

library(sf)
eez.boundaries <- st_read("~/Documents/Data/GIS/MarineRegions/World_EEZ_v11_20191118/eez_boundaries_v11.shp")
class(eez.boundaries)
eez <- st_read("~/Documents/Data/GIS/MarineRegions/World_EEZ_v11_20191118/eez_v11.shp")
# crop to med extent because too large to plot
eez.med <- st_crop(eez, xmin=-13, xmax=42, ymin=25,ymax=50)

med_eez <- ggplot() + 
  geom_sf(data = eez.med, fill = NA) +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = "black") + 
  coord_sf(xlim=c(-8,37), ylim=c(29.5,46.5)) +
  labs(x= "Longitude (°)", y = "Latitude (°)") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))  # add black border on top again

ggsave(med_eez, filename= "maps/eez_shapes_med.pdf", width = 13, height = 7)

med_eez.b <- ggplot() + 
  geom_sf(data = eez.boundaries, color = "gray47") +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = "gray47") + 
  coord_sf(xlim=c(-8,37), ylim=c(29.5,46.5)) +
  labs(x= "Longitude (°)", y = "Latitude (°)") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))  # add black border on top again

ggsave(med_eez.b, filename= "maps/eez_boundaries_med.pdf", width = 13, height = 7)
 
# add sampling
sampling_eez <- med_eez.b +
  geom_point(data = map_data, aes(x=Longitude, y = Latitude, shape = species), cex = 3) +
  scale_shape_discrete(labels=c("both species", "Diplodus sargus", "Mullus surmuletus")) +
  theme(legend.position = c(0.07, 0.07),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank())
ggsave(sampling_eez, filename= "maps/sampling_both_shape_eez.pdf", width = 13, height = 7)

#### Marine Ecoregions of the World ####
# source: https://www.worldwildlife.org/publications/marine-ecoregions-of-the-world-a-bioregionalization-of-coastal-and-shelf-areas
library(sf)
library(forcats)
meow <- st_read("~/Documents/Data/GIS/MarineRegions/MEOW/meow_ecos.shp")
class(meow)
# crop to med extent because too large to plot
meow.med <- st_crop(meow, xmin=-13, xmax=42, ymin=25,ymax=50)
# remove non-med seas
meow.med.bis <- meow.med %>% filter(ECOREGION %in% c("Alboran Sea", "Western Mediterranean",
                                                     "Adriatic Sea", "Ionian Sea", "Tunisian Plateau/Gulf of Sidra",
                                                     "Aegean Sea", "Levantine Sea"))
meow.med.bis$ECOREGION <- fct_recode(meow.med.bis$ECOREGION,'Tunisian Plateau' = "Tunisian Plateau/Gulf of Sidra")

meow.med.bis <- cbind(meow.med.bis, st_coordinates(st_centroid(meow.med.bis))) # add centroid for lables
# coordinates almeria oran front
aof <- data.frame(city = c("Almeria", "Oran"),lon=c(-2.4597400,-0.6416700), lat = c(36.8381400,35.6911100))

med_meow <- ggplot() + 
  geom_sf(data = meow.med.bis, aes(fill = ECOREGION), alpha = 0.6) +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = "black") + 
  coord_sf(xlim=c(-8,37), ylim=c(29.5,46.5)) +
  scale_fill_brewer() +
  geom_text(data = meow.med.bis, aes(X, Y, label = ECOREGION), size = 3, fontface = "italic",
            angle = c(-40, 0, 0, 0, 0, 0, 0),
            nudge_x = c(0, 0, 0, 0, -1, -0.5, 0),
            nudge_y = c(-1,0,2,0,-2.1,0.9,-1)) +
  #geom_path(data = aof, aes(x=lon, y = lat), linetype = 2) +
  labs(x= "Longitude (°)", y = "Latitude (°)") +
  theme(legend.position = "none",
    panel.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))  # add black border on top again

ggsave(med_meow, filename= "maps/ecoregions_med.pdf", width = 13, height = 7)

# adjust colours
# ecoregion palette
colregion <- c("#FF0000","#5FB7FF","#1A8C18","#8D0000","#34638D","#FFA600","#99CF1C")
med_meow_cadj <- ggplot() + 
  geom_sf(data = meow.med.bis, aes(fill = ECOREGION), alpha = 0.8) +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = "black") + 
  coord_sf(xlim=c(-8,37), ylim=c(29.5,46.5)) +
  scale_fill_manual(values = colregion) +
  geom_text(data = meow.med.bis, aes(X, Y, label = ECOREGION), size = 3, fontface = "italic",
            angle = c(-40, 0, 0, 0, 0, 0, 0),
            nudge_x = c(0, 0, 0, 0, -1, -0.5, 0),
            nudge_y = c(-1,0,2,0,-2.1,0.9,-1)) +
  #geom_path(data = aof, aes(x=lon, y = lat), linetype = 2) +
  labs(x= "Longitude (°)", y = "Latitude (°)") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))  # add black border on top again

ggsave(med_meow_cadj, filename= "maps/ecoregions_med_coladj.pdf", width = 13, height = 7)

# add sampling
# add sampling
sampling_meow <- med_meow_cadj +
  #geom_point(data = map_data, aes(x=Longitude, y = Latitude, shape = species), cex = 3) +
  #scale_shape_discrete(labels=c("both species", "Diplodus sargus", "Mullus surmuletus")) +
  geom_text(data=map_data, aes(x=Longitude, y=Latitude, label=STRATA))
ggsave(sampling_meow, filename= "maps/sampling_both_shape_meow_coladj.pdf", width = 13, height = 7)
