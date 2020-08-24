# create maps of sampling cell with size indicating sample size

# libraries
library(maps)
library(mapdata)
library(dplyr)
library(ggplot2)
library(scales)

library(stringr)
library(reshape)

# import stamps 
data_pic <- list.files(path = "../../00-Misc/stamps/", pattern="*.png",recursive = FALSE)
source_pic <- paste0("../../00-Misc/stamps/", data_pic)
pic <- lapply(source_pic,readPNG)
names(pic) <- str_sub(data_pic,1, -5)

# coord data
coord <- read.table("data/coord_seaconnect_tous.txt", sep = "\t", head = TRUE) 
pop_dip <- read.table("data/dip_297ind_population_map.txt", sep = "\t", head = TRUE)
pop_mul <- read.table("data/mul_424ind_population_map.txt", sep = "\t", head = TRUE)
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


# map in ggplot2 ----
# tutorial: http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html

wH <- map_data("worldHires",  xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea
# further subset dataset so don't plot whole polygons
# wH_sub <- wH[wH$long<37 & wH$long > c(-8) & wH$lat < 47 & wH$lat > 29.5,] # creates weird margins around map. rather set limits in ggplot

med_base <- ggplot() + 
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) + 
  coord_fixed( xlim=c(-8,37), ylim=c(29.5,46.5), ratio = 1.3) +
  #theme_nothing() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())

