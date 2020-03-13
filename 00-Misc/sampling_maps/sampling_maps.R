# create maps of sampling cell with size indicating sample size

# libraries
library(maps)
library(mapdata)
library(dplyr)
library(ggplot2)
library(scales)

library(stringr)
library(reshape)

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
points(coord_size$Longitude, coord_size$Latitude, pch=19, col="#053061", cex=rescale(coord_size$n_dip,to = c(1,4))) 
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


# map in ggplot2
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

