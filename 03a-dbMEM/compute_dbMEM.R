# ## Code from Laura Benestan: https://github.com/laurabenestan/Moran-Eigenvector-Maps-MEMs

# Moran's eigenvectors maps (MEM) are complete and easy-to-use mathematical objects 
# that aims to partition the spatial distribution across samples.

# load libraries
library(codep)
library(adespatial)
library(adegraphics)
library(vegan)
library(car)
library(dplyr)
library(data.table)
library(ggplot2)
library(sf)
library(tidyr)
library(patchwork)

# load data
# spatial coordinates with ind labels, latitude and longitude information
   # diplodus
geo.dip <- read.table("../00-Data/20-envdata-EB/dip_individual_envData.csv", sep = ",", head = T, row.names = 1)
geo.dip <- geo.dip[,c("Individual", "Longitude", "Latitude")]
   # mullus
geo.mul <- read.table("../00-Data/20-envdata-EB/mul_individual_envData.csv", sep = ",", head = T, row.names = 1)
geo.mul <- geo.mul[,c("Individual", "Longitude", "Latitude")]

geo <- geo.dip

# keep latitude and longitude in this order as for the function gcd.hf, latitude needs to be first

coor <- geo[, 3:2]
coor.xy <- geo[,2:3]

# look at the spatial distribution
plot(coor.xy, asp=1)

# compute spatial distances among sites accounting for the earth curvature
DistSpatial=gcd.hf(coor.xy)

# alternative: Dowload geographic distance 
# sea distance between individuals
# calculated in 21_GEODIST
dip_seadist <- as.dist(read.csv("../00-Data/21-GEODIST/output/seadist_diplodus_297ind.csv", row.names = 1))
mul_seadist <- as.dist(read.csv("../00-Data/21-GEODIST/output/seadist_mullus_467ind.csv", row.names = 1))

#Compute MEM by keeping default setting for truncation (length of the longest edge 
# of the minimum spanning tree will be used as the threshold) and just positive MEM.
# use sea distance

dbmem.dip <- dbmem(dip_seadist)
dbmem.mul <- dbmem(mul_seadist)

# add individual ID
rownames(dbmem.dip) <- labels(dip_seadist)
rownames(dbmem.mul) <- labels(mul_seadist)

summary(dbmem.dip)
summary(dbmem.mul)

# visualise MEM neighbourhood graph
s.label(geo.dip[,2:3], nb = attr(dbmem.dip, "listw"))
s.label(geo.mul[,2:3], nb = attr(dbmem.mul, "listw"))

# visualise MEM values
s.value(geo.dip[,2:3], dbmem.dip[,1]) # first dbMEM, represents large spatial scale
s.value(geo.mul[,2:3], dbmem.mul[,1])

# add geographic information
dbmem.geo.dip <- cbind(dbmem.dip, geo.dip)
dbmem.geo.mul <- cbind(dbmem.mul, geo.mul)

# export
write.csv(dbmem.geo.dip, file = "output/dbmem_seaDist_dip.csv")
write.csv(dbmem.geo.mul, file = "output/dbmem_seaDist_mul.csv")


# plot most informative MEMs from redundancy analysis ----
# with ggplot

# load dbMEMs
dbmem.geo.dip <- read.csv(file = "output/dbmem_seaDist_dip.csv", row.names = 1)
dbmem.geo.mul <- read.csv(file = "output/dbmem_seaDist_mul.csv", row.names = 1)

# change from wide to long format
dbmem_long.dip <- gather(dbmem.geo.dip, MEM, Value, MEM1:(ncol(dbmem.geo.dip)-3))
dbmem_long.mul <- gather(dbmem.geo.mul, MEM, Value, MEM1:(ncol(dbmem.geo.mul)-3))

# Calculate an MEM average value for each GPS point.
library(stringr)
dbmem_gps.dip <- dbmem_long.dip %>% 
  mutate(cell = gsub("i.*","",  .$Individual)) %>% 
  group_by(cell, Latitude, Longitude, MEM)%>%
  dplyr::summarise(mem_mean = mean(Value))
dbmem_wide.dip <- spread(dbmem_gps.dip, MEM, mem_mean)

dbmem_gps.mul <- dbmem_long.mul %>%
  mutate(cell = gsub("i.*","",  .$Individual)) %>% 
  group_by(cell, Latitude, Longitude, MEM)%>%
  dplyr::summarise(mem_mean = mean(Value))
dbmem_wide.mul <- spread(dbmem_gps.mul, MEM, mem_mean)

#  Download a high resolution map with the sf package
library(maps)
library(mapdata)

wH <- map_data("worldHires",  xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea

# Create a ggmap showing the most significant MEMs vectors (from RDA OrdiStep)
mycol <- c("#FF0000", "#FF8534", "#FFB63E", "#FFE95A", "#E8F171","#C1DA90", "#80CBA8","#4FB0BB","#00A1D0")

# function to map s values 
map_svalue <- function(df, scale, species) {
  library(dplyr)
  filter(df,MEM ==scale) %>% 
    ggplot() +
    geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
    coord_fixed(xlim = c(-8,37), ylim=c(29.5,47), ratio=1.2)+
    geom_point(aes(x=Longitude, y=Latitude, fill=mem_mean),shape = 21, size=4) +
    #facet_wrap(~MEM, nrow = 2, ncol = 1)+
    theme_bw() +
    labs(y="Latitude", x="Longitude", title = paste0(species," | ", scale)) +
    theme(#legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x=element_text(colour="black"),
      axis.text.y=element_text(colour="black")) +
    scale_fill_gradientn(colours = mycol, name = "S value") }

## Diplodus
map_svalue(dbmem_gps.dip, scale = "MEM1", species = "Diplodus sargus")

# Select only MEM6 and MEM17 for dip
dip.mem6  <- map_svalue(dbmem_gps.dip, scale = "MEM6",  species = "Diplodus sargus")
dip.mem16 <- map_svalue(dbmem_gps.dip, scale = "MEM16", species = "Diplodus sargus")
dip.mem7  <- map_svalue(dbmem_gps.dip, scale = "MEM7",  species = "Diplodus sargus")
dip.mem9  <- map_svalue(dbmem_gps.dip, scale = "MEM9",  species = "Diplodus sargus")

graph.dip <- dip.mem6 + dip.mem16 + dip.mem7 + dip.mem9
#ggsave(graph.dip, file = "figures/meanMEM_mostInformative_dip_adaptive.pdf")

## Mullus
# Select only MEM12 and MEM5 for mul.
mul.mem5  <- map_svalue(dbmem_gps.mul, scale = "MEM5",  species = "Mullus surmuletus")
mul.mem1  <- map_svalue(dbmem_gps.mul, scale = "MEM1",  species = "Mullus surmuletus")
mul.mem12 <- map_svalue(dbmem_gps.mul, scale = "MEM12", species = "Mullus surmuletus")
mul.mem17 <- map_svalue(dbmem_gps.mul, scale = "MEM17", species = "Mullus surmuletus")

graph.mul <- mul.mem5 + mul.mem1 + mul.mem12 + mul.mem17
#ggsave(graph.mul, file = "figures/meanMEM_mostInformative_mul_adaptive.pdf")

#all together
dip.mul.mem <- (dip.mem6 + dip.mem9) / (mul.mem5 + mul.mem12)
ggsave(dip.mul.mem, file = "figures/meanMEM_bothSP_mostInformative.pdf")
