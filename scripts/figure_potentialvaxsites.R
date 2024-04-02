library(tidyverse)
library(sf)
library(ggmap)
library(ggsn)
## REMEMBER TO DELETE THESE WHEN UPLOADING TO GITHUB***
register_google(key="AIzaSyBrwTLG8vz3RPeoVUrM4e6wcWMGOimtTrs")
register_stadiamaps(key="07107bab-7f12-480f-812b-6ed99f0cb3e2")

##-----------------------------------------------------------------------------
# 1. Load and clean data 
##-----------------------------------------------------------------------------

# Load ASA optimization data
setwd("~/Optimization")
load("data/asa_optimization_data.rda")

# Re-project houses and supply points (vax sites)
houses <- st_transform(houses, crs = 4326)
supply <- st_transform(supply, crs = 4326)
asa <- st_transform(asa, crs = 4326)

##-----------------------------------------------------------------------------
# 2. Create map of ASA with potential vaccination sites
##-----------------------------------------------------------------------------

# Create basemap for ASA 
bbox <- st_bbox(asa)
pad <- 0.005
pad2 <- 0.01
map_borders <- c(bottom = as.numeric(bbox$ymin) - pad, 
                 top = as.numeric(bbox$ymax) + pad, 
                 left = as.numeric(bbox$xmin) - pad2, 
                 right = as.numeric(bbox$xmax) + pad2)

# Get basemap -- Stamen terrain 
asa_terrain <- get_stadiamap(bbox = map_borders, zoom = 14, maptype = "stamen_terrain")
basemap <- ggmap(asa_terrain)

# Save different color versions

pdf("_figures manuscript/R_output/potentialsites.pdf", height = 6, width = 5.5)
basemap +
  geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", lwd = 0.9) + #
 # geom_sf(data = houses, inherit.aes = FALSE, color = "#bf7271", alpha = 0.1, size = 0.5) +
  geom_sf(data = houses, inherit.aes = FALSE, color = "#a88d4a", alpha = 0.2, size = 0.3) +
  geom_sf(data = supply, inherit.aes = FALSE, size = 2, color = "black", shape = 23,
          fill = "#C00000") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()

pdf("figures/R_output/figure1_potentialtents_lightred.pdf", height = 6, width = 5.5)
basemap +
  geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
          lwd = 1) +
  geom_sf(data = supply, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "#f7cdcd") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
              dist = 1, transform = TRUE, model = 'WGS84')
dev.off()

pdf("figures/R_output/figure1_potentialtents_darkred.pdf", height = 6, width = 5.5)
basemap +
  geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
          lwd = 1) +
  geom_sf(data = supply, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "#bf4141") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()


