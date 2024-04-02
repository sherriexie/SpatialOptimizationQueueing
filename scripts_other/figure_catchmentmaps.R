library(tidyverse)
library(sf)
library(ggmap)
library(ggsn)
library(viridis)
library(scales)
library(RColorBrewer)

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

GetArrivalsDist <- function(pm, houses){
  
  # Get site assignment and probability of arrival for reach house
  probability_byhouse <- apply(pm, MARGIN = 1, max)
  assignment_byhouse <- apply(pm, MARGIN = 1, which.max)
  names <- colnames(pm)
  assignment_names <- factor(assignment_byhouse, labels = names, ordered = T)
  
  # Add variables to houses sf with site assignments/probabilities of arrival
  houses$assignment <- assignment_names
  houses$assignmentprob <- probability_byhouse
  
  # Make summary data frame of expected number of arrivals at each site
  # Note: raw arrivals are scaled by habitability (0.57), dog ownership rate 
  # (0.4), and number of dogs per hh among hh's with dogs (1.86)
  arrivals <- st_drop_geometry(houses) %>%
    group_by(assignment) %>%
    summarise(arrivals_raw = sum(assignmentprob)) %>%
    mutate(num_arrivals = arrivals_raw*0.57*0.4*1.86) %>%
    arrange(desc(num_arrivals)) %>%
    dplyr::select(assignment, num_arrivals)
  arrivals$label <- as.factor(1:ncol(pm))
  
  return(arrivals)
}

MakeArrivalsMap <- function(supplyset, pm, houses, cols, title, lims, 
                            legend = F, grey = F, mix = F){
  
  # Find closest vax site for each house and save info using vax site ID's
  names <- colnames(pm)
  assignment_byhouse <- apply(pm, MARGIN = 1, which.max)
  assignment_names <- factor(assignment_byhouse, labels = names, ordered = T)
  #table(assignment_names)  # Check labels assigned as expected
  houses$assignment <- assignment_names
  
  # Get expected number of arrivals at each vax site using GetArrivalsDist fxn
  # and merge with supply points
  arrivals <- GetArrivalsDist(pm, houses) %>%
    rename(VaccPoint = assignment)
  supplyset <- left_join(supplyset, arrivals, by = "VaccPoint")
  
  # Add text label for the number of arrivals
  supplyset$num_label <- round(supplyset$num_arrivals)
  
  map <- ggplot() +
    geom_sf(data = houses, aes(col = assignment), inherit.aes = FALSE,
            size = 1, alpha = 0.3) +
    geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
            lwd = 0.9) +    scale_color_manual(values= cols) +
    geom_sf(data = supplyset, aes(size = num_arrivals), fill = "white", shape = 21,
            color = "black", inherit.aes = FALSE) +
    scale_size_area(limits = lims, max_size = 15) +
    geom_sf_text(data = supplyset, aes(label = num_label), inherit.aes = FALSE) +
    ggtitle(title) +
    scalebar(location = "bottomright", data = houses, dist_unit = "km",
             dist = 1, transform = TRUE, model = 'WGS84') +
    theme_void() +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  
  if(legend == FALSE){
    map <- map + 
      theme(legend.position = "none")
  }
  
  if(grey == TRUE){  # plot all grey circles
    map <- ggplot() +
      geom_sf(data = houses, aes(col = assignment), inherit.aes = FALSE,
              size = 1, alpha = 0.3) +
      geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
              lwd = 1) +    scale_color_manual(values= cols) +
      geom_sf(data = supplyset, aes(size = num_arrivals), fill = "grey", shape = 21,
              color = "black", inherit.aes = FALSE) +
      scale_size_area(limits = lims, max_size = 15) +
      geom_sf_text(data = supplyset, aes(label = num_label), inherit.aes = FALSE) +
      ggtitle(title) +
      scalebar(location = "bottomright", data = houses, dist_unit = "km",
               dist = 1, transform = TRUE, model = 'WGS84') +
      theme_void() +
      theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  }
  
  if(mix == TRUE){  # plot mix of grey/white circles depending on whether site 
                    # appears in the original solution (grey) or not (white)
    map <- ggplot() +
      geom_sf(data = houses, aes(col = assignment), inherit.aes = FALSE,
              size = 1, alpha = 0.3) +
      geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
              lwd = 1) +    
      scale_color_manual(values= cols) +
      geom_sf(data = supplyset, aes(size = num_arrivals, fill = status), shape = 21,
              color = "black", inherit.aes = FALSE) +
      scale_fill_manual(values = c("white", "grey")) +
      scale_size_area(limits = lims, max_size = 15) +
      geom_sf_text(data = supplyset, aes(label = num_label), inherit.aes = FALSE) +
      ggtitle(title) +
      scalebar(location = "bottomright", data = houses, dist_unit = "km",
               dist = 1, transform = TRUE, model = 'WGS84') +
      theme_void() +
      theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))  
    
  }
  
  return(map)
  
}

##-----------------------------------------------------------------------------
# 1. Load data
##-----------------------------------------------------------------------------

# Load ASA optimization data
setwd("~/Optimization")
load("data/asa_optimization_data.rda")

# Re-project houses and supply points (vax sites)
houses <- st_transform(houses, crs = 4326)
supply <- st_transform(supply, crs = 4326)

# Indices of vax sites optimized by the queue-naive, the queue-conscious/low
# attrition, and queue-conscious/high attrition algorithms
v_naive  <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 41, 45, 47, 
              57, 61, 63)
v_low <- c(2, 3, 5, 9, 15, 19, 20, 21, 23, 26, 28, 29, 30, 31, 36, 37, 41, 49, 
           54, 70)
v_high <- c(1, 2, 5, 6, 9, 12, 15, 19, 21, 26, 28, 30, 37, 42, 44, 45, 57, 66, 
            68, 70)

# Use vax site indices to select optimized vax sites
s_naive<- supply[v_naive,]
s_low <- supply[v_low,]
s_high <- supply[v_high,]

# Subset probability matrix for each group
pm_naive<- probmatrix[,v_naive]
pm_low <- probmatrix[,v_low]
pm_high <- probmatrix[,v_high]

# Names of selected sites
s_naive$VaccPoint
s_low$VaccPoint
s_high$VaccPoint

# Create color palette
display.brewer.pal(12, "Paired")
colpal <- brewer.pal(12, "Paired")
# - tweak colors
colpal[c(1:2, 5:6, 8)] <- c("#cdeafa", "skyblue3", "#f09292", "#C00000", 
                            "#e37910")

##-----------------------------------------------------------------------------
# 2. Make catchment map for queue-naive solution
##-----------------------------------------------------------------------------

# Create color palette
col_index_naive <- c(1,2,6,7,9:12, 3:6, 8, 10:11, 2:3, 11, 7:8)
cols <- colpal[col_index_naive]

# Make catchment map
pdf("_figures manuscript/R_output/catchmentmap_queuenaive.pdf", height = 5.5, 
    width = 5.5)
MakeArrivalsMap(s_naive, pm_naive, houses, cols, "Queue-naive", 
                lims = c(107.9, 501.5), legend = T)
dev.off()

##-----------------------------------------------------------------------------
# 3. Make catchment map for low attrition solution
##-----------------------------------------------------------------------------

# Create color palette
col_index_low <- c(2:4, 7, 9, 11:12, 1, 3:11, 4:5, 6)
cols <- colpal[col_index_low]

# Make catchment map
pdf("_figures manuscript/R_output/catchmentmap_lowattrition.pdf", height = 5.5, 
    width = 5.5)
MakeArrivalsMap(s_low, pm_low, houses, cols, "Low attrition", 
                lims = c(107.9, 501.5), legend = T)
dev.off()

##-----------------------------------------------------------------------------
# 4. Make catchment map for high attrition solution
##-----------------------------------------------------------------------------

# Create color palette
col_index_high <- c(1:2, 4:5, 7:9, 11, 1, 4:5, 7, 10, 12, 6, 2, 11, 9:10, 6)
cols <- colpal[col_index_high]

# Make catchment map
pdf("_figures manuscript/R_output/catchmentmap_highattrition.pdf", height = 5.5, 
    width = 5.5)
MakeArrivalsMap(s_high, pm_high, houses, cols, "High attrition", 
                lims = c(107.9, 501.5), legend = T)
dev.off()

