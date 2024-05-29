library(tidyverse)

setwd("~/Optimization/")


##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

# Bin data by 30 m
BinData <- function(df){
  
  # Make 30m bin ref table
  BinWidth <- 30 # set bin width
  MaxSeq <- ceiling(max(df$distance)/BinWidth) # Number of bins
  bins <- data.frame(Bin = seq(1:MaxSeq),  # Make ref table
                     StartDist= BinWidth*(seq(1:MaxSeq) - 1),
                     EndDist = BinWidth*seq(1:MaxSeq)) #ref table
  
  # Assign a bin number based on ref table
  df$Bin <- NA
  for(i in 1:nrow(bins)){
    row = bins[i,]
    df <- df %>%
      dplyr::mutate(Bin = ifelse(distance >= row$StartDist[[1]] & distance < row$EndDist[[1]],
                                 row$Bin[[1]], Bin))
  }
  
  # Find average distance for each bin and year 
  vac_binavg <- df %>%
    group_by(Bin, year) %>%
    summarise(mean_distance = mean(distance))
  
  # Find the number of houses in each bin by vaccine status and year and merge
  # other bin variables
  df_out <- data.frame(plyr::count(df%>%select(vac_status, Bin, year)))%>%
    tidyr::pivot_wider(names_from = vac_status, values_from = freq)%>%
    rename(num_not_vac = '0',
           num_vac ='1')%>%
    mutate(num_vac = tidyr::replace_na(num_vac, 0),
           num_not_vac = tidyr::replace_na(num_not_vac, 0),
           num_of_houses = num_not_vac+num_vac,
           VaccFreq = num_vac/num_of_houses) %>%
    arrange(year, Bin) %>%
    left_join(., vac_binavg, by = c("year", "Bin")) %>%
    left_join(., bins)
  
  return(df_out)
}

##-----------------------------------------------------------------------------
# 1. Load data and fit Poisson regression model
##-----------------------------------------------------------------------------

# Bin raw data
load("optimization_data.rda")
vac_binned <- BinData(survey)

# Fit Poisson model
poisreg<-glm(num_vac~mean_distance+offset(log(num_of_houses)),
             family=poisson(link=log), data=vac_binned)
summary(poisreg)

# Calculate regression curve
coef(poisreg)
#   (Intercept) mean_distance 
# -0.2890246742 -0.0006483417 
summary(vac_binned$mean_distance)  # range from 13-1068 m
xvals <- seq(13, 1068, by = 1)
yvals_pois <- exp(-0.2890246742 - 0.0006483417*xvals)
xy_pois <- bind_rows(x = xvals, y = yvals_pois, type = "Poisson")


##-----------------------------------------------------------------------------
# 2. Plot survey data and Poisson regression curve
##-----------------------------------------------------------------------------

# Colors for plotting
cols <- c("#cab2d6", "#fcb353", "#b2df8a", "skyblue3")
pdf("_figures manuscript/R_output/Poissonregression.pdf", width = 5, height = 4)
ggplot()+
  theme_classic()+
  geom_point(data=vac_binned, aes(x=mean_distance, y=VaccFreq, size=num_of_houses,
                                  color=as.factor(year)))+
  scale_color_manual(values = cols) +
  geom_line(data = xy_pois, linewidth = 1.1,
            aes(x = x, y = y, linetype = type)) +
  xlab("Binned mean walking distance, m") + ylab("MDVC participation probability")
dev.off()

