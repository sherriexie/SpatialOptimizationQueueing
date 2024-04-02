library(tidyverse)
library(splines)
library(scales)
library(gridExtra)

setwd("~/Optimization")

#------------------------------------------------------------------------------
# 0. User-defined functions
#------------------------------------------------------------------------------

SummarizeData <- function(df){
  
  reps <- length(unique(df$rep))
  
  # Summarize number vaxed and number lost to attrition for each 4-day campaign
  # at an MDVC site
  simulation_df <- df %>%
    group_by(n, rep) %>%
    summarise(num_vaxed = sum(num_vaxed),
              num_attrition = sum(num_attrition)) %>%
    mutate(num_arrivals = n*4) %>%
    select(rep, num_arrivals, num_vaxed, num_attrition)
  
  # Summarize number vaxed (median and IQR) for each number of arrivals n
  vax_summary <- simulation_df %>%
    group_by(num_arrivals) %>%
    summarise(median = median(num_vaxed),
              mean = mean(num_vaxed),
              lower_25 = sort(num_vaxed)[.25*reps],
              upper_75 = sort(num_vaxed)[.75*reps])
  
  # Summarize number lost to attrition (median and IQR) for each number of 
  # arrivals n
  attrition_summary <- simulation_df %>%
    group_by(num_arrivals) %>%
    summarise(median = median(num_attrition),
              mean = mean(num_attrition),
              lower_25 = sort(num_attrition)[.25*reps],
              upper_75 = sort(num_attrition)[.75*reps])
  
  return(list(simulation_df = simulation_df, vax_df = vax_summary, 
              attrition_df = attrition_summary))
  
}

nrep_df2 <- nrep_df %>%
  group_by(n, rep) %>%
  summarise(num_vaxed = sum(num_vaxed),
            num_attrition = sum(num_attrition)) %>%
  mutate(num_arrivals = n*4) %>%
  select(rep, num_arrivals, num_vaxed, num_attrition)

PlotVax <- function(df, title = NULL, margin = 0){
  ggplot(df, aes(x = num_arrivals)) +
    geom_line(aes(y = median), linetype = "dashed", color = "skyblue3", lwd = 1.1) +
    geom_ribbon(aes(ymin = lower_25, ymax = upper_75), fill = "skyblue3", alpha = 0.3) +
    xlab("Number of dogs arrived") +
    ylab("Number vaccinated") +
    ggtitle(title) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_bw() +
    theme(axis.text=element_text(size=22),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24, face='bold', hjust = 0.5),
          plot.margin = unit(rep(margin, 4), "cm")) 
}

PlotAttrition <- function(df, title = NULL, margin = 0){
  ggplot(df, aes(x = num_arrivals)) +
    geom_line(aes(y = median), linetype = "dashed", color = "#C00000", lwd = 1.1) +
    geom_ribbon(aes(ymin = lower_25, ymax = upper_75), fill = "#C00000", alpha = 0.3) +
    xlab("Number of dogs arrived") +
    ylab("Number lost to attrition") +
    ggtitle(title) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_bw() +
    theme(axis.text=element_text(size=22),
          axis.title=element_text(size=24),
          plot.title=element_text(size=24, face='bold', hjust = 0.5),
          plot.margin = unit(rep(margin, 4), "cm")) 
}

#------------------------------------------------------------------------------
# 1. Load simulation outputs
#------------------------------------------------------------------------------

a_low <- read_rds("_output manuscript/simulations_arrdensityA_lowattrition.rds")
a_high <- read_rds("_output manuscript/simulations_arrdensityA_highattrition.rds")
b_low <- read_rds("_output manuscript/simulations_arrdensityB_lowattrition.rds")
b_high <- read_rds("_output manuscript/simulations_arrdensityB_highattrition.rds")
b_low <- read_rds("_output manuscript/simulations_arrdensityB_lowattrition.rds")
b_high <- read_rds("_output manuscript/simulations_arrdensityB_highattrition.rds")
c_low <- read_rds("_output manuscript/simulations_arrdensityC_lowattrition.rds")
c_high <- read_rds("_output manuscript/simulations_arrdensityC_highattrition.rds")
d_low <- read_rds("_output manuscript/simulations_arrdensityD_lowattrition.rds")
d_high <- read_rds("_output manuscript/simulations_arrdensityD_highattrition.rds")

#------------------------------------------------------------------------------
# 2. Plot simulation results
#------------------------------------------------------------------------------

# Plot simulation results for densities + low-attrition parameters
v_a <- PlotVax(SummarizeData(a_low)$vax_df, "Density A, low attrition", 1)
v_b <- PlotVax(SummarizeData(b_low)$vax_df, "Density B, low attrition", 1)
v_c <- PlotVax(SummarizeData(c_low)$vax_df, "Density C, low attrition", 1)
v_d <- PlotVax(SummarizeData(d_low)$vax_df, "Density D, low attrition", 1)
a_a <- PlotAttrition(SummarizeData(a_low)$attrition_df, "Density A, low attrition", 1)
a_b <- PlotAttrition(SummarizeData(b_low)$attrition_df, "Density B, low attrition", 1)
a_c <- PlotAttrition(SummarizeData(c_low)$attrition_df, "Density C, low attrition", 1)
a_d <- PlotAttrition(SummarizeData(d_low)$attrition_df, "Density D, low attrition", 1)

# Plot simulation results for densities + high-attrition parameters
v_ah <- PlotVax(SummarizeData(a_low)$vax_df, "Density A, high attrition", 1)
v_bh <- PlotVax(SummarizeData(b_low)$vax_df, "Density B, high attrition", 1)
v_ch <- PlotVax(SummarizeData(c_low)$vax_df, "Density C, high attrition", 1)
v_dh <- PlotVax(SummarizeData(d_low)$vax_df, "Density D, high attrition", 1)
a_ah <- PlotAttrition(SummarizeData(a_low)$attrition_df, "Density A, high attrition", 1)
a_bh <- PlotAttrition(SummarizeData(b_low)$attrition_df, "Density B, high attrition", 1)
a_ch <- PlotAttrition(SummarizeData(c_low)$attrition_df, "Density C, high attrition", 1)
a_dh <- PlotAttrition(SummarizeData(d_low)$attrition_df, "Density D, high attrition", 1)

pdf("_figures manuscript/R_output/sa_simulationresults.pdf", width = 24, height = 24)
grid.arrange(v_a, a_a, v_ah, a_ah,
             v_b, a_b, v_bh, a_bh,
             v_c, a_c, v_ch, a_ch,
             v_d, a_d, v_dh, a_dh, ncol = 4)
dev.off()

#------------------------------------------------------------------------------
# 3a. Distribution A - low
#------------------------------------------------------------------------------

output <- SummarizeData(a_low)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated", 
     main = "Density A, low attrition")
points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density A, low attrition")
points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), 
          "_output manuscript/splines_a_low.rds")

#------------------------------------------------------------------------------
# 3b. Distribution A - high
#------------------------------------------------------------------------------

output <- SummarizeData(a_high)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated",
     main = "Density A, high attrition")
points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density A, high attrition")
points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), 
          "_output manuscript/splines_a_high.rds")

#------------------------------------------------------------------------------
# 4a. Distribution B - low
#------------------------------------------------------------------------------

output <- SummarizeData(b_low)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated",
     main = "Density B, low attrition")
points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density B, low attrition")
points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), 
          "_output manuscript/splines_b_low.rds")

#------------------------------------------------------------------------------
# 4b. Distribution B - high
#------------------------------------------------------------------------------

output <- SummarizeData(b_high)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated",
     main = "Density B, high attrition")

points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density B, high attrition")

points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), "_output manuscript/splines_b_high.rds")

#------------------------------------------------------------------------------
# 5a. Distribution C - low
#------------------------------------------------------------------------------

output <- SummarizeData(c_low)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated",
     main = "Density C, low attrition")
points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density C, low attrition")
points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), 
          "_output manuscript/splines_c_low.rds")

#------------------------------------------------------------------------------
# 5b. Distribution C - high
#------------------------------------------------------------------------------

output <- SummarizeData(c_high)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated",
     main = "Density C, high attrition")
points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density C, high attrition")

points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), 
          "_output manuscript/splines_c_high.rds")

#------------------------------------------------------------------------------
# 6a. Distribution D - low
#------------------------------------------------------------------------------

output <- SummarizeData(d_low)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated",
     main = "Density D, low attrition")

points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density D, low attrition")

points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), 
          "_output manuscript/splines_d_low.rds")

#------------------------------------------------------------------------------
# 6b. Distribution D - high
#------------------------------------------------------------------------------

output <- SummarizeData(d_high)
PlotVax(output$vax_df)
PlotAttrition(output$attrition_df)

# Fit spline model to data
spline_vax <- glm(median ~ns(num_arrivals, df = 7), data = output$vax_df)
spline_att <- glm(median ~ns(num_arrivals, df = 7), data = output$attrition_df)

# Check fit for vax spline
plot(output$vax_df$num_arrivals, output$vax_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number vaccinated",
     main = "Density D, high attrition")
points(output$vax_df$num_arrivals, 
       predict(spline_vax, data.frame(num_arrivals = output$vax_df$num_arrivals)), 
       type = "l", add = T, col = "skyblue3", lwd = 3)

# Check fit for attrition spline
plot(output$attrition_df$num_arrivals, output$attrition_df$median, pch = 16, 
     xlab = "Number of dogs arrived", ylab = "Number lost to attrition",
     main = "Density D, high attrition")
points(output$attrition_df$num_arrivals, 
       predict(spline_att, data.frame(num_arrivals = output$attrition_df$num_arrivals)), 
       type = "l", add = T, col = "#C00000", lwd = 3)

# Save splines
write_rds(list(spline_vax = spline_vax, spline_att = spline_att), 
          "_output manuscript/splines_d_high.rds")


