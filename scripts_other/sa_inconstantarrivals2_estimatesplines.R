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

PlotVax <- function(df, spline, title = NULL, margin = 0){
  
  df$pred <- predict(spline, data.frame(num_arrivals = df$num_arrivals))
  
  ggplot(df, aes(x = num_arrivals)) +
    geom_line(aes(y = median), color = "skyblue3", lwd = 1.1) +
    geom_ribbon(aes(ymin = lower_25, ymax = upper_75), fill = "skyblue3", alpha = 0.3) +
    geom_line(aes(y = pred), linetype = "dashed", color = "black", lwd = 1) +
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

PlotAttrition <- function(df, spline, title = NULL, margin = 0){
  
  df$pred <- predict(spline, data.frame(num_arrivals = df$num_arrivals))
  
  ggplot(df, aes(x = num_arrivals)) +
    geom_line(aes(y = median), linetype = "dashed", color = "#C00000", lwd = 1.1) +
    geom_ribbon(aes(ymin = lower_25, ymax = upper_75), fill = "#C00000", alpha = 0.3) +
    geom_line(aes(y = pred), linetype = "dashed", color = "black", lwd = 1) +
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
# 2. Create spline functions
#------------------------------------------------------------------------------

# Create splines for densities a-d + low-attrition parameters
spline_a_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(a_low)$vax_df)
spline_a_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(a_low)$attrition_df)
spline_b_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(b_low)$vax_df)
spline_b_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(b_low)$attrition_df)
spline_c_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(c_low)$vax_df)
spline_c_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(c_low)$attrition_df)
spline_d_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(d_low)$vax_df)
spline_d_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(d_low)$attrition_df)

# Create splines for densities a-d + high-attrition parameters
spline_ah_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(a_high)$vax_df)
spline_ah_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(a_high)$attrition_df)
spline_bh_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(b_high)$vax_df)
spline_bh_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(b_high)$attrition_df)
spline_ch_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(c_high)$vax_df)
spline_ch_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(c_high)$attrition_df)
spline_dh_vax <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(d_high)$vax_df)
spline_dh_att <- glm(median ~ns(num_arrivals, df = 7), data = SummarizeData(d_high)$attrition_df)

#------------------------------------------------------------------------------
# 3. Plot results
#------------------------------------------------------------------------------

# Plot simulation results and splines for densities a-d + low-attrition parameters
v_a <- PlotVax(SummarizeData(a_low)$vax_df, spline_a_vax, "Density A, low attrition", 1)
v_b <- PlotVax(SummarizeData(b_low)$vax_df, spline_b_vax, "Density B, low attrition", 1)
v_c <- PlotVax(SummarizeData(c_low)$vax_df, spline_c_vax, "Density C, low attrition", 1)
v_d <- PlotVax(SummarizeData(d_low)$vax_df, spline_d_vax, "Density D, low attrition", 1)
a_a <- PlotAttrition(SummarizeData(a_low)$attrition_df, spline_a_att, "Density A, low attrition", 1)
a_b <- PlotAttrition(SummarizeData(b_low)$attrition_df, spline_b_att, "Density B, low attrition", 1)
a_c <- PlotAttrition(SummarizeData(c_low)$attrition_df, spline_c_att, "Density C, low attrition", 1)
a_d <- PlotAttrition(SummarizeData(d_low)$attrition_df, spline_d_att, "Density D, low attrition", 1)

# Plot simulation results and splines for densities a-d + high-attrition parameters
v_ah <- PlotVax(SummarizeData(a_high)$vax_df, spline_ah_vax, "Density A, high attrition", 1)
v_bh <- PlotVax(SummarizeData(b_high)$vax_df, spline_bh_vax, "Density B, high attrition", 1)
v_ch <- PlotVax(SummarizeData(c_high)$vax_df, spline_ch_vax, "Density C, high attrition", 1)
v_dh <- PlotVax(SummarizeData(d_high)$vax_df, spline_dh_vax, "Density D, high attrition", 1)
a_ah <- PlotAttrition(SummarizeData(a_high)$attrition_df, spline_ah_att, "Density A, high attrition", 1)
a_bh <- PlotAttrition(SummarizeData(b_high)$attrition_df, spline_bh_att, "Density B, high attrition", 1)
a_ch <- PlotAttrition(SummarizeData(c_high)$attrition_df, spline_ch_att, "Density C, high attrition", 1)
a_dh <- PlotAttrition(SummarizeData(d_high)$attrition_df, spline_dh_att, "Density D, high attrition", 1)

pdf("_figures manuscript/R_output/sa_simulationresults.pdf", width = 24, height = 24)
grid.arrange(v_a, a_a, v_ah, a_ah,
             v_b, a_b, v_bh, a_bh,
             v_c, a_c, v_ch, a_ch,
             v_d, a_d, v_dh, a_dh, ncol = 4)
dev.off()