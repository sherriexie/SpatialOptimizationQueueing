library(tidyverse)
library(scales)

setwd("~/Optimization")

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

Expr1 <- function(lambda, alpha, beta, mu, n){
  # This expression is used in calculating p0 for α>0/β>0
  num <- lambda^n * exp(-alpha*n*(n-1)/(2*mu)) * gamma(mu/beta)
  den <- beta^n * gamma(n+mu/beta)
  return(num/den)
}

Expr2 <- function(lambda, alpha, beta, mu, p0, n){
  # This expression is used in calculating first term of V (effective arrival
  # rate or lambda-bar) for α>0/β>0
  # -Specify p_n
  num <- lambda^n * exp(-alpha*n*(n-1)/(2*mu)) * gamma(mu/beta)
  den <- beta^n * gamma(n+mu/beta)
  p_n <- (num/den)*p0
  # -Specify lambda_n 
  lambda_n <- exp(-alpha*n/mu)*lambda
  return(p_n * lambda_n)
}

Expr3 <- function(lambda, alpha, beta, mu, p0, n){
  # This expression is used in calculating the second term of V (average rate
  # of reneging) for α>0/β>0
  # -Specify p_n
  num <- lambda^n * exp(-alpha*n*(n-1)/(2*mu)) * gamma(mu/beta)
  den <- beta^n * gamma(n+mu/beta)
  p_n <- (num/den)*p0
  return(p_n * (n-1) * beta)
}

Calculate_p0 <- function(lambda, alpha, beta, mu){
  # Pass n = 1, 2, .., 99 to Expr 1
  vec <- 1:99
  Expr1_out <- sapply(vec, function(x) Expr1(lambda, alpha, beta, mu, x))
  # Sum over all summands and find p0 
  p0 <- 1/(1+sum(Expr1_out, na.rm = T))
  return(p0)
}

GetArrivalsDF <- function(supplyset, probmatrix, alpha, beta, mu, t){
  
  # Get probability matrix with chosen supply points
  pm_set <- probmatrix[, supplyset]
  
  # Calculate number of houses that are expected to arrive
  site_index <- apply(pm_set, MARGIN = 1, which.max)
  arrive_prob <- apply(pm_set, MARGIN = 1, max)
  arrivals_bysite <- numeric()
  for(i in 1:ncol(pm_set)){
    arrivals_bysite[i] <- sum(arrive_prob[site_index == i])
  }
  
  # Scale number by habitability (0.57), dog ownership rate (0.4), and number of 
  # dogs per households (1.86) in ASA
  num_arrivals <- arrivals_bysite*0.57*0.4*1.86
  
  # Calculate lambdas (number of arrivals per minute at each site)
  lambda <- num_arrivals/t
  
  # Calculate p0 which we'll need to calculate p_n (which is in both terms of V)
  # -- note there is a different p0 for each value of lambda
  p0 <- numeric()
  for(i in 1:length(lambda)){
    p0[i] <- Calculate_p0(lambda[i], alpha, beta, mu)
  }
  # Calculate effective arrivals by passing n = 0, 1, 2, ..., 99 to Expr 2
  # Then calculate balking rate by subtracting eff arrival rate from lambda
  vec1 <- 0:99
  eff_arrival_rate <- balking_rate <- numeric()
  for (i in 1:length(lambda)){
    Expr2_out <- sapply(vec1, function(x) Expr2(lambda[i], alpha, beta, mu, p0[i], x))
    eff_arrival_rate[i] <- sum(Expr2_out, na.rm = T)
    balking_rate[i] <- lambda[i] - eff_arrival_rate[i]  
  }
  
  # Calculate average reneging rate by passing n = 1, 2, ..., 99 to Expr 3
  vec2 <- 1:99
  reneging_rate <- numeric()
  for(i in 1:length(lambda)){
    Expr3_out <- sapply(vec2, function(x) Expr3(lambda[i], alpha, beta, mu, p0[i], x))
    reneging_rate[i] <- sum(Expr3_out, na.rm = T)
  }
  
  # Calculate vaccination rate (V)
  V <- eff_arrival_rate - reneging_rate
  
  # Calculate number vaccinated, number balked, and number reneged
  num_vaxed <- V*t
  num_balked <- balking_rate*t
  num_reneged <- reneging_rate*t
  
  # Format output as a data frame
  wide_df <- data.frame(num_arrivals = num_arrivals, Balked = num_balked, 
                        Reneged = num_reneged, Vaccinated = num_vaxed) %>%
    arrange(desc(num_arrivals)) %>%
    dplyr::select(Balked, Reneged, Vaccinated)
  wide_df$site <- as.factor(1:length(num_arrivals))
  
  # Convert from wide to long format
  out_df <- gather(wide_df, outcome, number, Balked:Vaccinated, factor_key = T)
  
  # Print number and % of arrivals lost to attrition
  tot_arrivals <- round(sum(out_df$number))
  tot_attrition <- round(sum(out_df$number[out_df$outcome %in% c("Balked", 
                                                                 "Reneged")]))
  pct_attrition <- round(tot_attrition/tot_arrivals*100, 1)
  print(noquote(paste0(tot_attrition, " (", pct_attrition, "%) of ", tot_arrivals,
                       " arrivals lost to attrition")))
  
  # Print total vaccination coverage
  tot_dogs <- nrow(probmatrix)*0.57*0.4*1.86
  vax_dogs <- sum(out_df$number[out_df$outcome == "Vaccinated"])
  vax_coverage <- round(vax_dogs/tot_dogs*100,1)
  print(noquote(paste0("Vaccination coverage = ", vax_coverage, "%")))
  
  return(out_df)
}

PlotArrivals <- function(df, ymax, grey = F){
  plot <- ggplot(df) +
    geom_bar(aes(x = site, y = number, fill = outcome), stat = "identity", 
             width = 1, color = "black") +
    #scale_fill_manual(values = c("#C00000", "#db5151", "skyblue3", "#a6cee3", "#1f78b4"), 
    #                  name = "Outcome") +
    scale_fill_manual(values = c("#C00000", "#db5151", "skyblue3"), 
                      name = "Outcome") +
    xlab("Vaccination site") +
    ylab("Number of arrivals") +
    scale_y_continuous(breaks = pretty_breaks(), limits = c(0, ymax)) +
]    theme_minimal() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22, face = "bold"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22, face = "bold"),
          legend.spacing = unit(1, "cm"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    guides(fill = guide_legend(byrow = TRUE))
  
  if (grey == T){
    plot <- ggplot(df) +
      geom_bar(aes(x = site, y = number, fill = outcome), stat = "identity", 
               width = 1, color = "black") +
      scale_fill_manual(values = c("grey", "grey", "grey"), 
                        name = "Outcome") +
      xlab("Vaccination site") +
      ylab("Number of arrivals") +
      scale_y_continuous(breaks = pretty_breaks(), limits = c(0, ymax)) +
      theme_minimal() +
      theme(axis.text=element_text(size=20),
            axis.title=element_text(size=22, face = "bold"),
            legend.text=element_text(size=20),
            legend.title=element_text(size=22, face = "bold"),
            legend.spacing = unit(1, "cm"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) +
      guides(fill = guide_legend(byrow = TRUE))
  }
  
  return(plot)
}

##-----------------------------------------------------------------------------
# 1. Load data and get arrivals data frames
##-----------------------------------------------------------------------------

setwd("~/Optimization/")

# Indices of vax sites for queue-naive solution and solutions obtained assuming 
# low or high queue attrition rates
naive <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 41, 45, 47, 57, 
           61, 63)
low <- c(2, 3, 5, 9, 15, 19, 20, 21, 23, 26, 28, 29, 30, 31, 36, 37, 41, 49, 
         54, 70)
high <- c(1, 2, 5, 6, 9, 12, 15, 19, 21, 26, 28, 30, 37, 42, 44, 45, 57, 66, 
          68, 70)

# Load ASA optimization data
load("data/asa_optimization_data.rda")

##-----------------------------------------------------------------------------
# 2. Make arrivals histograms -- LOW BALKING
##-----------------------------------------------------------------------------

# Figures for slides
optim_df <- GetArrivalsDF(low, probmatrix, alpha = 0.01, beta = 0.02, mu = 0.5,
                          t = 960)
PlotArrivals(optim_df, ymax = 502)
PlotArrivals(optim_df, ymax = 502, grey = T)

naive_df <- GetArrivalsDF(naive, probmatrix, alpha = 0.01, beta = 0.02, mu = 0.5,
                           t = 960)
PlotArrivals(naive_df, ymax = 502)
PlotArrivals(naive_df, ymax = 502, grey = T)

# Figures for manuscript
pdf("_figures manuscript/R_output/arrivalshistograms_lowattrition.pdf", 
    width = 7, height = 6)
PlotArrivals(optim_df, ymax = 502)
PlotArrivals(naive_df, ymax = 502)
dev.off()

##-----------------------------------------------------------------------------
# 3. Make arrivals histograms -- HIGH BALKING
##-----------------------------------------------------------------------------

# Figures for slides
optim_df <- GetArrivalsDF(high, probmatrix, alpha = 0.1, beta = 0.1, mu = 0.5,
                          t = 960)
PlotArrivals(optim_df, ymax = 502)
PlotArrivals(optim_df, ymax = 502, grey = T)

naive_df <- GetArrivalsDF(naive, probmatrix, alpha = 0.1, beta = 0.1, mu = 0.5,
                           t = 960)
PlotArrivals(naive_df, ymax = 502)
PlotArrivals(naive_df, ymax = 502, grey = T)

# Figures for manuscript
pdf("_figures manuscript/R_output/arrivalshistograms_highattrition.pdf", 
    width = 7, height = 6)
PlotArrivals(optim_df, ymax = 502)
PlotArrivals(naive_df, ymax = 502)
dev.off()

