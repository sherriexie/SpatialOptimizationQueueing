library(tidyverse)
library(scales)

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

GetResults <- function(supplyset, probmatrix, alpha, beta, mu, t){
  
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
  
  # Calculate number vaccinated, number balked, and number reneged for each site
  num_vaxed <- V*t
  num_balked <- balking_rate*t
  num_reneged <- reneging_rate*t
  
  # Calculate total vaccinated and total attrition across all sites, as well as
  # overall vaccination coverage
  tot_vaxed <- sum(num_vaxed)
  tot_attrition <- sum(num_balked) + sum(num_reneged)
  tot_dogs <- nrow(probmatrix)*0.57*0.4*1.86
  vax_coverage <- tot_vaxed/tot_dogs*100
  
  # Print total attrition, total vaccinated and vaccination coverage
  print(noquote(paste0("Number lost to attrition = ", tot_attrition)))
  print(noquote(paste0("Number of dogs vaccinated = ", tot_vaxed)))
  print(noquote(paste0("Vaccination coverage = ", vax_coverage, "%")))
  
  return(list(tot_vaxed = tot_vaxed, tot_attrition = tot_attrition, 
              vax_coverage = vax_coverage))
}

##-----------------------------------------------------------------------------
# 1. Load data 
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
# 2. Calculate number vaxed and lost to attrition for 4 scenarios and 3 
#    solutions
##-----------------------------------------------------------------------------

# Scenario 1: low attrition (low balking + reneging)
# - Implement low attrition solution
low_low <- GetResults(low, probmatrix, alpha = 0.01, beta = 0.02, mu = 0.5,
                          t = 960)
# - Implement high attrition solution
low_high <- GetResults(high, probmatrix, alpha = 0.01, beta = 0.02, mu = 0.5,
                          t = 960)
# - Implement queue-naive solution
low_naive <- GetResults(naive, probmatrix, alpha = 0.01, beta = 0.02, mu = 0.5,
                          t = 960)

# Scenario 2: high attrition (high balking + reneging)
high_low <- GetResults(low, probmatrix, alpha = 0.1, beta = 0.1, mu = 0.5,
                          t = 960)
high_high <- GetResults(high, probmatrix, alpha = 0.1, beta = 0.1, mu = 0.5,
                            t = 960)
high_naive <- GetResults(naive, probmatrix, alpha = 0.1, beta = 0.1, mu = 0.5,
                           t = 960)

# Scenario 3: low balking and high reneging
lowhigh_low <- GetResults(low, probmatrix, alpha = 0.01, beta = 0.1, mu = 0.5,
                             t = 960)
lowhigh_high <- GetResults(high, probmatrix, alpha = 0.01, beta = 0.1, mu = 0.5,
                           t = 960)
lowhigh_naive <- GetResults(naive, probmatrix, alpha = 0.01, beta = 0.1, mu = 0.5,
                            t = 960)

# Scenario 4: high balking and low reneging
highlow_low <- GetResults(low, probmatrix, alpha = 0.1, beta = 0.02, mu = 0.5,
                             t = 960)
highlow_high <- GetResults(high, probmatrix, alpha = 0.1, beta = 0.02, mu = 0.5,
                              t = 960)
highlow_naive <- GetResults(naive, probmatrix, alpha = 0.1, beta = 0.02, mu = 0.5,
                               t = 960)

# Compile data frame of results
results_df <- data.frame(scenario = c(rep("low attrition", 2), 
                                      rep("high attrition", 2),
                                      rep("low-high", 2), rep("high-low", 2)),
                         solution = rep(c("low", "high"), 4),
                         vaxed = c(low_low$tot_vaxed - low_naive$tot_vaxed,
                                   low_high$tot_vaxed - low_naive$tot_vaxed,
                                   high_low$tot_vaxed - high_naive$tot_vaxed,
                                   high_high$tot_vaxed - high_naive$tot_vaxed,
                                   lowhigh_low$tot_vaxed - lowhigh_naive$tot_vaxed,
                                   lowhigh_high$tot_vaxed - lowhigh_naive$tot_vaxed,
                                   highlow_low$tot_vaxed - highlow_naive$tot_vaxed,
                                   highlow_high$tot_vaxed - highlow_naive$tot_vaxed),
                         attrition = c(low_low$tot_attrition - low_naive$tot_attrition,
                                       low_high$tot_attrition - low_naive$tot_attrition,
                                       high_low$tot_attrition - high_naive$tot_attrition,
                                       high_high$tot_attrition - high_naive$tot_attrition,
                                       lowhigh_low$tot_attrition - lowhigh_naive$tot_attrition,
                                       lowhigh_high$tot_attrition - lowhigh_naive$tot_attrition,
                                       highlow_low$tot_attrition - highlow_naive$tot_attrition,
                                       highlow_high$tot_attrition - highlow_naive$tot_attrition))

# Convert to ordered factor variable to control order in plotting
results_df$scenario <- factor(results_df$scenario, ordered = T,
                              levels = c("low attrition" , "high attrition",
                                         "low-high", "high-low"))
results_df$solution <- factor(results_df$solution, ordered = T, 
                              levels = c("low", "high"))
##-----------------------------------------------------------------------------
# 3. Plot additional dogs vaccinated beyond the number achieved under the 
#    queue-naive solution for all 4 scenarios
##-----------------------------------------------------------------------------

pdf("_figures manuscript/R_output/sa_balkingrenegingparameter_dogsvaxed.pdf",
    height = 6, width = 12)
ggplot(results_df) +
  geom_bar(aes(x = scenario, y = vaxed, fill = solution), stat = "identity", 
           width = .6, color = "black", position = position_dodge(.7)) +
  scale_fill_manual(values = c("#93c2db", "#5794b5"), 
                    name = "MDVC sites \nemployed") +
  geom_hline(yintercept = 0, size = 1.1) +
  xlab("Scenario") +
  ylab("Additional dogs vaccinated*") +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22, face = "bold"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=22, face = "bold"),
        legend.spacing = unit(1, "cm")) +
  guides(fill = guide_legend(byrow = TRUE))
dev.off()

##-----------------------------------------------------------------------------
# 4. Plot reduction in the number of dogs lost to attrition
##-----------------------------------------------------------------------------

pdf("_figures manuscript/R_output/sa_balkingrenegingparameter_attrition.pdf",
    height = 6, width = 12)
ggplot(results_df) +
  geom_bar(aes(x = scenario, y = attrition, fill = solution), stat = "identity", 
           width = .6, color = "black", position = position_dodge(.7)) +
  scale_fill_manual(values = c("#eb7878", "#C00000"), 
                    name = "MDVC sites \nemployed") +
  geom_hline(yintercept = 0, size = 1.1) +
  xlab("Scenario") +
  ylab("Reduction in attrition*") +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22, face = "bold"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=22, face = "bold"),
        legend.spacing = unit(1, "cm")) +
  guides(fill = guide_legend(byrow = TRUE))

dev.off()
