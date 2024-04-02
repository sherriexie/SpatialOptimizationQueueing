library(tidyverse)
library(splines)
library(scales)

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

GetResults <- function(supplyset, probmatrix, splines){
  
  spline_vax <- splines$spline_vax
  spline_att <- splines$spline_att
  
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
  
  # Use splines to obtian the number of dogs vaccination and the number lost
  # to attrition
  vaxed_bysite <- predict(spline_vax, data.frame(num_arrivals = num_arrivals))
  attrition_bysite <- predict(spline_att, data.frame(num_arrivals = num_arrivals))
  
  num_vaxed <- sum(vaxed_bysite)
  num_attrition <- sum(attrition_bysite)
  
  # Print number vaccinated and number lost to attrition
  print(noquote(paste0("Number of dogs vaccinated = ", num_vaxed)))
  print(noquote(paste0("Number lost to attrition = ", num_attrition)))

  return(list(num_vaxed = num_vaxed, num_attrition = num_attrition))
}


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

GetResults2 <- function(supplyset, probmatrix, alpha, beta, mu, t){
  
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
  num_vaxed <- sum(num_vaxed)
  num_attrition <- sum(num_balked) + sum(num_reneged)
  tot_dogs <- nrow(probmatrix)*0.57*0.4*1.86

  # Print number vaccinated and number lost to attrition
  print(noquote(paste0("Number of dogs vaccinated = ", num_vaxed)))
  print(noquote(paste0("Number lost to attrition = ", num_attrition)))
  
  return(list(num_vaxed = num_vaxed, num_attrition = num_attrition))
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

# Load vax and attrition regression splines for arrival densities A-D and 
# assuming low and high attrition rates
splines_a_low <- read_rds("_output manuscript/splines_a_low.rds")
splines_a_high <- read_rds("_output manuscript/splines_a_high.rds")
splines_b_low <- read_rds("_output manuscript/splines_b_low.rds")
splines_b_high <- read_rds("_output manuscript/splines_b_high.rds")
splines_c_low <- read_rds("_output manuscript/splines_c_low.rds")
splines_c_high <- read_rds("_output manuscript/splines_c_high.rds")
splines_d_low <- read_rds("_output manuscript/splines_d_low.rds")
splines_d_high <- read_rds("_output manuscript/splines_d_high.rds")

##-----------------------------------------------------------------------------
# 2. Get results for each of the 8 scenarios with the queue-naive, low-, and 
#    high-attrition solutions
##-----------------------------------------------------------------------------

# A - low attrition
a_low_low <- GetResults(low, probmatrix, splines_a_low)
a_low_high <- GetResults(high, probmatrix, splines_a_low)
a_low_naive <- GetResults(naive, probmatrix, splines_a_low)

# A - high attrition
a_high_low <- GetResults(low, probmatrix, splines_a_high)
a_high_high <- GetResults(high, probmatrix, splines_a_high)
a_high_naive <- GetResults(naive, probmatrix, splines_a_high)

# B - low attrition
b_low_low <- GetResults(low, probmatrix, splines_b_low)
b_low_high <- GetResults(high, probmatrix, splines_b_low)
b_low_naive <- GetResults(naive, probmatrix, splines_b_low)

# B - high attrition
b_high_low <- GetResults(low, probmatrix, splines_b_high)
b_high_high <- GetResults(high, probmatrix, splines_b_high)
b_high_naive <- GetResults(naive, probmatrix, splines_b_high)

# C - low attrition
c_low_low <- GetResults(low, probmatrix, splines_c_low)
c_low_high <- GetResults(high, probmatrix, splines_c_low)
c_low_naive <- GetResults(naive, probmatrix, splines_c_low)

# C - high attrition
c_high_low <- GetResults(low, probmatrix, splines_c_high)
c_high_high <- GetResults(high, probmatrix, splines_c_high)
c_high_naive <- GetResults(naive, probmatrix, splines_c_high)

# D - low attrition
d_low_low <- GetResults(low, probmatrix, splines_d_low)
d_low_high <- GetResults(high, probmatrix, splines_d_low)
d_low_naive <- GetResults(naive, probmatrix, splines_d_low)

# D - high attrition
d_high_low <- GetResults(low, probmatrix, splines_d_high)
d_high_high <- GetResults(high, probmatrix, splines_d_high)
d_high_naive <- GetResults(naive, probmatrix, splines_d_high)

# Constant arrival rate - low attrition
const_low_low <- GetResults2(low, probmatrix, alpha = 0.01, beta = 0.02, mu = 0.5,
                             t = 960)
const_low_high <- GetResults2(high, probmatrix, alpha = 0.01, beta = 0.02, mu = 0.5,
                              t = 960)

# Constant arrival rate - high attrition
const_high_low <- GetResults2(low, probmatrix, alpha = 0.1, beta = 0.1, mu = 0.5,
                             t = 960)
const_high_high <- GetResults2(high, probmatrix, alpha = 0.1, beta = 0.1, mu = 0.5,
                              t = 960)

##-----------------------------------------------------------------------------
# 3. Assemble data frames for plotting
##-----------------------------------------------------------------------------

# Results data frame comparing number vaxed and lost to attrition for the 
# different inconstant arrival densities vs. constant arrivals
resultsvsconst_df <- data.frame(scenario = c(rep("A - low attrition", 2), 
                                             rep("A - high attrition", 2),
                                             rep("B - low attrition", 2), 
                                             rep("B - high attrition", 2),
                                             rep("C - low attrition", 2), 
                                             rep("C - high attrition", 2),
                                             rep("D - low attrition", 2), 
                                             rep("D - high attrition", 2)),
                                solution = rep(c("low", "high"), 8),
                                vaxed = c(a_low_low$num_vaxed - const_low_low$num_vaxed,
                                          a_low_high$num_vaxed - const_low_high$num_vaxed,
                                          a_high_low$num_vaxed - const_high_low$num_vaxed,
                                          a_high_high$num_vaxed - const_high_high$num_vaxed,
                                          b_low_low$num_vaxed - const_low_low$num_vaxed,
                                          b_low_high$num_vaxed - const_low_high$num_vaxed,
                                          b_high_low$num_vaxed - const_high_low$num_vaxed,
                                          b_high_high$num_vaxed - const_high_high$num_vaxed,
                                          c_low_low$num_vaxed - const_low_low$num_vaxed,
                                          c_low_high$num_vaxed - const_low_high$num_vaxed,
                                          c_high_low$num_vaxed - const_high_low$num_vaxed,
                                          c_high_high$num_vaxed - const_high_high$num_vaxed,
                                          d_low_low$num_vaxed - const_low_low$num_vaxed,
                                          d_low_high$num_vaxed - const_low_high$num_vaxed,
                                          d_high_low$num_vaxed - const_high_low$num_vaxed,
                                          d_high_high$num_vaxed - const_high_high$num_vaxed),
                                attrition = c(a_low_low$num_attrition - const_low_low$num_attrition,
                                              a_low_high$num_attrition - const_low_high$num_attrition,
                                              a_high_low$num_attrition - const_high_low$num_attrition,
                                              a_high_high$num_attrition - const_high_high$num_attrition,
                                              b_low_low$num_attrition - const_low_low$num_attrition,
                                              b_low_high$num_attrition - const_low_high$num_attrition,
                                              b_high_low$num_attrition - const_high_low$num_attrition,
                                              b_high_high$num_attrition - const_high_high$num_attrition,
                                              c_low_low$num_attrition - const_low_low$num_attrition,
                                              c_low_high$num_attrition - const_low_high$num_attrition,
                                              c_high_low$num_attrition - const_high_low$num_attrition,
                                              c_high_high$num_attrition - const_high_high$num_attrition,
                                              d_low_low$num_attrition - const_low_low$num_attrition,
                                              d_low_high$num_attrition - const_low_high$num_attrition,
                                              d_high_low$num_attrition - const_high_low$num_attrition,
                                              d_high_high$num_attrition - const_high_high$num_attrition))

# Convert to ordered factor variable to control order in plotting
resultsvsconst_df$scenario <- factor(resultsvsconst_df$scenario, ordered = T,
                                    levels = c("A - low attrition" , 
                                               "B - low attrition" , 
                                               "C - low attrition" , 
                                               "D - low attrition" , 
                                               "A - high attrition",
                                               "B - high attrition",
                                               "C - high attrition",
                                               "D - high attrition"))
resultsvsconst_df$solution <- factor(resultsvsconst_df$solution, ordered = T, 
                                     levels = c("low", "high"))


# Results data frame comparing queue-conscious vs. queue-naive solutions for 
# the eight scenarios
resultsvsnaive_df <- data.frame(scenario = c(rep("A - low attrition", 2), 
                                             rep("A - high attrition", 2),
                                             rep("B - low attrition", 2), 
                                             rep("B - high attrition", 2),
                                             rep("C - low attrition", 2), 
                                             rep("C - high attrition", 2),
                                             rep("D - low attrition", 2), 
                                             rep("D - high attrition", 2)),
                         solution = rep(c("low", "high"), 8),
                         vaxed = c(a_low_low$num_vaxed - a_low_naive$num_vaxed,
                                   a_low_high$num_vaxed - a_low_naive$num_vaxed,
                                   a_high_low$num_vaxed - a_high_naive$num_vaxed,
                                   a_high_high$num_vaxed - a_high_naive$num_vaxed,
                                   b_low_low$num_vaxed - b_low_naive$num_vaxed,
                                   b_low_high$num_vaxed - b_low_naive$num_vaxed,
                                   b_high_low$num_vaxed - b_high_naive$num_vaxed,
                                   b_high_high$num_vaxed - b_high_naive$num_vaxed,
                                   c_low_low$num_vaxed - c_low_naive$num_vaxed,
                                   c_low_high$num_vaxed - c_low_naive$num_vaxed,
                                   c_high_low$num_vaxed - c_high_naive$num_vaxed,
                                   c_high_high$num_vaxed - c_high_naive$num_vaxed,
                                   d_low_low$num_vaxed - d_low_naive$num_vaxed,
                                   d_low_high$num_vaxed - d_low_naive$num_vaxed,
                                   d_high_low$num_vaxed - d_high_naive$num_vaxed,
                                   d_high_high$num_vaxed - d_high_naive$num_vaxed),
                         attrition = c(a_low_low$num_attrition - a_low_naive$num_attrition,
                                       a_low_high$num_attrition - a_low_naive$num_attrition,
                                       a_high_low$num_attrition - a_high_naive$num_attrition,
                                       a_high_high$num_attrition - a_high_naive$num_attrition,
                                       b_low_low$num_attrition - b_low_naive$num_attrition,
                                       b_low_high$num_attrition - b_low_naive$num_attrition,
                                       b_high_low$num_attrition - b_high_naive$num_attrition,
                                       b_high_high$num_attrition - b_high_naive$num_attrition,
                                       c_low_low$num_attrition - c_low_naive$num_attrition,
                                       c_low_high$num_attrition - c_low_naive$num_attrition,
                                       c_high_low$num_attrition - c_high_naive$num_attrition,
                                       c_high_high$num_attrition - c_high_naive$num_attrition,
                                       d_low_low$num_attrition - d_low_naive$num_attrition,
                                       d_low_high$num_attrition - d_low_naive$num_attrition,
                                       d_high_low$num_attrition - d_high_naive$num_attrition,
                                       d_high_high$num_attrition - d_high_naive$num_attrition))
resultsvsnaive_df$scenario <- factor(resultsvsnaive_df$scenario, ordered = T,
                                     levels = c("A - low attrition" , 
                                                "B - low attrition" , 
                                                "C - low attrition" , 
                                                "D - low attrition" , 
                                                "A - high attrition",
                                                "B - high attrition",
                                                "C - high attrition",
                                                "D - high attrition"))
resultsvsnaive_df$solution <- factor(resultsvsnaive_df$solution, ordered = T, 
                                     levels = c("low", "high"))

##-----------------------------------------------------------------------------
# 4. Plot difference in number of dogs vaccinated compared to the constant 
#    arrival rate scenario
##-----------------------------------------------------------------------------

pdf("_figures manuscript/R_output/sa_constvsinconstantarrivalrate_dogsvaxed.pdf",
    height = 6, width = 12)
ggplot(resultsvsconst_df) +
  geom_bar(aes(x = scenario, y = vaxed, fill = solution), stat = "identity", 
           width = .6, color = "black", position = position_dodge(.7)) +
  scale_fill_manual(values = c("#93c2db", "#5794b5"), 
                    name = "MDVC sites \nemployed") +
  geom_hline(yintercept = 0, size = 1.1) +
  xlab("Scenario") +
  ylab("Difference in number of dogs vaccinated") +
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
# 5. Plot reduction in the number of dogs lost to attrition
##-----------------------------------------------------------------------------

pdf("_figures manuscript/R_output/sa_constvsinconstantarrivalrate_attrition.pdf",
    height = 6, width = 12)
ggplot(resultsvsconst_df) +
  geom_bar(aes(x = scenario, y = attrition, fill = solution), stat = "identity", 
           width = .6, color = "black", position = position_dodge(.7)) +
  scale_fill_manual(values = c("#eb7878", "#C00000"), 
                    name = "MDVC sites \nemployed") +
  geom_hline(yintercept = 0, size = 1.1) +
  xlab("Scenario") +
  ylab("Difference in number lost to attrition") +
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
# 6. Plot additional dogs vaccinated beyond the number achieved under the 
#    queue-naive solution for all 8 scenarios
##-----------------------------------------------------------------------------

pdf("_figures manuscript/R_output/sa_inconstantarrivalrate_dogsvaxed.pdf",
    height = 6, width = 12)
ggplot(resultsvsnaive_df) +
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
# 7. Plot reduction in the number of dogs lost to attrition
##-----------------------------------------------------------------------------

pdf("_figures manuscript/R_output/sa_inconstantarrivalrate_attrition.pdf",
    height = 6, width = 12)
ggplot(resultsvsnaive_df) +
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
