library(tidyverse)

setwd("~/Optimization/_output manuscript/")

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

GetBestSolution <- function(soln_list){
  
  # Obtain best score and associated vax coverage and set of vax sites
  best_index <- which.max(soln_list$scores)
  best_score <- soln_list$scores[best_index]
  best_vc <- soln_list$vaxcoverage[best_index]
  best_soln <- sort(soln_list$supplysets[[best_index]])
  
  # Print output
  print(noquote(paste0("Best score: ", best_score)))
  print(noquote(paste0("Max vax coverage: ", best_vc, "%")))
  print(noquote(paste0("Best solution: ", paste(sort(best_soln),
                                                collapse = ", "))))
  
  # Histogram of exp vax coverage
  hist(soln_list$vaxcoverage, breaks = 100, xlab = "Expected vax coverage")
  
}

##-----------------------------------------------------------------------------
# 1. Check stopping condition for queue-naive optimization, k = 20
##-----------------------------------------------------------------------------

k20_naive <- read_rds("asaT16_queuenaive_round1.rds")
k20_naive2 <- read_rds("asaT16_queuenaive_round2.rds")
k20_naive3 <- read_rds("asaT16_queuenaive_round3.rds")

GetBestSolution(k20_naive)  # Best score: 6824.01426381027 (63.222%)
GetBestSolution(k20_naive2)  # Best score: 6824.01426381027 (63.222%)
GetBestSolution(k20_naive3)  # Best score: 6824.01426381027 (63.222%)

##-----------------------------------------------------------------------------
# 2. Check stopping condition for low attrition optimization, k = 20
##-----------------------------------------------------------------------------

k20_low <- read_rds("asaT16_a0.01_b0.02_round1.rds")
k20_low2 <- read_rds("asaT16_a0.01_b0.02_round2.rds")
k20_low3 <- read_rds("asaT16_a0.01_b0.02_round3.rds")
k20_low4 <- read_rds("asaT16_a0.01_b0.02_round4.rds")
k20_low5 <- read_rds("asaT16_a0.01_b0.02_round5.rds")

GetBestSolution(k20_low)   # Best score: 6171.50907163532
GetBestSolution(k20_low2)  # Best score: 6173.91700878312 (57.199%)
GetBestSolution(k20_low3)  # Best score: 6174.45600700063 (57.204%)
GetBestSolution(k20_low4)  # Best score: 6174.45600700063 (57.204%)
GetBestSolution(k20_low5)  # Best score: 6174.45600700063 (57.204%)

##-----------------------------------------------------------------------------
# 3. Check stopping condition for high attrition optimization, k = 20
##-----------------------------------------------------------------------------

k20_high <- read_rds("asaT16_a0.1_b0.1_round1.rds")
k20_high2 <- read_rds("asaT16_a0.1_b0.1_round2.rds")
k20_high3 <- read_rds("asaT16_a0.1_b0.1_round3.rds")
k20_high4 <- read_rds("asaT16_a0.1_b0.1_round4.rds")

GetBestSolution(k20_high)   # Best score: 5180.58042488482 (47.996%)
GetBestSolution(k20_high2)  # Best score: 5182.81232497812 (48.017%)
GetBestSolution(k20_high3)  # Best score: 5182.81232497812 (48.017%)
GetBestSolution(k20_high4)  # Best score: 5182.81232497812 (48.017%)
