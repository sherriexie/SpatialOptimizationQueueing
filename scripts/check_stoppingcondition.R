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

##-----------------------------------------------------------------------------
# 4. Check stopping condition for minimum sites needed - low attrition
##-----------------------------------------------------------------------------

# k = 25
k25_low <- read_rds("asaT16_lowatt_k25_round1.rds")
GetBestSolution(k25_low)   # Best score: 6475.8681141456 (59.997%)

k25_low2 <- read_rds("asaT16_lowatt_k25_round2.rds")
GetBestSolution(k25_low2)   # Best score: 6476.03254892421 (59.998%)

##-----------------------------------------------------------------------------
# 5. Check stopping condition for minimum sites needed - high attrition
##-----------------------------------------------------------------------------

# k = 25
k25_high <- read_rds("asaT16_highatt_k25_round1.rds")
GetBestSolution(k25_high)   # Best score: 5569.20007311315 (51.597%)
k25_high2 <- read_rds("asaT16_highatt_k25_round2.rds")
GetBestSolution(k25_high2)  # Best score: 5569.20007311315 (51.597%)

# k = 30
k30_high <- read_rds("asaT16_highatt_k30_round1.rds")
GetBestSolution(k30_high)   # Best score: 5847.73759082076 (54.178%)
k30_high2 <- read_rds("asaT16_highatt_k30_round2.rds")
GetBestSolution(k30_high2)  # Best score: 5848.08285387748 (54.181%)

# k = 35
k35_high <- read_rds("asaT16_highatt_k35_round1.rds")
GetBestSolution(k35_high)   # Best score: 6054.10868783491 (56.089%)
k35_high2 <- read_rds("asaT16_highatt_k35_round2.rds")
GetBestSolution(k35_high2)  # Best score: 6054.10868783491 (56.089%)

# k = 40
k40_high <- read_rds("asaT16_highatt_k40_round1.rds")
GetBestSolution(k40_high)   # Best score: 6197.4266961641 (57.417%)
k40_high2 <- read_rds("asaT16_highatt_k40_round2.rds")
GetBestSolution(k40_high2)  # Best score: 6197.4266961641 (57.417%)

# k = 45
k45_high <- read_rds("asaT16_highatt_k45_round1.rds")
GetBestSolution(k45_high)   # Best score: 6314.8933355322 (58.505%)
k45_high2 <- read_rds("asaT16_highatt_k45_round2.rds")
GetBestSolution(k45_high2)  # Best score: 6314.8933355322 (58.505%)

# k = 50
k50_high <- read_rds("asaT16_highatt_k50_round1.rds")
GetBestSolution(k50_high)   # Best score: 6394.84928232098 (59.246%)
k50_high2 <- read_rds("asaT16_highatt_k50_round2.rds")
GetBestSolution(k50_high2)  # Best score: 6394.84928232098 (59.246%)

# k = 55
k55_high <- read_rds("asaT16_highatt_k55_round1.rds")
GetBestSolution(k55_high)   # Best score: 6454.40526230431 (59.798%)
k55_high2 <- read_rds("asaT16_highatt_k55_round2.rds")
GetBestSolution(k55_high2)  # Best score: 6454.40526230431 (59.798%)


##-----------------------------------------------------------------------------
# 5. Old minimum number of sites needed to reach 60% vax covearge
##-----------------------------------------------------------------------------

# Low attrition, k = 21
k21_low <- read_rds("asaT16_lowatt_k21_round1.rds")
GetBestSolution(k21_low)   # Best score: 6244.52087453462 (57.853%)

# Low attrition, k = 22
k22_low <- read_rds("asaT16_lowatt_k22_round1.rds")
GetBestSolution(k22_low)   # Best score: 6308.91651435958 (58.450%)

# Low attrition, k = 23
k23_low <- read_rds("asaT16_lowatt_k23_round1.rds")
GetBestSolution(k23_low)   # Best score: 6372.31931177635 (59.0375%)

# Low attrition, k = 24
k24_low <- read_rds("asaT16_lowatt_k24_round1.rds")
GetBestSolution(k24_low)   # Best score: 6427.42872413128 (59.548%)

# High attrition, k = 56
k56_high <- read_rds("asaT16_highatt_k56_round1.rds")
GetBestSolution(k56_high)   # Best score: 6464.12983203772 (59.888%)

# High attrition, k = 57
k57_high <- read_rds("asaT16_highatt_k57_round1.rds")
GetBestSolution(k57_high)   # Best score: 6472.28976497437 (59.964%) 

# High attrition, k = 58
k58_high <- read_rds("asaT16_highatt_k58_round1.rds")
GetBestSolution(k58_high)   # Best score: 6478.50516213892 (60.021%)