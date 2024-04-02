library(tidyverse)
library(sf)
library(analyze.stuff)
library(tictoc)

#-----------------------------------------------------------------------------
# 1. Load data and functions
##-----------------------------------------------------------------------------

load("asa_optimization_data.rda")
source("functions_recursiveinterchange.R")

##-----------------------------------------------------------------------------
# 2. Optimize placement while incorporating balking and reneging
##-----------------------------------------------------------------------------

set.seed(123)

# Get 1000 sets of solutions using 1000 different seeds
n_iterations <- 1000
rand_seed <- sample(9999999, n_iterations)
solutions <- list()

# First 200 runs
tic()
for (i in 1:200){  # Run algorithm over iterations
  print(paste0("============ Iteration ", i, " of ", n_iterations, " ============"))
  set.seed(rand_seed[i])
  solutions[[i]] <- RunOptimizer(probmatrix, k = 20, alpha = 0.01, beta = 0.02, 
                                 mu = 0.5, t = 960, habitability = 0.57, 
                                 hh_withdogs = 0.4)
}
toc()

# Format solutions and save
soln_out <- FormatSolutions(solutions)
saveRDS(soln_out, file = "asaT16_a0.01_b0.02_1to200.rds")

# Runs 201-500
tic()
for (i in 201:500){  # Run algorithm over iterations
  print(paste0("============ Iteration ", i, " of ", n_iterations, " ============"))
  set.seed(rand_seed[i])
  solutions[[i]] <- RunOptimizer(probmatrix, k = 20, alpha = 0.01, beta = 0.02, 
                                 mu = 0.5, t = 960, habitability = 0.57, 
                                 hh_withdogs = 0.4)
}
toc()

# Format solutions and save
soln_out <- FormatSolutions(solutions)
saveRDS(soln_out, file = "asaT16_a0.01_b0.02_1to500.rds")

# Runs 501-1000
tic()
for (i in 501:n_iterations){  # Run algorithm over iterations
  print(paste0("============ Iteration ", i, " of ", n_iterations, " ============"))
  set.seed(rand_seed[i])
  solutions[[i]] <- RunOptimizer(probmatrix, k = 20, alpha = 0.01, beta = 0.02, 
                                 mu = 0.5, t = 960, habitability = 0.57, 
                                 hh_withdogs = 0.4)
}
toc()

# Format solutions and save
soln_out <- FormatSolutions(solutions)
saveRDS(soln_out, file = "asaT16_a0.01_b0.02_round1.rds")