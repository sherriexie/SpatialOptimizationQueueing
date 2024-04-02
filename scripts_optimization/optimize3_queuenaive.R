library(tidyverse)
library(sf)
library(analyze.stuff)
library(tictoc)

#-----------------------------------------------------------------------------
# 1. Load data and functions
##-----------------------------------------------------------------------------

# Load optimization data for ASA
load("asa_optimization_data.rda")

# Load GA and RI functions
source("functions_recursiveinterchange.R")
source("functions_geneticalgorithm.R")

# Load solutions from first round of RI optimization
solns <- read_rds("solns/asaT16_queuenaive_round2.rds")
supplysets <- solns$supplysets
scores <- solns$scores

# Set random seed
set.seed(345)

##-----------------------------------------------------------------------------
# 2a. Create 500 offspring sets (to serve as starting sets) by mating pairs of 
#     optimized parents
##-----------------------------------------------------------------------------

newsupplysets1 <- list()
for(i in 1:500){
  newsupplysets1[[i]] <- PerformMating(supplysets[[2*i-1]], supplysets[[2*i]],
                                       zones, totalsites = 70, k = 20, 
                                       prob = 0.1)
}

##-----------------------------------------------------------------------------
# 2b. Create 500 offspring sets (to serve as starting sets) by mating an 
#     optimized parent with a random parent
##-----------------------------------------------------------------------------

# Create random parents using the RandomSelection function
randomsupplysets <- list()
for(i in 1:500){
  randomsupplysets[[i]] <- RandomSelection(zones, minsites = 4, k = 20, 
                                           totalsites = 70)
}

# Create offspring with one optimized parent and one randomly selected parent
newsupplysets2 <- list()
for(i in 1:500){
  newsupplysets2[[i]] <- PerformMating(supplysets[[i]], randomsupplysets[[i]],
                                       zones, totalsites = 70, k = 20, 
                                       prob = 0.1)
}

##-----------------------------------------------------------------------------
# 3. Perform recursive interchange optimization using 1000 offspring sets as
#    the initial starting set 
##-----------------------------------------------------------------------------

offspring <- append(newsupplysets1, newsupplysets2)

# Get 1000 sets of solutions using 1000 different seeds
n_iterations <- 1000
rand_seed <- sample(9999999, n_iterations)
solutions <- list()

# Run algorithm over all iterations
tic()
for (i in 1:n_iterations){  
  print(paste0("============ Iteration ", i, " of ", n_iterations, " ============"))
  set.seed(rand_seed[i])
  solutions[[i]] <- RunOptimizerNaiveSP(probmatrix, startsupply = offspring[[i]], 
                                        mu = 0.5, t = 960, habitability = 0.57, 
                                        hh_withdogs = 0.4)
}
toc()

# Format solutions and save
soln_out <- FormatSolutions(solutions)
saveRDS(soln_out, file = "asaT16_queuenaive_round3.rds")