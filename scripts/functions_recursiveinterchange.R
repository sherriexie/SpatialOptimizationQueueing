# Running recursive interchange optimizer for the case where alpha and beta > 0
# and beginning with a randomly selected starting set of vax sites
RunOptimizer <- function(probmatrix, k, alpha, beta, mu, t, habitability, 
                         hh_withdogs) {
  
  # Obtain number of candidate vax sites and total dogs in the study area
  # **Note the ASA estimate for the dogs per dog-owning household is 1.86**
  n_supply <- ncol(probmatrix)
  tot_dogs <- nrow(probmatrix)*habitability*hh_withdogs*1.86
  
  # Generate random sample of k vax sites and calculate expected number 
  # vaccinated
  supplyset <- sample(c(1:n_supply), k)
  pm_set <- probmatrix[, supplyset]  # Prob matrix with selected vax sites
  num_vaxed <- EstimateExpVaxed(pm_set, habitability, hh_withdogs, alpha, beta, 
                                mu, t)
  vax_coverage <- num_vaxed/tot_dogs*100
  print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
               "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  
  # For loop across each selected site
  for (i in 1:length(supplyset)){
    
    # Other sites to try are all those that are not in the selected set
    othersites <- c(1:n_supply)[! c(1:n_supply) %in% supplyset]
    
    # For loop for the i'th site to rotate through all other sites to try
    for (j in othersites){
      # Replace i'th site with site j
      supplyset_tmp <- supplyset
      supplyset_tmp[i] <- j  
      pm_tmp <- probmatrix[, supplyset_tmp]  # Prob matrix with new set
      
      # Calculate expected number vaccinated for this new set
      numvaxed_tmp <- EstimateExpVaxed(pm_tmp, habitability, hh_withdogs, 
                                       alpha, beta, mu, t)
      
      # If replacing i'th site with j increases exp. num vaxed then replace
      if(numvaxed_tmp > num_vaxed){  
        num_vaxed <- numvaxed_tmp
        vax_coverage <- num_vaxed/tot_dogs*100
        supplyset <- supplyset_tmp
      }
    }
    # Print the optimal supply set and expected number vaxed for each iteration 
    # of the algorithm
    print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
                 "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  }
  
  # Save the names of the optimized set of vaccination points
  supplynames <- colnames(probmatrix)[supplyset]
  
  # Return a list of: the optimal supply set given by index (supply) and name 
  # (supplynames), as well as the expected number vaxed for the optimal set
  return(list(supply = supplyset, supplynames = supplynames, score = num_vaxed,
              vax_coverage = vax_coverage))
}

# Same as RunOptimizer but specify starting set of vax sites (used in 
# conjunction with the genetic algorithm)
RunOptimizerSelectedPoints <- function(probmatrix, startsupply, alpha, beta, 
                                       mu, t, habitability, hh_withdogs) {
  
  # Obtain number of candidate vax sites and total dogs in the study area
  n_supply <- ncol(probmatrix)
  tot_dogs <- nrow(probmatrix)*habitability*hh_withdogs*1.86
  
  # Start with the given set of vax sites
  supplyset <- startsupply
  pm_set <- probmatrix[, supplyset]  # Prob matrix with selected vax sites
  num_vaxed <- EstimateExpVaxed(pm_set, habitability, hh_withdogs, alpha, beta, 
                                mu, t)
  vax_coverage <- num_vaxed/tot_dogs*100
  print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
               "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  
  # For loop across each selected site
  for (i in 1:length(supplyset)){
    
    # Other sites to try are all those that are not in the selected set
    othersites <- c(1:n_supply)[! c(1:n_supply) %in% supplyset]
    
    # For loop for the i'th site to rotate through all other sites to try
    for (j in othersites){
      # Replace i'th site with site j
      supplyset_tmp <- supplyset
      supplyset_tmp[i] <- j  
      pm_tmp <- probmatrix[, supplyset_tmp]  # Prob matrix with new set
      
      # Calculate expected number vaccinated for this new set
      numvaxed_tmp <- EstimateExpVaxed(pm_tmp, habitability, hh_withdogs, 
                                       alpha, beta, mu, t)
      
      # If replacing i'th site with j increases exp. num vaxed then replace
      if(numvaxed_tmp > num_vaxed){  
        num_vaxed <- numvaxed_tmp
        vax_coverage <- num_vaxed/tot_dogs*100
        supplyset <- supplyset_tmp
      }
    }
    # Print the optimal supply set and expected number vaxed for each iteration 
    # of the algorithm
    print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
                 "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  }
  
  # Save the names of the optimized set of vaccination points
  supplynames <- colnames(probmatrix)[supplyset]
  
  # Return a list of: the optimal supply set given by index (supply) and name 
  # (supplynames), as well as the expected number vaxed for the optimal set
  return(list(supply = supplyset, supplynames = supplynames, score = num_vaxed,
              vax_coverage = vax_coverage))
}

# Same as RunOptimizer but draw a spatially structured starting set
RunOptimizerSpatial <- function(probmatrix, k, alpha, beta, mu, t, habitability, 
                               hh_withdogs) {
  
  # Obtain number of candidate vax sites and total dogs in the study area
  n_supply <- ncol(probmatrix)
  tot_dogs <- nrow(probmatrix)*habitability*hh_withdogs*1.86
  
  # Generate random sample of k vax sites and calculate expected number 
  # vaccinated
  supplyset <- RandomSelection(k = k)
  pm_set <- probmatrix[, supplyset]  # Prob matrix with selected vax sites
  num_vaxed <- EstimateExpVaxed(pm_set, habitability, hh_withdogs, alpha, beta, 
                                mu, t)
  vax_coverage <- num_vaxed/tot_dogs*100
  print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
               "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  
  # For loop across each selected site
  for (i in 1:length(supplyset)){
    
    # Other sites to try are all those that are not in the selected set
    othersites <- c(1:n_supply)[! c(1:n_supply) %in% supplyset]
    
    # For loop for the i'th site to rotate through all other sites to try
    for (j in othersites){
      # Replace i'th site with site j
      supplyset_tmp <- supplyset
      supplyset_tmp[i] <- j  
      pm_tmp <- probmatrix[, supplyset_tmp]  # Prob matrix with new set
      
      # Calculate expected number vaccinated for this new set
      numvaxed_tmp <- EstimateExpVaxed(pm_tmp, habitability, hh_withdogs, 
                                       alpha, beta, mu, t)
      
      # If replacing i'th site with j increases exp. num vaxed then replace
      if(numvaxed_tmp > num_vaxed){  
        num_vaxed <- numvaxed_tmp
        supplyset <- supplyset_tmp
        vax_coverage <- num_vaxed/tot_dogs*100
      }
    }
    # Print the optimal supply set and expected number vaxed for each iteration 
    # of the algorithm
    print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
                 "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  }
  
  # Save the names of the optimized set of vaccination points
  supplynames <- colnames(probmatrix)[supplyset]
  
  # Return a list of: the optimal supply set given by index (supply) and name 
  # (supplynames), as well as the expected number vaxed for the optimal set
  return(list(supply = supplyset, supplynames = supplynames, score = num_vaxed,
              vax_coverage = vax_coverage))
}

# Perform queue-naive optimization (recursive interchange without balking or reneging)
RunOptimizerNaive <- function(probmatrix, k, mu, t, habitability, hh_withdogs){
  
  # Obtain number of candidate vax sites and total dogs in the study area
  n_supply <- ncol(probmatrix)
  tot_dogs <- nrow(probmatrix)*habitability*hh_withdogs*1.86
  
  # Generate random sample of k vax sites and calculate expected number 
  # vaccinated
  supplyset <- sample(c(1:n_supply), k)  
  pm_set <- probmatrix[, supplyset]  # Probability matrix with selected vax points
  num_vaxed <- NaiveEstExpVaxed(pm_set, habitability, hh_withdogs, mu, t)
  vax_coverage <- num_vaxed/tot_dogs*100
  print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
              "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  
  # For loop across each selected site
  for (i in 1:length(supplyset)){
    
    # Other sites to try are all those that are not in the selected set
    othersites <- c(1:n_supply)[! c(1:n_supply) %in% supplyset]
    
    # For loop for the i'th site to rotate through all other sites to try
    for (j in othersites){
      # Replace i'th site with site j
      supplyset_tmp <- supplyset
      supplyset_tmp[i] <- j  
      pm_tmp <- probmatrix[, supplyset_tmp]  # Prob matrix with new set
      
      # Calculate expected number vaccinated for this new set
      numvaxed_tmp <- NaiveEstExpVaxed(pm_tmp, habitability, hh_withdogs, 
                                       mu, t)
      
      # If replacing i'th site with j increases exp. num vaxed then replace
      if(numvaxed_tmp > num_vaxed){  
        num_vaxed <- numvaxed_tmp
        vax_coverage <- num_vaxed/tot_dogs*100
        supplyset <- supplyset_tmp
      }
    }
    # Print the optimal supply set and expected number vaxed for each iteration 
    # of the algorithm
    print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
                 "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  }
  
  # Save the names of the optimized set of vaccination points
  supplynames <- colnames(probmatrix)[supplyset]
  
  # Return a list of: the optimal supply set given by index (supply) and name 
  # (supplynames), as well as the expected number vaxed for the optimal set
  return(list(supply = supplyset, supplynames = supplynames, score = num_vaxed,
              vax_coverage = vax_coverage))
}

# Same as RunOptimizerNaive but specify starting set of vax sites (used in 
# conjunction with the genetic algorithm)
RunOptimizerNaiveSP <- function(probmatrix, startsupply, mu, t, habitability, 
                                hh_withdogs){
  
  # Obtain number of candidate vax sites and total dogs in the study area
  n_supply <- ncol(probmatrix)
  tot_dogs <- nrow(probmatrix)*habitability*hh_withdogs*1.86
  
  # Start with the given set of vax sites
  supplyset <- startsupply
  pm_set <- probmatrix[, supplyset]  # Prob matrix with selected vax sites
  num_vaxed <- NaiveEstExpVaxed(pm_set, habitability, hh_withdogs, mu, t)
  vax_coverage <- num_vaxed/tot_dogs*100
  print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
               "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  
  # For loop across each selected site
  for (i in 1:length(supplyset)){
    
    # Other sites to try are all those that are not in the selected set
    othersites <- c(1:n_supply)[! c(1:n_supply) %in% supplyset]
    
    # For loop for the i'th site to rotate through all other sites to try
    for (j in othersites){
      # Replace i'th site with site j
      supplyset_tmp <- supplyset
      supplyset_tmp[i] <- j  
      pm_tmp <- probmatrix[, supplyset_tmp]  # Prob matrix with new set
      
      # Calculate expected number vaccinated for this new set
      numvaxed_tmp <- NaiveEstExpVaxed(pm_tmp, habitability, hh_withdogs, 
                                       mu, t)
      
      # If replacing i'th site with j increases exp. num vaxed then replace
      if(numvaxed_tmp > num_vaxed){  
        num_vaxed <- numvaxed_tmp
        vax_coverage <- num_vaxed/tot_dogs*100
        supplyset <- supplyset_tmp
      }
    }
    # Print the optimal supply set and expected number vaxed for each iteration 
    # of the algorithm
    print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
                 "E(num vaxed): ", num_vaxed, " (", vax_coverage, "%)"))
  }
  
  # Save the names of the optimized set of vaccination points
  supplynames <- colnames(probmatrix)[supplyset]
  
  # Return a list of: the optimal supply set given by index (supply) and name 
  # (supplynames), as well as the expected number vaxed for the optimal set
  return(list(supply = supplyset, supplynames = supplynames, score = num_vaxed,
              vax_coverage = vax_coverage))
}

# Select a spatially structured random set (at least 4 sites in each geographic 
# zone)
RandomSelection <- function(minsites = 4, k, totalsites = 70){
  
  # Define zones for potential vaccination sites
  z1 <- c(11, 13, 14, 15, 16, 19, 23, 24, 25, 26, 27, 28, 61, 62, 63, 70)
  z2 <- c(10, 12, 17, 18, 20, 21, 22, 29, 52, 54, 55, 56, 57, 58, 59, 60, 68, 69)
  z3 <- c(1, 2, 4, 5, 6, 7, 30, 31, 32, 33, 34, 35, 37, 47, 67)
  z4 <- c(3, 8, 9, 36, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 53, 
          64, 65, 66)
  
  # Select a minimum number of sites per zone
  sz1 <- sample(z1, minsites)
  sz2 <- sample(z2, minsites)
  sz3 <- sample(z3, minsites)
  sz4 <- sample(z4, minsites)
  initial_set <- c(sz1, sz2, sz3, sz4)
  
  # Sample from remaining sites
  allsites <- 1:totalsites
  remainingsites <- allsites[! allsites %in% initial_set]
  add_set <- sample(remainingsites, k-length(initial_set))
  final_set <- append(initial_set, add_set)
  
  # Return final, shuffled set
  sample(final_set, k)
  return(final_set)
}

# Find catchment size - number of dogs that are expected to show up with their
# owners at each vaccination site
FindCatchment <- function(pm_set, habitability, hh_withdogs){
  
  # Calculate number of houses that are expected to arrive at each site
  # - Find the closest site for each house (i.e., the site for which they have 
  #   highest arrival prob)
  site_index <- apply(pm_set, MARGIN = 1, which.max)
  # - Obtain the probability of arriving to the closest site for each house
  arrive_prob <- apply(pm_set, MARGIN = 1, max)
  # - For each site, sum the arrival probabilities for all houses in its 
  #   catchment and thus the expected number of arrivals
  arrivals_bysite <- numeric()
  for(i in 1:ncol(pm_set)){
    arrivals_bysite[i] <- sum(arrive_prob[site_index == i])
  }
  names(arrivals_bysite) <- colnames(pm_set)
  
  # Scale number by habitability, dog ownership rate, and number of dogs per
  # households (NOTE ASA ESTIMATE FOR DOGS/HH IS 1.86****)
  num_dogs_arriving <- arrivals_bysite*habitability*hh_withdogs*1.86
  return(num_dogs_arriving)
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

CalculateV <- function(lambda, alpha, beta, mu){
  # Calculate p0 which we'll need to calculate p_n (which is in both terms of V)
  p0 <- Calculate_p0(lambda, alpha, beta, mu)
  
  # Calculate effective arrivals by passing n = 0, 1, 2, ..., 99 to Expr 2
  vec1 <- 0:99
  Expr2_out <- sapply(vec1, function(x) Expr2(lambda, alpha, beta, mu, p0, x))
  
  # Calculate average reneging rate by passing n = 1, 2, ..., 99 to Expr 3
  vec2 <- 1:99
  Expr3_out <- sapply(vec2, function(x) Expr3(lambda, alpha, beta, mu, p0, x))
  
  # V is equal to the effective arrivals - those that renege
  V <- sum(Expr2_out, na.rm = T) - sum(Expr3_out, na.rm = T)
  return(V)
  
}

# Calculate expected number vaxed for the optimizer
EstimateExpVaxed <- function(pm_set, habitability, hh_withdogs, alpha, beta, mu, t){
  
  # Use probability matrix with chosen supply set to generate NUMBER OF DOGS
  # arriving at each site
  catchment <- FindCatchment(pm_set, habitability, hh_withdogs)
  
  # Use catchments to calculate lambda (/min) for each vax site assuming 
  # uniform arrivals
  lambdas <- catchment/t
  
  # Calculate expected vaccination rate V for each catchment
  V_catchments <- sapply(lambdas, CalculateV, alpha = alpha, beta = beta, mu = mu)
  
  # Sum over lambda bars to get overall vaccination rate (/min)
  totvaxrate <- sum(V_catchments)
  
  # Return the total vaccinated over total operational time
  return(totvaxrate*t)
  
}

# Same as EstimateExpVaxed but for the queue-naive optimizer
NaiveEstExpVaxed <- function(pm_set, habitability, hh_withdogs, mu, t){
  
  # Use probability matrix with chosen supply set to generate NUMBER OF DOGS
  # arriving at each site
  catchment <- FindCatchment(pm_set, habitability, hh_withdogs)
  
  # There is no attrition in the naive optimizer so all arrivals are assumed
  # to get vaxed
  return(sum(catchment))
  
}

FormatSolutions <- function(solutions){
  scores <- numeric()
  supplysets <- list()
  supplynames <- list()
  vaxcoverage <- numeric()
  for (i in 1:length(solutions)){  # Parse solutions
    scores[i] <- solutions[[i]]$score
    supplysets[[i]] <- solutions[[i]]$supply
    supplynames[[i]] <- solutions[[i]]$supplynames
    vaxcoverage[i] <- solutions[[i]]$vax_coverage
  }
  out <- list(solutions = solutions, scores = scores, supplysets = supplysets, 
              supplynames = supplynames, vaxcoverage = vaxcoverage)
  return(out)
}