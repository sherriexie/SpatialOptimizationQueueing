# Define vax site zones to set spatial structure for crossover and random selection
z1 <- c(11, 13, 14, 15, 16, 19, 23, 24, 25, 26, 27, 28, 61, 62, 63, 70)
z2 <- c(10, 12, 17, 18, 20, 21, 22, 29, 52, 54, 55, 56, 57, 58, 59, 60, 68, 69)
z3 <- c(1, 2, 4, 5, 6, 7, 30, 31, 32, 33, 34, 35, 37, 47, 67)
z4 <- c(3, 8, 9, 36, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 53, 
        64, 65, 66)
zones <- list(z1 = z1, z2 = z2, z3 = z3, z4 = z4)

# For a given supply set and a list of zone divisions, divide sites into zones.
# Output is a list of sets grouped into vectors (z1, z2, z3, z4), each 
# containing all the sites for a given zone.
SplitZones <- function(set, zones){
  # Split set into zones
  sz1 <- sort(set[set %in% zones$z1])
  sz2 <- sort(set[set %in% zones$z2])
  sz3 <- sort(set[set %in% zones$z3])
  sz4 <- sort(set[set %in% zones$z4])
  
  # Return as a named list
  return(list(z1 = sz1, z2 = sz2, z3 = sz3, z4 = sz4))
}

PerformCrossover <- function(set1, set2, zones){
  # Split sets into zones 
  set1_zones <- SplitZones(set1, zones)
  set2_zones <- SplitZones(set2, zones)
  
  # Sample 2 of the 4 zones
  sz <- sample(4, 2)
  crossover_zones <- append(set1_zones[sz], set2_zones[-sz])
  
  # Return all chosen elements as a single vector
  return(sort(unname(unlist(crossover_zones))))
}

AddDelete <- function(set, totalsites, k){
  if(length(set) > k){  # Delete at random if set has too many sites
    newset <- sample(set, k)
  }else if(length(set) < k){  # Add at random if set has too few sites
    allsites <- 1:totalsites
    remainingsites <- allsites[! allsites %in% set]
    addset <- sample(remainingsites, k - length(set))
    newset <- append(set, addset)
  }else{  # If set has the right number of sites return the set
    newset <- set
  }
  return(sort(newset))
}

PerformMutation <- function(set, prob, totalsites){
  # Indicator variable for whether each site is mutated
  mutate <- runif(length(set)) < prob
  
  # If there is at least one site undergoing mutation sample from remaining 
  # sites to perform mutation
  num_mutate <- sum(mutate)
  if(num_mutate > 0){
    mutate_index <- which(mutate)
    allsites <- 1:totalsites
    for (msite in mutate_index){
      remainingsites <- allsites[! allsites %in% set]
      set[msite] <- sample(remainingsites, 1)
    }
  }
  return(sort(set))  
}

# This is the main function used to combine two parent sets to obtain an 
# an offspring. Inputs are as folows:
# -> set1, set2 = these are the two parent sets
# -> zones = zone assignments of all candidate vax sites in list format
# -> totalsites = total number of candidate sites
# -> p = number of selected sites
# -> prob = mutation probability
PerformMating <- function(set1, set2, zones, totalsites, k, prob){
  # 1. Perform crossover
  crossover_set <- PerformCrossover(set1, set2, zones) 
  # 2. Add/delete sites to ensure there are the right number of total sites
  if(length(crossover_set) != k){
    crossover_set <- AddDelete(crossover_set, totalsites, k)
  }
  # 3. Perform mutation
  mutated_set <- PerformMutation(crossover_set, prob, totalsites)
  # 4. Shuffle sites to prepare for recursive interchange algorithm
  final_set <- sample(mutated_set, k)
  return(final_set)
}

# Randomly select a set of vax sites while ensuring adequate spatial spread (at
# least a minimum number of sites per zone)
RandomSelection <- function(zones, minsites, k, totalsites){
  
  # Select a minimum number of site per zone
  sz1 <- sample(zones$z1, minsites)
  sz2 <- sample(zones$z2, minsites)
  sz3 <- sample(zones$z3, minsites)
  sz4 <- sample(zones$z4, minsites)
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