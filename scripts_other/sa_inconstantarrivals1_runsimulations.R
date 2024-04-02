library(tidyverse)

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

SimulateBalkingReneging <- function(arrdensity, n, mu, alpha, beta, t=240, 
                                    fulloutput = F){
  
  lambdas <- arrdensity*n
  
  # Generate times of arrivals for 30 minute intervals
  at_list <- NA
  # Ensure the number of dogs arriving is ~n
  while(length(unlist(at_list)) < 0.95*n | length(unlist(at_list)) > 1.05*n){
    at_list <- sapply(lambdas, function(x) GenerateArrivalTimes(lambda = x, total_t = 30))
  }
  arrival_times <- c(at_list[[1]], at_list[[2]] + 30, at_list[[3]] + 60, 
                     at_list[[4]] + 90, at_list[[5]] + 120, at_list[[6]] + 150, 
                     at_list[[7]] + 180, at_list[[8]] + 210)
 
  # Generate service times
  # - Service times can be drawn from exp(mu) -- note not all times will be used
  #   due to balking
  service_times <- rexp(length(arrival_times), mu)
  
  # Keep a vector of time that events occur, initially populated only by the 
  # first arrival
  times_of_events <- arrival_times[1]
  
  # Keep a vector that tracks the type of event that occurs at each time in
  # times_of_events. We know that the first arrival will join the queue.
  types_of_events <- factor(NA, levels = c("joins", "balks", "reneges", "served"))
  types_of_events[1] <- "joins"
  
  # Keep a vector that tracks the queue length after each event occurs. Note the 
  # queue length is 1 after the first customer joins the queue
  queue_length <- numeric()
  queue_length[1] <- 1
  
  # Initialize trackers
  ct <- times_of_events  # ct (current time) tracks time in the while loop
  i <- 1  # i tracks number of events (we start with 1 because we already have
  # a customer who has joined the queue)
  j <- 1  # j tracks number of arrivals
  k <- 0  # k tracks number of households served
  tt_next_arrival <- arrival_times[j+1] - ct  # time 'til next arrival
  tt_next_service <- service_times[k+1]  # time 'til next service
  renege_time <- NA  # renege_time only calculated when queue_length > 1
  
  # Each loop represents one event (customer joins, balks, reneges, or served)
  # Continue loop until the current time exceeds the total time t AND the 
  # renege_time, tt_next_arrival, and tt_next_service are not ALL NA
  while (ct < t & sum(is.na(c(renege_time, tt_next_arrival, tt_next_service))) < 3){
    
    # If the renege time exists and is smaller than the time 'til next arrival 
    # or service then update accordingly
    if(which.min(c(renege_time, tt_next_arrival, tt_next_service)) == 1){
      i <- i+1
      times_of_events[i] <- ct + renege_time
      ct <- ct + renege_time
      types_of_events[i] <- "reneges"
      queue_length[i] <- queue_length[i-1] - 1
      tt_next_arrival <- tt_next_arrival - renege_time
      if(queue_length[i] == 0){  # if no one left in queue - service time resets
        tt_next_service <- NA
      }else{
        tt_next_service <- tt_next_service - renege_time
      }
      if(queue_length[i] > 1){
        renege_time <- rexp(1, (queue_length[i]-1)*beta)
      }else{
        renege_time <- NA
      }
    # Check to see if next event is an arrival and if so, whether the arriving
    # owner JOINS or BALKS  
    }else if(which.min(c(renege_time, tt_next_arrival, tt_next_service)) == 2){
      # Flip a coin with p = exp(-alpha*n/mu) to determine if arrival joins
      # or balks
      if(runif(1) < exp(-alpha*queue_length[i]/mu)){  # customer JOINS
        i <- i+1
        j <- j+1
        times_of_events[i] <- ct + tt_next_arrival
        ct <- ct + tt_next_arrival
        types_of_events[i] <- "joins"
        queue_length[i] <- queue_length[i-1] + 1
        if(is.na(tt_next_service)){
          tt_next_service <- service_times[k]
        }else{
          tt_next_service <- tt_next_service - tt_next_arrival
        }
        tt_next_arrival <- arrival_times[j+1] - ct
        if(queue_length[i] > 1){  # Update renege time to reflect longer line
          renege_time <- rexp(1, (queue_length[i]-1)*beta)
        }else{
          renege_time <- NA
        }
      }else{  # customer BALKS
        i <- i+1
        j <- j+1
        times_of_events[i] <- ct + tt_next_arrival
        ct <- ct + tt_next_arrival
        types_of_events[i] <- "balks"
        queue_length[i] <- queue_length[i-1]
        if(queue_length[i] == 0){  # if no one left in queue - service time resets
          tt_next_service <- NA
        }else{
          tt_next_service <- tt_next_service - tt_next_arrival
        }
        tt_next_arrival <- arrival_times[j+1] - ct
        #if(!is.na(renege_time) & queue_length[i] > 1){  # Update renege time
        #  renege_time <- renege_time - tt_next_arrival
        #}else if(is.na(renege_time) & queue_length[i] > 1){
        #  renege_time <- rexp(1, (queue_length[i]-1)*beta)
        if(queue_length[i] > 1){ # Update renege time to reflect shorter queue length
          renege_time <- rexp(1, (queue_length[i]-1)*beta)
        }else{
          renege_time <- NA
        }
      }
      
    # Check to see if next event is a SERVICE and update accordingly
    }else if(which.min(c(renege_time, tt_next_arrival, tt_next_service)) == 3){
      i <- i+1
      k <- k+1
      times_of_events[i] <- ct + tt_next_service
      ct <- ct + tt_next_service
      types_of_events[i] <- "served"
      queue_length[i] <- queue_length[i-1] - 1
      tt_next_arrival <- tt_next_arrival - tt_next_service
      if(queue_length[i] == 0){  # if no one left in queue - service time resets
        tt_next_service <- NA
      }else{
        tt_next_service <- service_times[k]
      }
      #if(!is.na(renege_time) & queue_length[i] > 1){  # Update renege time
      #  renege_time <- renege_time - tt_next_arrival
      #}else if(is.na(renege_time) & queue_length[i] > 1){
      #  renege_time <- rexp(1, (queue_length[i]-1)*beta)
      if(queue_length[i] > 1){ # Update renege time to reflect shorter queue length
        renege_time <- rexp(1, (queue_length[i]-1)*beta)
      }else{
        renege_time <- NA
      }

    }
    
  }  
  
  if(tail(times_of_events, 1) > t){  # Remove elements that occur after time t
    if(tail(types_of_events, 1) == "served"){
      k <- k-1
    }
    times_of_events <- head(times_of_events, -1)
    types_of_events <- head(types_of_events, -1)
    queue_length <- head(queue_length, -1)
  }
  
  # Calculate number lost to attrition by summing the number balked and reneged
  num_attrition <- table(types_of_events)[2] + table(types_of_events)[3]
  
  return(list(num_vaxed = k, num_attrition = num_attrition))
}

GenerateArrivalTimes <- function(lambda, total_t, draws = 100){
  # Obtain interarrival times for rate = lambda
  intarr_times <- rexp(draws, lambda)
  arr_times <- cumsum(intarr_times)
  
  # Truncate arrival times to the set time duration (total_t)
  arr_times <- arr_times[arr_times < total_t]
  return(arr_times)
}

##-----------------------------------------------------------------------------
# 1. Load arrival density distributions
##-----------------------------------------------------------------------------

#setwd("~/Optimization/")
#arrdensity_list <- read_rds("data/arrivalratedistributions.rds")

setwd("/project/cricardolab")
arrdensity_list <- read_rds("arrivalratedistributions.rds")

arrdensity_a <- arrdensity_list$dist_a
arrdensity_b <- arrdensity_list$dist_b
arrdensity_c <- arrdensity_list$dist_c
arrdensity_d <- arrdensity_list$dist_d

##-----------------------------------------------------------------------------
# 2. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION A - low attrition
##-----------------------------------------------------------------------------

# More support for smaller n values  
n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
length(n_range)  # 90
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                      length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_a, n = n_i, mu = 0.5, 
                                       alpha = 0.01, beta = 0.02)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityA_lowattrition.rds")
2+2
 
##-----------------------------------------------------------------------------
# 3. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION B - low attrition
##-----------------------------------------------------------------------------

# More support for smaller n values
n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
length(n_range)  # 90
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_b, n = n_i, mu = 0.5, 
                                       alpha = 0.01, beta = 0.02)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityB_lowattrition.rds")
2+2

##-----------------------------------------------------------------------------
# 4. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION C - low attrition
##-----------------------------------------------------------------------------

n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
length(n_range)  # 90
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_c, n = n_i, mu = 0.5, 
                                       alpha = 0.01, beta = 0.02)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityC_lowattrition.rds")
2+2

##-----------------------------------------------------------------------------
# 5. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION D - low attrition
##-----------------------------------------------------------------------------

n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_d, n = n_i, mu = 0.5, 
                                       alpha = 0.01, beta = 0.02)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityD_lowattrition.rds")
2+2


##-----------------------------------------------------------------------------
# 6. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION A - high attrition
##-----------------------------------------------------------------------------

n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_a, n = n_i, mu = 0.5, 
                                       alpha = 0.1, beta = 0.1)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityA_highattrition.rds")
2+2

##-----------------------------------------------------------------------------
# 7. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION B - high attrition
##-----------------------------------------------------------------------------


n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_b, n = n_i, mu = 0.5, 
                                       alpha = 0.1, beta = 0.1)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityB_highattrition.rds")
2+2

##-----------------------------------------------------------------------------
# 8. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION C - high attrition
##-----------------------------------------------------------------------------

n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_c, n = n_i, mu = 0.5, 
                                       alpha = 0.1, beta = 0.1)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityC_highattrition.rds")
2+2

##-----------------------------------------------------------------------------
# 9. Run over many iterations to estimate trend in behavior
#    DISTRIBUTION D - high attrition
##-----------------------------------------------------------------------------

n_range <- c(seq(2, 100, by = 2), seq(104, 200, by = 4), seq(220, 500, by = 20))
reps <- 200

# Note each rep is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be reps*4 runs per n value
nrep_df <- data.frame(n = rep(n_range, each = reps*4), 
                      day = rep(1:(4*reps), length(n_range)),
                      rep = rep(rep(1:reps, each = 4), # 4 days per rep
                                length(n_range)),  
                      num_vaxed = NA, num_attrition = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*reps)) {
    out_day <- SimulateBalkingReneging(arrdensity = arrdensity_d, n = n_i, mu = 0.5, 
                                       alpha = 0.1, beta = 0.1)  # Simulate a single day of the campaign
    nrep_df$num_vaxed[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_vaxed
    nrep_df$num_attrition[nrep_df$n == n_i & nrep_df$day == day] <- out_day$num_attrition
  }
}

nrep_df <- write_rds(nrep_df, "simulations_arrdensityD_highattrition.rds")
2+2