library(tidyverse)
library(scales)
library(gridExtra)
library(ggpubr)

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

SimulateBalkingReneging <- function(n, mu, alpha, beta, t, fulloutput = F){
  
  # Generate times of arrivals and service times
  # - Times of arrival can be drawn from a uniform distribution U(0, t)
  arrival_times <- sort(runif(n, 0, t))
  arrival_times <- append(arrival_times, 999999)  # add very large value for sake 
  # of ending the while loop
  # - Service times can be drawn from exp(mu) -- note not all times will be used
  #   due to balking
  service_times <- rexp(n, mu)
  
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
  while (ct < t){
    
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
      
      # Check to see if next event is an arrival and update accordingly
    }else if(is.na(tt_next_service) | tt_next_arrival < tt_next_service){
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
      
      # Check to see if next event is a service     
    }else if(tt_next_service < tt_next_arrival){
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
  
  # When fulloutput option is F (default), only return the number of customer
  # served. Otherwise, return event times/types and queue length vectors, along
  # with a figure illustrating queue lengths and balking over time. 
  if(fulloutput == F){
    return(k)  # default is to return just the total customers serviced
  }else{  # if full output option is chosen then...
    # Graph queue length over time
    # - Make data frame of queue length over time (add 0,0 as first row)
    qlength_df <- data.frame(time = c(0, times_of_events),
                             qlength = c(0, queue_length))
    # - Make data frame of balking and reneging times and choose y-coord that is 
    #   5% below 0 relative to the max queue length
    vax_times <- times_of_events[types_of_events == "served"]
    vax_df <- data.frame(vtimes = vax_times, y = -0.6*max(queue_length))
    if(sum(types_of_events == "balks") == 0){
      balk_df <- data.frame(btimes = as.numeric(NA), y = as.numeric(NA))
    }else{
      balk_times <- times_of_events[types_of_events == "balks"]
      balk_df <- data.frame(btimes = balk_times, y = -0.2*max(queue_length))
    }
    if(sum(types_of_events == "reneges") == 0){
      renege_df <- data.frame(rtimes = as.numeric(NA), y = as.numeric(NA))
    }else{
      renege_times <- times_of_events[types_of_events == "reneges"]
      renege_df <- data.frame(rtimes = renege_times, y = -0.4*max(queue_length))   
    }
    
    # - Make plot
    plot <- ggplot() +
      geom_step(data = qlength_df, mapping = aes(x = time, y = qlength), 
                lwd = 1.05) +
      geom_point(data = balk_df, aes(x = btimes, y = y), color = "#C00000", 
                 fill = "#C00000", alpha = 0.5, shape = 24, size = 6) +
      geom_point(data = renege_df, aes(x = rtimes, y = y), color = "#C00000", 
                 fill = "#C00000", shape = 23, alpha = 0.5, size = 6) +
      geom_point(data = vax_df, aes(x = vtimes, y = y), color = "skyblue3",
                 fill = "skyblue3", alpha = 0.5, shape = 21, size = 6) +
      xlab("Time, minutes") +
      ylab("Queue length") +
      scale_y_continuous(breaks = pretty_breaks(), 
                         limits = c(-0.7*max(queue_length), max(queue_length))) +
      theme_pubr() +
      theme(axis.text=element_text(size=20),
            axis.title=element_text(size=22, face = "bold")) 
    
    print(plot)
    print(table(types_of_events))
    
    # Return full output
    return(list(num_served = k, event_times = times_of_events, 
                event_types = types_of_events, queue_length = queue_length,
                plot = plot))
  }
}

# Create lookup table function for given values of mu, alpha, beta and T
CreateLookup <- function(mu, alpha, beta, t, min, max){
  num_arrivals <- seq(min, max, by = 1)
  num_served <- sapply(num_arrivals, function(x) EstimateCustomersServed(x, mu, alpha, beta, t))
  lu_table <- data.frame(num_arrivals = num_arrivals, num_served = num_served)
  return(lu_table)
}

# Function calculates the estimated number of customers served given the total
# number of arrivals (n), service rate (mu), balking/reneging parameters (α, β),
# and total operation time (t)
EstimateCustomersServed <- function(n, mu, alpha, beta, t){
  
  # Calculate lambda assuming the rate of arrivals is constant over time
  lambda <-n/t
  
  # Calculate vaccination rate 
  V <- CalculateV(lambda, alpha, beta, mu)
  
  # Calculate total customers served by multiplying vax rate by the total time
  return(V*t)
  
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


##-----------------------------------------------------------------------------
# 2. Run over multiple iterations to estimate trend in behavior
#    LOW balking
##-----------------------------------------------------------------------------

n_range <- seq(2, 200, by = 2)
length(n_range)  # 100

# First run 1000 iterations for each n value ----
# Note each iteration is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be 1000*4 reps per n value
n_iterations <- 1000
nrep_df <- data.frame(n = rep(n_range, each = n_iterations*4), 
                      day = rep(1:(4*n_iterations), length(n_range)),
                      iteration = rep(rep(1:n_iterations, each = 4), # 4 days per iteration
                                      length(n_range)),  
                      num_served = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*n_iterations)) {
    out_day <- SimulateBalkingReneging(n = n_i, mu = 0.5, alpha = 0.01, beta = 0.02, 
                                       t = 240)  # Simulate a single day of the campaign
    nrep_df$num_served[nrep_df$n == n_i & nrep_df$day == day] <- out_day
  }
}

# The above was run on a node in the HPC we'll save as RDS and reload locally
#saveRDS(nrep_df, "queueingsimulationreps_a0.01b0.02_10.6.23.rds")
nrep_df <- readRDS("output_hpc/queueingsimulationreps_a0.01b0.02_10.6.23.rds")

# Summarize number served per 4-day campaign
nrep_df2 <- nrep_df %>%
  group_by(n, iteration) %>%
  summarise(num_served = sum(num_served)) %>%
  mutate(num_arrivals = n*4) %>%
  select(iteration, num_arrivals, num_served)

# Make summary of output for each n value and convert for plotting
n_summary <- nrep_df2 %>%
  group_by(num_arrivals) %>%
  summarise(median = median(num_served),
            mean = mean(num_served),
            lower_25 = sort(num_served)[.25*n_iterations],
            upper_75 = sort(num_served)[.75*n_iterations])
theoretical <- CreateLookup(mu = 0.5, alpha = 0.01, beta = 0.02, t = 960,
                            min = min(n_range)*4, max = max(n_range)*4)
# Plot output
ggplot(n_summary, aes(x = num_arrivals)) +
  geom_line(data = theoretical, aes(x = num_arrivals, y = num_served), lwd = 1.1) +
  geom_line(aes(y = median), linetype = "dashed", color = "#C00000", lwd = 1.1) +
  geom_ribbon(aes(ymin = lower_25, ymax = upper_75), fill = "#C00000", alpha = 0.3) +
  xlab("Number of dogs arrived") +
  ylab("Number of dogs vaccinated") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 455)) +
  scale_x_continuous(limits = c(0, 600)) +
  theme_bw() +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24)) 

##-----------------------------------------------------------------------------
# 2. Run over multiple iterations to estimate trend in behavior
#    HIGH balking
##-----------------------------------------------------------------------------

n_range <- seq(2, 200, by = 2)
length(n_range)  # 100

# First run 1000 iterations for each n value ----
# Note each iteration is composed of 4 vax "days" lasting 4 hours. Thus, there 
# will be 1000*4 reps per n value
n_iterations <- 1000
nrep_df <- data.frame(n = rep(n_range, each = n_iterations*4), 
                      day = rep(1:(4*n_iterations), length(n_range)),
                      iteration = rep(rep(1:n_iterations, each = 4), # 4 days per iteration
                                      length(n_range)),  
                      num_served = NA)
for (i in 1:length(n_range)) {
  n_i <- n_range[i]
  print(paste0("Running simulations for n = ", n_i))
  for (day in 1:(4*n_iterations)) {
    out_day <- SimulateBalkingReneging(n = n_i, mu = 0.5, alpha = 0.1, beta = 0.1, 
                                       t = 240)  # Simulate a single day of the campaign
    nrep_df$num_served[nrep_df$n == n_i & nrep_df$day == day] <- out_day
  }
}

# The above was run on a node in the HPC we'll save as RDS and reload locally
#saveRDS(nrep_df, "output_hpc/queueingsimulationreps_a0.1b0.1_10.6.23.rds")
nrep_df_high <- readRDS("output_hpc/queueingsimulationreps_a0.1b0.1_10.6.23.rds")

# Summarize number served per 4-day campaign
nrep_df2_high <- nrep_df_high %>%
  group_by(n, iteration) %>%
  summarise(num_served = sum(num_served)) %>%
  mutate(num_arrivals = n*4) %>%
  select(iteration, num_arrivals, num_served)

# Make summary of output for each n value and convert for plotting
n_summary_high <- nrep_df2_high %>%
  group_by(num_arrivals) %>%
  summarise(median = median(num_served),
            mean = mean(num_served),
            lower_25 = sort(num_served)[.25*n_iterations],
            upper_75 = sort(num_served)[.75*n_iterations])
theoretical_high <- CreateLookup(mu = 0.5, alpha = 0.1, beta = 0.1, t = 960,
                            min = min(n_range)*4, max = max(n_range)*4)
# Plot output
ggplot(n_summary_high, aes(x = num_arrivals)) +
  geom_line(data = theoretical_high, aes(x = num_arrivals, y = num_served), lwd = 1.1) +
  geom_line(aes(y = median), linetype = "dashed", color = "#C00000", lwd = 1.1) +
  geom_ribbon(aes(ymin = lower_25, ymax = upper_75), fill = "#C00000", alpha = 0.2) +
  xlab("Number of dogs arrived") +
  ylab("Number of dogs vaccinated") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 455)) +
  scale_x_continuous(limits = c(0, 600)) +
  theme_bw() +
  #theme_pubr() +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24)) 
