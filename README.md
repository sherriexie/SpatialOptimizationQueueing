# Incorporating Queueing Theory into a Spatial Optimization Framework to Improve Mass Vaccination Programs

This repository contains scripts used in the analysis presented in the manuscript titled _Incorporating Queueing Theory into a Spatial Optimization Framework to Improve Mass Vaccination Programs_. A description of the repository's contents are provided below.


## Optimization script files

Scripts used to optimize the locations of vaccination sites via the hybrid recursive interchange-genetic algorithm are located in the **scripts_optimization** folder. A description of these scripts are provided below:
* **functions_recursiveinterchange.R** - contains all the functions needed to run the recursive interchange component of the hybrid algorithm
* **functions_geneticalgorithm.R** - contains all the functions needed to run the genetic algorithm component of the hybrid algorithm
* **optimize...R** - each script file runs a single round of the hybrid algorithm
  - The values of alpha and beta are specified in the portion of the script name that comes after the underscore "_" (alpha = beta = 0 are specified as "queuenaive")
  - The round of the optimization is given by the number following optimize
  - So as an example, optimize2_alpha0.1beta0.1.R was used to run the second round of optimization for the high-attrition scenario (alpha = beta = 0.1)
* **check_stoppingcondition.R** - determines at the conclusion of each round of optimization whether the stopping condition had been met


## Other script files

Scripts used to run the queueing simulations, generate figures, and conduct the sensitivity analyses are located in the **scripts_other** folder. A description of these scripts are provided below:
* **simulate_vsequilibriumsolutions.R** - compares results of the stochastic queueing model to the closed-form equations
* **figure_simulations.R** - runs single trials of the queueing model and plots the results
* **figure_potentialvaxsites.R** - maps the locations of all potential vaccination sites
* **figure_arrivalhistograms.R** - generates the arrivals histograms that were presented in the manuscript
* **figure_catchmentmaps.R** - maps the locations of optimized vaccination sites and numbers of expected arrivals that were presented in the manuscript
* **sfig_poissonregression.R** - fits and plots the MDVC participation probability function and the post-MDVC survey data that was used in the regression
* **sa_balkingreneging.R** - evaluates sensitivity of results to misspecification of alpha and beta
* **sa_inconstantarrivals....R** - evaluates sensitivity of results to the assumption of a constant arrival rate

## Data files

Non-sensitive data used in the analysis are provided as a single RDA file, which can be found in the **data** folder. Note that the locations of houses in the study area and the probability-distance matrix (containing the probability of each household participating in the campaign for all possible vaccination sites) have been removed to protect privacy.
* **optimization_data.rda**
  - asa: shapefile of Alto Selva Algre (study area)
  - queuedata: times of arrivals collected during an MDVC
  - supplycoord: coordinates of the potential vaccination sites in UTM projection
  - survey: survey data used to fit the MDVC participation probability function
  
