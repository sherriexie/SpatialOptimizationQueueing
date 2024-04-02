# Incorporating Queueing Theory into a Spatial Optimization Framework to Improve Mass Vaccination Programs

This repository contains scripts used in the analysis presented in the manuscript titled _Incorporating Queueing Theory into a Spatial Optimization Framework to Improve Mass Vaccination Programs_. A description of the repository's contents are provided below.

## Optimization script files

Scripts used to optimize the locations of vaccination sites via the hybrid recursive interchange-genetic algorithm are located in the *scripts_optimization* folder. A description of these scripts are provided below:
* **functions_recursiveinterchange.R** - contains all the functions needed to run the recursive interchange component of the hybrid algorithm
* **functions_geneticalgorithm.R** - contains all the functions needed to run the genetic algorithm component of the hybrid algorithm
* **optimize...R** - each script file runs a single round of the hybrid algorithm
  - The values of alpha and beta are specified in the portion of the script name that comes after the underscore "_" (alpha = beta = 0 are specified as "queuenaive")
  - The round of the optimization is given by the number following optimize
  - So as an example, optimize2_alpha0.1beta0.1.R was used to run the second round of optimization for the high-attrition scenario (alpha = beta = 0.1)
* **check_stoppingcondition.R** -  
