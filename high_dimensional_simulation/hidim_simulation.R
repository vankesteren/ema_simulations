# High-dimensional simulation supporting the manuscript
# "Exploratory Mediation Analysis with Many Potential Mediators"
# Last edited: 14/05/2018
# (c) 2017 - 2018 Erik-Jan van Kesteren
# This work was supported by NWO Talent Grant 406.17.057

# For outline in RStudio, press ctrl+shift+o

# Load packages ----
library(tidyverse)
library(cmfilter)
library(regsem)
library(HIMA)
library(glmnet)
library(MASS)
library(Matrix)
library(parallel)
library(pbapply)

# Load other R files ----
source("hidim_simulation-helper.R")

# Parameters ----
pars <- list(
  noiseProp = 1,                      # Ratio of noise-true mediators
  p.value = 0.1,                      # p-value to operate at
  rSquared = 0.5,                     # Proportion of variance explained in y
  N = 100L,                           # Sample size
  P = 1000L,                          # Number of variables
  Ptrue = 10L,                        # Number of true mediators
  noisePaths = 0.7,                   # Noise in the a and b paths
  paths = 0.5,                        # Size of the alpha and beta paths
  dir = 0                             # Presence of a direct effect
)

# Run the simulation ----
clus <- makeCluster(detectCores(), outfile = file.path("logs", "perf.log"))
pcks <- clusterEvalQ(clus, {
  library(cmfilter)
  library(regsem)
  library(HIMA)
  library(glmnet)
  library(MASS)
  library(Matrix)
  library(parallel)
  source("hidim_simulation-helper.R")
})

result <- pblapply(1:500, FUN = perfSim, pars = pars, cl = clus)

stopCluster(clus)

# Save the results ----
save(result, file = file.path("output", "high_dimensional_results.Rdata"))

# Set cutoff to 0.075 ----
for (i in seq_along(result)) 
  result[[i]]$result$cmf <- setCutoff(result[[i]]$result$cmf, .075)
save(result, file = file.path("output", "high_dimensional_results_cutoff.Rdata"))
