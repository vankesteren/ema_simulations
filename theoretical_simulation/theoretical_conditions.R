# Theoretical Simulations supporting the manuscript
# "Exploratory Mediation Analysis with Many Potential Mediators"
# Last edited: 11/09/2018
# (c) 2017 - 2018 Erik-Jan van Kesteren
# This work was supported by NWO Talent Grant 406.17.057

# For outline in RStudio, press ctrl+shift+o / cmd+shift+o

# Initialisation ----
# install packages
# source("tools/install.R")

# filter functions
library(cmfilter)
library(regsem)
library(HIMA)
source("helpers.R")

# simulation libraries
library(MASS)
library(parallel)
library(pbapply)

# Sim params
nIter <- 100
nCores <- detectCores() - 1
sampleSizes <- c(400, 400, 400, 600)
outFolder <- "Paper_Results_20180807"

if (!dir.exists(file.path("logs", outFolder))) 
  dir.create(file.path("logs", outFolder))
if (!dir.exists(file.path("output", outFolder))) 
  dir.create(file.path("output", outFolder))

# Simulation 1: Suppression ----

# The data-generating model
S <- matrix(c( 1, -.6, -.6, 1), 2)
apaths <- c(-0.4, 0.4)
bpaths <- c(0.4, 0.24)
rsquared <- bpaths %*% S %*% bpaths

# Set up the cluster
set.seed(142857)
clus <- makeCluster(nCores, outfile = file.path("logs", outFolder, "sup.log"))
clusterEvalQ(clus, {
  library(cmfilter)
  library(regsem)
  library(HIMA)
  source("helpers.R")
})
clusterExport(clus, c("S", "apaths", "bpaths", "rsquared", "sampleSizes"))

sup_olist <- pblapply(1:nIter, function(x) {
  pid <- Sys.getpid()
  set.seed(pid + x)
  seed <- round(runif(1)*100000)
  message("PID: ", pid, " | Iter: ", x, " | seed: ", seed)
  set.seed(seed)
  
  # generate data
  d <- generateMed(n = sampleSizes[1], 
                   a = apaths,
                   b = bpaths,
                   Sigma = S,
                   r2y = rsquared,
                   empirical = TRUE)
  
  # proper Sobel test with multiple mediator model
  message(pid, " | Sobel") 
  sobel <- lavSel(df = d, p.value = 0.1, dir = FALSE)
  
  # filter
  message(pid, " | filter_prodcoef") 
  filter_prodcoef <- filterSel(d$x, d[, -c(1,ncol(d))], d$y, 
                               decisionFunction = cmfilter:::prodCoef,
                               p.value = 0.1,
                               dir = TRUE)
  
  message(pid, " | filter_difcoef") 
  filter_difcoef <- filterSel(d$x, d[, -c(1,ncol(d))], d$y, 
                              decisionFunction = cmfilter:::corMinusPartCor,
                              p.value = 0.1)
  
  # xmed
  message(pid, " | xmed") 
  xmed_res <- xmedSel(df = d, dir = FALSE, fitLav = TRUE, cv = TRUE)
  
  # hima
  message(pid, " | hima") 
  hima_res <- himaSel(df = d, p.value = 0.1)
  
  # cmf
  message(pid, " | cmf_pc_res") 
  cmf_pc_res <- cmf(d, decisionFunction = "prodcoef", nCores = 1,
                    nStarts = 1000, p.value = 0.1, pb = FALSE)
  
  message(pid, " | cmf_dc_res") 
  cmf_dc_res <- cmf(d, decisionFunction = "causalsteps", nCores = 1,
                    nStarts = 1000, p.value = 0.1, pb = FALSE)
  
  return(
    list(
      sobel = sobel,
      filter_prodcoef = filter_prodcoef, 
      filter_difcoef = filter_difcoef,
      xmed = xmed_res,
      hima = hima_res, 
      cmf_prodcoef = cmf_pc_res$selection,
      cmf_difcoef = cmf_dc_res$selection
    )
  )
}, cl = clus)

stopCluster(clus)

# Save the output list
save(sup_olist, file = file.path("output", outFolder, "sup_olist.Rdata"))


# Simulation 2: Noise a ----

# The data-generating model

# generate a random covariance matrix with small but nonzero covariance
# which will serve as the residual covariance of M 
# https://stats.stackexchange.com/a/215647/116878
# to tune the rate parameter, run source("tools/tune_rate.R")
set.seed(142857)
p <- 16
P <- qr.Q(qr(matrix(rnorm(p^2), p))) # eigenvectors
rate <- 1.1
e <- (rate^(p:1)/rate*p)/sum(rate^(p:1)/rate) # eigenvalues sum to p
S <- cov2cor(crossprod(P, P * e))

apaths <- c(0.3, sign(S)[-1,1]*c(rep(0.8, 3), rep(0.4, 12)))
bpaths <- c(0.3, rep(0, 15))

Sigma <- diag(1 - apaths^2)
S <- S * tcrossprod(diag(Sigma))
diag(S) <- 0
Sigma <- Sigma + S

rsquared <- 0.5 # bpaths %*% Sigma %*% bpaths

# Set up the cluster

clus <- makeCluster(nCores, outfile = file.path("logs", outFolder, "noia.log"))
clusterEvalQ(clus, {
  library(cmfilter)
  library(regsem)
  library(HIMA)
  source("helpers.R")
})
clusterExport(clus, c("Sigma", "apaths", "bpaths", "rsquared", "sampleSizes"))


noia_olist <- pblapply(1:nIter, function(x) {
  pid <- Sys.getpid()
  set.seed(pid + x)
  seed <- round(runif(1)*100000)
  message("PID: ", pid, " | Iter: ", x, " | seed: ", seed)
  set.seed(seed)
  
  # generate data
  d <- generateMed(n = sampleSizes[2], 
                   a = apaths,
                   b = bpaths,
                   Sigma = Sigma,
                   residual = TRUE,
                   r2y = rsquared,
                   empirical = TRUE)
  
  # proper Sobel test with multiple mediator model
  message(pid, " | Sobel") 
  sobel <- lavSel(df = d, p.value = 0.1, dir = FALSE)
  
  # filter 
  message(pid, " | filter_prodcoef") 
  filter_prodcoef <- filterSel(d$x, d[, -c(1,ncol(d))], d$y, 
                               decisionFunction = cmfilter:::prodCoef,
                               p.value = 0.1,
                               dir = TRUE)
  
  message(pid, " | filter_difcoef") 
  filter_difcoef <- filterSel(d$x, d[, -c(1,ncol(d))], d$y, 
                              decisionFunction = cmfilter:::corMinusPartCor,
                              p.value = 0.1)
  
  # xmed
  message(pid, " | xmed") 
  xmed_res <- xmedSel(df = d, dir = FALSE, fitLav = TRUE, cv = TRUE)
  
  # hima
  message(pid, " | hima") 
  hima_res <- himaSel(df = d, p.value = 0.1)
  
  # cmf
  message(pid, " | cmf_pc_res") 
  cmf_pc_res <- cmf(d, decisionFunction = "prodcoef", nCores = 1,
                    nStarts = 1000, p.value = 0.1, 
                    pb = FALSE)
  
  message(pid, " | cmf_dc_res") 
  cmf_dc_res <- cmf(d, decisionFunction = "causalsteps", nCores = 1,
                    nStarts = 1000, p.value = 0.1, 
                    pb = FALSE)
  
  return(
    list(
      sobel = sobel,
      filter_prodcoef = filter_prodcoef, 
      filter_difcoef = filter_difcoef,
      xmed = xmed_res,
      hima = hima_res, 
      cmf_prodcoef = cmf_pc_res$selection,
      cmf_difcoef = cmf_dc_res$selection
    )
  )
}, cl = clus)

stopCluster(clus)

# Save the output list
save(noia_olist, file = file.path("output", outFolder, "noia_olist.Rdata"))


# Simulation 3: Noise b ----

# The data-generating model

# generate a random covariance matrix with small but nonzero covariance
# which will serve as the residual covariance of M 
# https://stats.stackexchange.com/a/215647/116878
# to tune the rate parameter, run source("tools/tune_rate.R")
set.seed(142857)
p <- 16
P <- qr.Q(qr(matrix(rnorm(p^2), p))) # eigenvectors
rate <- 1.1
e <- (rate^(p:1)/rate*p)/sum(rate^(p:1)/rate) # eigenvalues sum to p
S <- cov2cor(crossprod(P, P * e))

apaths <- c(0.3, rep(0, 15))
bpaths <- c(0.3, sign(S)[-1,1]*c(rep(0.8, 3), rep(0.4, 12)))

Sigma <- diag(1-apaths^2)
S <- S * tcrossprod(diag(Sigma))
diag(S) <- 0
Sigma <- Sigma + S

rsquared <- 0.5

# Set up the cluster

clus <- makeCluster(nCores, outfile = file.path("logs", outFolder, "noib.log"))
clusterEvalQ(clus, {
  library(cmfilter)
  library(regsem)
  library(HIMA)
  source("helpers.R")
})
clusterExport(clus, c("Sigma", "apaths", "bpaths", "rsquared", "sampleSizes"))


noib_olist <- pblapply(1:nIter, function(x) {
  pid <- Sys.getpid()
  set.seed(pid + x)
  seed <- round(runif(1)*100000)
  message("PID: ", pid, " | Iter: ", x, " | seed: ", seed)
  set.seed(seed)
  
  # generate data
  d <- generateMed(n = sampleSizes[3], 
                   a = apaths,
                   b = bpaths,
                   Sigma = Sigma,
                   residual = TRUE,
                   r2y = rsquared,
                   empirical = TRUE)
  
  # proper Sobel test with multiple mediator model
  message(pid, " | Sobel") 
  sobel <- lavSel(df = d, p.value = 0.1, dir = FALSE)
  
  
  message(pid, " | filter_prodcoef") 
  filter_prodcoef <- filterSel(d$x, d[, -c(1,ncol(d))], d$y, 
                               decisionFunction = cmfilter:::prodCoef,
                               p.value = 0.1,
                               dir = TRUE)
  message(pid, " | filter_difcoef")
  filter_difcoef <- filterSel(d$x, d[, -c(1,ncol(d))], d$y, 
                              decisionFunction = cmfilter:::corMinusPartCor,
                              p.value = 0.1)
  
  # xmed
  message(pid, " | xmed_res")
  xmed_res <- xmedSel(df = d, dir = FALSE, fitLav = TRUE, cv = TRUE)
  
  # hima
  message(pid, " | hima_res")
  hima_res <- himaSel(df = d, p.value = 0.1)
  
  # cmf
  message(pid, " | cmf_pc_res")
  cmf_pc_res <- cmf(d, decisionFunction = "prodcoef", nCores = 1,
                    nStarts = 1000, p.value = 0.1, 
                    pb = FALSE)
  
  message(pid, " | cmf_dc_res")
  cmf_dc_res <- cmf(d, decisionFunction = "causalsteps", nCores = 1,
                    nStarts = 1000, p.value = 0.1, 
                    pb = FALSE)
  
  return(
    list(
      sobel = sobel,
      filter_prodcoef = filter_prodcoef, 
      filter_difcoef = filter_difcoef,
      xmed = xmed_res,
      hima = hima_res, 
      cmf_prodcoef = cmf_pc_res$selection,
      cmf_difcoef = cmf_dc_res$selection
    )
  )
}, cl = clus)

stopCluster(clus)

# Save the output list
save(noib_olist, file = file.path("output", outFolder, "noib_olist.Rdata"))

# Simulation 4: Noise + Suppression ----
# The data-generating model

# generate a random covariance matrix with small but nonzero covariance
# which will serve as the residual covariance of M 
# https://stats.stackexchange.com/a/215647/116878
# to tune the rate parameter, run source("tools/tune_rate.R")
set.seed(142857)
p <- 15
rate <- 1.1
e <- (rate^(p:1)/rate*p)/sum(rate^(p:1)/rate) # eigenvalues sum to p

Pa <- qr.Q(qr(matrix(rnorm(p^2), p))) # eigenvectors of noise a
Sa <- cov2cor(crossprod(Pa, Pa * e)) # covmat of noise a

Pb <- qr.Q(qr(matrix(rnorm(p^2), p))) # eigenvectors of noise b
Sb <- cov2cor(crossprod(Pb, Pb * e)) # covmat of noise b

Sm <- matrix(c(1,-.44, -.44, 1))

S <- diag(32)

S[1:2,1:2] <- Sm
S[3:17,3:17] <- Sa
S[18:32,18:32] <- Sb
diag(S) <- 1-apaths^2

apaths <- c(-0.4, 0.4, rep(0.8, 3), rep(0.4, 12), rep(0, 15))
bpaths <- c(0.4, 0.24, rep(0, 15), rep(0.8, 3), rep(0.4, 12))
rsquared <- 0.5

# Set up the cluster

clus <- makeCluster(nCores, outfile = file.path("logs", outFolder, "nosu.log"))
clusterEvalQ(clus, {
  library(cmfilter)
  library(regsem)
  library(HIMA)
  source("helpers.R")
})
clusterExport(clus, c("S", "apaths", "bpaths", "rsquared", "sampleSizes"))

nosu_olist <- pblapply(1:nIter, function(x) {
  pid <- Sys.getpid()
  set.seed(pid + x)
  seed <- round(runif(1)*100000)
  message("PID: ", pid, " | Iter: ", x, " | seed: ", seed)
  set.seed(seed)
  
  # generate data
  d <- generateMed(n = sampleSizes[4], 
                   a = apaths,
                   b = bpaths,
                   Sigma = S,
                   r2y = rsquared,
                   empirical = TRUE,
                   residual = TRUE)
  
  
  # proper Sobel test with multiple mediator model
  message(pid, " | Sobel") 
  sobel <- lavSel(df = d, p.value = 0.1, dir = FALSE)
  
  # filter
  message(pid, " | filter_prodcoef") 
  filter_prodcoef <- filterSel(d$x, d[, -c(1,ncol(d))], d$y, 
                               decisionFunction = cmfilter:::prodCoef,
                               p.value = 0.1,
                               dir = TRUE)
  
  # xmed
  message(pid, " | xmed") 
  xmed_res <- xmedSel(df = d, dir = FALSE, fitLav = TRUE, cv = TRUE)
  
  # hima
  message(pid, " | hima") 
  hima_res <- himaSel(df = d, p.value = 0.1)
  
  # cmf
  message(pid, " | cmf_pc_res") 
  cmf_pc_res <- cmf(d, decisionFunction = "prodcoef", nCores = 1, p.value = 0.2,
                    nStarts = 1500, pb = FALSE)

  
  return(
    list(
      sobel = sobel,
      filter_prodcoef = filter_prodcoef, 
      xmed = xmed_res,
      hima = hima_res, 
      cmf_prodcoef = cmf_pc_res$selection
    )
  )
}, cl = clus)

stopCluster(clus)

# Save the output list
save(nosu_olist, file = file.path("output", outFolder, "nosu_olist.Rdata"))


# Done. See analysis.R for the analysis of the results.

