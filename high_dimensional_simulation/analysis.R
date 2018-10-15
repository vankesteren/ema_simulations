# Analysis of the results from the high-dimensional simulation
# R Script supporting the manuscript
# "Exploratory Mediation Analysis with Many Potential Mediators"
# Last edited: 11/09/2018
# (c) 2017 - 2018 Erik-Jan van Kesteren
# This work was supported by NWO Talent Grant 406.17.057

library(tidyverse)
library(firatheme) # devtools::install_github("vankesteren/firatheme")

load("output/high_dimensional_results_cutoff.Rdata")

# CMF
cmfpow <- sapply(result, function(res) {
  sum(res$result$cmf$selection[1:10])/10
})

cmferr <- sapply(result, function(res) {
  sum(res$result$cmf$selection[-c(1:10)])/990
})

cmfppv <- sapply(result, function(res) {
  sum(res$result$cmf$selection[1:10])/sum(res$result$cmf$selection)
})

cmfres <- c(mean(cmfpow), mean(cmferr), mean(cmfppv, na.rm = TRUE))


# filter
filpow <- sapply(result, function(res) {
  sum(abs(res$result$filter[1:10]) > qnorm(0.95))/10
})

filerr <- sapply(result, function(res) {
  sum(abs(res$result$filter[-c(1:10)]) > qnorm(0.95))/990
})

filppv <- sapply(result, function(res) {
  sum(abs(res$result$filter[1:10]) > qnorm(0.95))/sum(abs(res$result$filter) > qnorm(0.95))
})

filres <- c(mean(filpow), mean(filerr), mean(filppv, na.rm = TRUE))


# HIMA
himpow <- sapply(result, function(res) {
  sum(paste0("M.",1:10) %in% rownames(res$result$hima)) / 10
})

himerr <- sapply(result, function(res) {
  sum(!paste0("M.", 1:10) %in% rownames(res$result$hima)) / 990
})

himppv <- sapply(result, function(res) {
  sum(paste0("M.",1:10) %in% rownames(res$result$hima)) / nrow(res$result$hima)
})

himres <- c(mean(himpow), mean(himerr), mean(himppv, na.rm = TRUE))


# Naive lasso
nlapow <- sapply(result, function(res) {
  sum(res$result$nlasso[1:10] != 0) / 10
})

nlaerr <- sapply(result, function(res) {
  sum(res$result$nlasso[-c(1:10)] != 0) / 990
})

nlappv <- sapply(result, function(res) {
  sum(res$result$nlasso[1:10] != 0) / sum(res$result$nlasso != 0)
})

nlares <- c(mean(nlapow), mean(nlaerr), mean(nlappv, na.rm = TRUE))

hidim_tab <- rbind(CMF           = cmfres, 
                   Filter        = filres, 
                   HIMA          = himres, 
                   `Naive Lasso` = nlares)
colnames(hidim_tab) <-  c("Power", "Type I Error", "PPV")
save(hidim_tab, file = "output/tables/hidim_tab.Rdata")

