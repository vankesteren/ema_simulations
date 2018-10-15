# Illustrative Simulations supporting the manuscript
# "Exploratory Mediation Analysis with Many Potential Mediators"
# Last edited: 27/06/2018
# (c) 2017 - 2018 Erik-Jan van Kesteren
# This work was supported by NWO Talent Grant 406.17.057

# preamble ----
# use tidyverse
library(tidyverse)

subfolder <- "Paper_Results_20180807"

# load the data
load(file.path("output", subfolder, "sup_olist.Rdata" ))
load(file.path("output", subfolder, "noia_olist.Rdata"))
load(file.path("output", subfolder, "noib_olist.Rdata"))
load(file.path("output", subfolder, "nosu_olist.Rdata"))

# utility functions
toTibble <- function(resultList, name) {
  # takes an output list and a method name and outputs a tibble
  # where rows indicate iterations and columns indicate mediators
  map(resultList, name) %>%
    sapply(identity) %>%
    t %>%
    as.tibble %>%
    set_names(paste0("M.", 1:length(resultList[[1]][[1]])))
}


# analysis ----
lapply(names(sup_olist[[1]]),  function(name) sup_olist  %>% toTibble(name)) %>%
  set_names(names(sup_olist[[1]])) %>% 
  map(function(x) colSums(x, na.rm = TRUE)) -> sup

lapply(names(noia_olist[[1]]), function(name) noia_olist %>% toTibble(name)) %>%
  set_names(names(noia_olist[[1]])) %>% 
  map(function(x) colSums(x, na.rm = TRUE)) -> noia

lapply(names(noib_olist[[1]]), function(name) noib_olist %>% toTibble(name)) %>%
  set_names(names(noib_olist[[1]])) %>% 
  map(function(x) colSums(x, na.rm = TRUE)) -> noib

lapply(names(nosu_olist[[1]]), function(name) nosu_olist %>% toTibble(name)) %>%
  set_names(names(nosu_olist[[1]])) %>% 
  map(function(x) colSums(x, na.rm = TRUE)) -> nosu

nIter <- length(sup_olist)

suppression <- 
  cbind(M1 = sup %>% map_dbl(function(x) x[1]/nIter), 
        M2 = sup %>% map_dbl(function(x) x[2]/nIter))

noise_a <- 
  cbind(Power = noia %>% map_dbl(function(x) x[1]/nIter), 
        FPR = noia %>% map_dbl(function(x) mean(x[-1])/nIter),
        PPV = noia %>% map_dbl(function(x) x[1]/sum(x)))

noise_b <- 
  cbind(Power = noib %>% map_dbl(function(x) x[1]/nIter), 
        FPR = noib %>% map_dbl(function(x) mean(x[-1])/nIter),
        PPV = noib %>% map_dbl(function(x) x[1]/sum(x)))

noise_and_suppression <-
  cbind(M1 = nosu %>% map_dbl(function(x) x[1]/nIter),
        M2 = nosu %>% map_dbl(function(x) x[2]/nIter),
        FPR = nosu %>% map_dbl(function(x) mean(x[-c(1,2)])/nIter),
        PPV = nosu %>% map_dbl(function(x) sum(x[1:2])/sum(x)))


# results ----

meths <- c("sobel", "filter_prodcoef", "xmed", "hima", "cmf_prodcoef")

# Suppression
suppression %>% 
  as.tibble(rownames = "met") %>% 
  filter(met %in% meths) %>% 
  mutate(met = c("SEM", "Filter", "XMed", "HIMA", "CMF")) %>% 
  set_names(c("Method", "$M_1$", "$M_2$")) ->
  supptable

save(supptable, file = "output/tables/suppression.Rdata")
  
# noise A
t(as.data.frame(noia)) %>% 
  as.tibble(rownames = "met") %>% 
  filter(met %in% meths) %>% 
  mutate(met = c("SEM", "Filter", "XMed", "HIMA", "CMF")) %>% 
  set_names(c("Method", "$M$", 2:16)) %>% 
  mutate_all(funs(replace(., . == 0, "."))) ->
  noiatable

save(noiatable, file = "output/tables/noise_a.Rdata")

# noise B
t(as.data.frame(noib)) %>% 
  as.tibble(rownames = "met") %>% 
  filter(met %in% meths) %>% 
  mutate(met = c("SEM", "Filter", "XMed", "HIMA", "CMF")) %>% 
  set_names(c("Method", "$M$", 2:16)) %>% 
  mutate_all(funs(replace(., . == 0, "."))) ->
  noibtable

save(noibtable, file = "output/tables/noise_b.Rdata")

# both
t(as.data.frame(nosu)) 

save(noistable, file = "output/tables/noise_sup.Rdata")

noise_and_suppression %>% 
  as.tibble(rownames = "met") %>% 
  mutate(met = c("SEM", "Filter", "XMed", "HIMA", "CMF")) %>% 
  set_names(c("Method", "Power $M_1$", "Power $M_2$", "FPR", "PPV")) ->
  noistable


save(noistable, file = "output/tables/noise_sup.Rdata")
