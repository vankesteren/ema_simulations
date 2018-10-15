# create plot
load("CMFresult.Rdata")
library(cmfilter)
library(tidyverse)
library(ggrepel)
library(firatheme)

ggdat <- data.frame(selrate = cmfres$selectionRate, 
                    name    = names(cmfres$selectionRate),
                    lab     = ifelse(cmfres$selection, names(cmfres$selectionRate), NA))

gg <- ggplot(ggdat, aes(y = selrate, x = name, label = lab)) +
  geom_bar(stat = "identity", width = 1, fill = firaCols[1]) + 
  geom_abline(slope = 0, intercept = cmfres$call$cutoff, lty = 3, size = 0.3) +
  geom_label_repel(nudge_y = 0.005, label.r = unit(0, "lines"), 
                   label.padding = unit(0.4, "lines"), segment.alpha = 0, 
                   force = 0.1) + 
  ylim(c(0, 0.12)) +
  labs(y = "Empirical Selection Rate",
       x = "Methylation Probe") +
  theme_fira() +
  theme(axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(), 
        panel.grid.major = element_blank())

gg

firaSave("img/selrate.pdf", width = 12, height = 6)
