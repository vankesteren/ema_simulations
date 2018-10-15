library(FDb.InfiniumMethylation.hg18)
library(cmfilter)

load("results/CMFresult.Rdata")

# annotate (get nearest gene)
hm450 <- get450k()
cpg <- names(cmfres$selection[cmfres$selection])
probes <- hm450[cpg]
annot <- getNearestTSS(probes)

annot$description <- c(
  "Involved in transcriptional regulation",
  "Associated with the age-at-onset of diabetes",
  "Associated with sexual development",
  "Involved in regulation of cell growth and survival",
  "Involved in regulation of mRNA"
)

save(annot, file = "results/annotation.Rdata")
