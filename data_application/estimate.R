library(lavaan)
load("results/annotation.Rdata")
load("results/preproc.Rdata")

d <- data.frame(x = x, MselAB[, rownames(annot)], y = y)

# Generate a model in lavaan syntax
lavmod <- ""
for (i in 1:nrow(annot)) {
  lavmod <- paste0(lavmod, "\n", rownames(annot)[i], " ~ a", i, "*x")
}
lavmod <- paste0(lavmod, "\n\ny ~ ")
for (i in 1:nrow(annot)) {
  if (i != 1)  lavmod <- paste0(lavmod, " + ")
  lavmod <- paste0(lavmod, "b", i, "*", rownames(annot)[i])
}
for (i in 1:nrow(annot)) {
  lavmod <- paste0(lavmod, "\n\ni", i, " := a", i, "*b", i, "\n")
}
for (i in 1:(nrow(annot) - 1)) {
  lavmod <- paste0(lavmod, "\n", rownames(annot)[i], " ~~ ", 
                   paste0(rownames(annot)[-c(1:i)], collapse = " + "))
}

# Estimate the lavaan model with bootstrap standard errors
lav <- sem(lavmod, scale(d), se = "boot")

# Save the results to a text file
sink("results/lavaan_summary.txt")
summary(lav)
cat("\n\n\nAnnotation -----------\n\n")
annot
sink()
