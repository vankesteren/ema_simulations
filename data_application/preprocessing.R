# Generate data frame from participant info
tab <- read.table("data/samples.tsv", header = TRUE, sep = "\t")


# x and y variables
x <- tab$Characteristics..total.score.on.the.childhood.trauma.questionnaire.
y <- tab$FactorValue..cortisol.stress.response.area.under.the.curve.auc.with.respect.to.the.increase.

# Covariates
Z <- cbind(rep(1, nrow(tab)), 
           tab$Characteristics..age., 
           tab$Characteristics..sex.)
H <- Z %*% solve(crossprod(Z)) %*% t(Z) # hat matrix


# Residualise x and y
x <- x - H %*% x
y <- y - H %*% y

# M variables
indic <- sapply(strsplit(as.character(tab$Source.Name), " "), function(x) x[1])
filename <- paste0("data/", indic[1], "_sample_table.txt")
mdat <- read.delim(filename)
M <- matrix(0.0, nrow(tab), nrow(mdat))
colnames(M) <- mdat$Reporter.Identifier

M[1,] <- mdat$VALUE
for (i in 2:length(indic)) {
  filename <- paste0("data/", indic[i], "_sample_table.txt")
  M[i,] <- read.delim(filename)$VALUE
}

# Transform to M value; see Du et al., 2010
M <- log2(M / (1 - M))

# residualize
M <- M - H %*% M

# preselect / marginal filtering the top 1000
omegaAB <- abs(as.numeric(crossprod(M, y)) * as.numeric(crossprod(M, x)))
idxAB <- order(omegaAB, decreasing = TRUE)[1:1000]
MselAB <- M[, idxAB]

# save result
save(x, y, MselAB, file = "results/preproc.Rdata")
