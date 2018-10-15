load("results/preproc.Rdata")

# apply CMF
library(cmfilter)
cmfres <- cmf(x, MselAB, y, nStarts = 100)

for (i in 1:100) {
  cmfres <- update(cmfres)
  # For convergence animation create an ani folder and uncomment the lines below
  # if (i < 10) {
  #   fname <- paste0("ani/00", i, ".png")
  # } else if (i < 100) {
  #   fname <- paste0("ani/0", i, ".png")
  # } else if (i < 1000) {
  #   fname <- paste0("ani/", i, ".png")
  # }
  # png(fname, width = 1920, height = 1080, res = 150)
  # plot(cmfres, xaxt = "n", ylim = c(0,0.15), main = i, border = NA)
  # dev.off()
}

cmfres <- setCutoff(cmfres, 0.075)

save(cmfres, file = "results/CMFresult.Rdata")


# apply HIMA for comparison
library(HIMA)
himres <- hima(x, y, MselAB)
himres <- himres[order(himres$BH.FDR), ]

save(himres, file = "results/HIMAresult.Rdata")


# and the filter method
filres <- sapply(1:ncol(MselAB), function(i) {
  m     <- MselAB[,i]
  n     <- length(x)
  cpx   <- crossprod(x)
  alpha <- solve(cpx, crossprod(x, m))
  res_m <- m - x * c(alpha)
  var_m <- as.numeric(crossprod(res_m)/(n - 1))
  var_a <- var_m/cpx
  mm    <- cbind(m)
  cpm   <- crossprod(mm)
  beta  <- solve(cpm, crossprod(mm, y))
  res_y <- y - mm %*% c(beta)
  var_y <- as.numeric(crossprod(res_y)/(n - 1))
  var_b <- diag(var_y * chol2inv(chol(cpm)))
  stat  <- alpha * beta
  se    <- sqrt(alpha^2 * var_b + beta^2 * var_a)
  return(abs(stat/se))
})

names(filres) <- colnames(MselAB)
filres <- filres[order(filres, decreasing = TRUE)]
save(filres, file = "results/FilterResult.Rdata")
