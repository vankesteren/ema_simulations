if (!require(manipulate)) {install.packages("manipulate"); library(manipulate)}

manipulate({
    corImg <- function(cor, ...) {
      p <- dim(cor)[1]
      if (sum(diag(cor)) != p) cor <- cov2cor(cor)
      plotS <- abs(cor)
      plotS[lower.tri(cor)] <- NA
      diag(plotS) <- NA
      
      plot(NULL, xlim = c(0,1), ylim = c(1,0), asp = 1, axes = FALSE, 
           ylab = "", ...)
      image(plotS, add = T, 
            col = colorRampPalette(c("darkblue", "#FF8300"))(200), 
            breaks = seq(0, 1, len = 201))
    }
    
    par(mfrow = c(2, 2))
    set.seed(`Eigenvector seed`)
    p <- `Number of variables`
    P <- qr.Q(qr(matrix(rnorm(p^2), p))) # eigenvectors
    rate <- exp(`log(Eigenvalue decay rate)`)
    
    # Eigenvalues
    e <- rate^(p:1)/rate
    e <- e*p/sum(e) # eigenvalues sum to p
    if (fixaxis) elim <- c(0, p) else elim <- c(0, max(e) + 1)
    plot(e, type = "b", xlab = "Eigenvalue #", ylab = "value", 
         main = paste0("Eigenvalues (Rate = ", round(rate, 3), ")"),
         ylim = elim)
    
    # Create cor
    S <- cov2cor(crossprod(P, P * e))
    Sinv <- crossprod(P, P * 1/e)
    PS <- Sinv/sqrt(tcrossprod(diag(Sinv))) # partial cormat
    
    # Density plot
    d <- density(abs(S[lower.tri(S, diag = FALSE)]))
    d$y <- d$y/max(d$y)
    plot(d, main = "Normalised density of correlations", 
         xlim = c(0, 1), ylim = c(0,1))
    
    # Heatmap of cor
    means <- round(mean(abs(S[lower.tri(S)])), 3)
    corImg(S, xlab = paste("Mean abs =", means), 
           main = "Correlations")
    
    # Heatmap of pcor
    meanps <- round(mean(abs(PS[lower.tri(PS)])), 3)
    corImg(PS, xlab = paste("Mean abs =", meanps), 
           main = "Partial correlations")
    par(mfrow = c(1, 1))
  }, 
  `log(Eigenvalue decay rate)` = slider(0.0, 4.9), 
  fixaxis = checkbox(label = "Fix eigenvalue y-axis"),
  `Number of variables` = slider(2, 64, initial = 16), 
  `Eigenvector seed` = slider(1, 25))
