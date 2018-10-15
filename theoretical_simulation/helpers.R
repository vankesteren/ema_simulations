# Filter Method
filterSel <- function(x, M, y, decisionFunction, ...) {
  meds <- names(M)
  xres <- x
  yres <- y
  
  msel <- numeric(length(meds))
  names(msel) <- meds
  
  for (med in meds) {
    # get mediator
    m <- as.matrix(M[med])
    
    if (decisionFunction(x, m, y, ...)) {
      msel[med] <- 1
    } else {
      msel[med] <- 0
    }
  }
  
  return(msel)
}

# RegSEM Method
xmedSel <- function(df, dir = FALSE, fitLav = TRUE, cv = TRUE, ...) {
  nmed <- ncol(df) - 2
  
  # create the model: x -> M -> y 
  lavMod <- "# Structural parameters\n"
  lavMod <- paste0(lavMod, paste(paste0("M.", 1:nmed, collapse = " + "), "~ x\n"))
  lavMod <- paste0(lavMod, paste("y ~", paste0("M.", 1:nmed, collapse = " + ")))
  # add direct effect
  if (dir) lavMod <- paste(lavMod, "+ x")
  
  lavMod <- paste(lavMod, "\n\n# Residual Covariances")
  for (i in 1:(nmed - 1)) {
    lavMod <- paste(lavMod, 
                    paste0("M.", i, " ~~ ", 
                           paste0("M.", (i + 1):nmed, collapse = " + ")),
                    sep = "\n")
  }
  
  # create lavaan object (optionally, do not fit)
  lav <- suppressWarnings(
    lavaan(lavMod, data = scale(df), auto.var = TRUE, do.fit = fitLav, 
           missing = "ml")
  )
  
  if (cv) {
    # fit with cross-validation BIC selection
    out <- capture.output(
      regsemfit <- try(cv_regsem(lav, pars_pen = 1:(2*nmed), 
                                 type = "lasso", ...),
                       silent = TRUE)
    )
    if (inherits(regsemfit, "try-error")) {
      warning("Regsem iteration on process ", Sys.getpid(), " failed.")
      return(rep(NA, nmed))
    }
    
    return(as.numeric(as.numeric(regsemfit$final_pars[1:nmed]) * 
             as.numeric(regsemfit$final_pars[(nmed + 1):(2*nmed)]) != 0))
    
  } else {
    # fit with a single lambda to be specified in the "..." argument
    regsemfit <- multi_optim(lav, pars_pen = 1:(2*nmed), alpha = 0, ...)
    
    return(as.numeric(as.numeric(regsemfit$out$pars[1:nmed]) *
                      as.numeric(regsemfit$out$pars[(nmed + 1):(2*nmed)]) != 0))
  }
}

# HIMA Method
himaSel <- function(df, p.value = 0.1) {
  selection <- numeric(ncol(df) - 2)
  names(selection) <- paste0("M.", 1:(ncol(df) - 2))
  himaRes <- hima(X = as.vector(df$x), 
                  Y = as.vector(df$y), 
                  M = as.matrix(df[,-c(1, ncol(df))]))
  selection[rownames(himaRes)] <- himaRes[["adjusted.p"]] < p.value
  return(selection)
}

# lavSel
lavSel <- function(df, p.value = 0.1, dir = FALSE, ...) {
  nmed <- ncol(df) - 2
  
  # create the model: x -> M -> y 
  lavMod <- "# Structural parameters\n"
  lavMod <- paste0(lavMod, paste(paste0("M.", 1:nmed, collapse = " + "), "~ x\n"))
  lavMod <- paste0(lavMod, paste("y ~", paste0("M.", 1:nmed, collapse = " + ")))
  # add direct effect
  if (dir) lavMod <- paste(lavMod, "+ x")
  
  lavMod <- paste(lavMod, "\n\n# Residual Covariances")
  for (i in 1:(nmed - 1)) {
    lavMod <- paste(lavMod, 
                    paste0("M.", i, " ~~ ", 
                           paste0("M.", (i + 1):nmed, collapse = " + ")),
                    sep = "\n")
  }
  
  
  lav <- suppressWarnings(
    lavaan(lavMod, data = scale(df), auto.var = TRUE, missing = "ml")
  )
  
  # sobel test
  stat <- lav@ParTable$est[1:nmed] * lav@ParTable$est[(nmed+1):(2*nmed)]
  se <- sqrt(lav@ParTable$est[1:nmed]^2 * lav@ParTable$se[(nmed+1):(2*nmed)]^2 +
               lav@ParTable$est[(nmed+1):(2*nmed)]^2 * lav@ParTable$se[1:nmed]^2)
  
  return(as.numeric(abs(stat/se) > qnorm(p.value/2, lower.tail = F)))
}
