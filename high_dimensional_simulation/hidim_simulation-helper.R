# Helper functions for the main simulation
# R Script supporting the manuscript
# "Exploratory Mediation Analysis with Many Potential Mediators"
# Last edited: 14/05/2018
# (c) 2017 - 2018 Erik-Jan van Kesteren
# This work was supported by NWO Talent Grant 406.17.057

# input: simulation parameters
#     vpars = varying simulation parameters
#     fpars = fixed simulation parameters
#     cpars = control parameters

# output: simulation result

perfSim <- function(i, pars) {
  # CONTROL STUFF ----
  pars$pid <- Sys.getpid()
  if (is.null(pars$seed)) {
    pars$seed <- pars$pid + as.numeric(Sys.time())
  } 
  set.seed(pars$seed)
  
  # GENERATE DATASET ----
  # X -> M and M -> Y paths
  pnoise <- ceiling(pars$Ptrue * pars$noiseProp) # noise vars per block
  
  
  # X -> ( M , A )
  apaths <- sparseVector(x = c(rep(pars$paths, pars$Ptrue), 
                               rep(pars$noisePaths, pnoise)),
                         i = 1:(pars$Ptrue + pnoise),
                         len = pars$P)
  
  # ( M , B ) -> Y
  bpaths <- sparseVector(x = c(rep(pars$paths, pars$Ptrue), 
                               rep(pars$noisePaths, pnoise)),
                         i = c(1:pars$Ptrue, 
                               (pars$Ptrue+pnoise+1):(pars$Ptrue+2*pnoise)),
                         len = pars$P)
  
  # small nonzero residual covariance matrix for M block
  SresM <- generateCovMat(pars$Ptrue, 1.05)
  diag(SresM) <- 1 - apaths[1:pars$Ptrue]^2
  
  # Small nonzero residual correlations for A & B blocks
  SresA <- generateCovMat(pnoise)
  SresB <- generateCovMat(pnoise)
  diag(SresA) <- rep(1 - pars$noisePaths^2, pnoise)
  
  # I block
  pI <- pars$P - pars$Ptrue - 2 * pnoise
  I <- sparseMatrix(1:pI, 1:pI, x = 1)
  
  Sres <- sparseBlockMatrix(SresM, SresA, SresB, I)
  
  d <- generateMed(
           n = pars$N, 
           a = apaths,
           b = bpaths,
         dir = pars$dir,
         r2y = pars$rSquared,
       Sigma = Sres,
    residual = TRUE
  )
  
  # CREATE RESULTS ----
  result <- applyMethods(d, pars)
  
  return(list(pars = pars, result = result))
}


generateCovMat <- function(p, rate = 1.1) {
  # source("tools/tune_rate.R")
  P <- qr.Q(qr(matrix(rnorm(p^2), p))) # eigenvectors
  e <- (rate^(p:1)/rate*p)/sum(rate^(p:1)/rate) # eigenvalues sum to p
  return(cov2cor(crossprod(P, P * e)))
}


sparseBlockMatrix <- function(...) {
  argg <- c(as.list(environment()), list(...))
  stopifnot(all(vapply(argg, function(x) ncol(x) == nrow(x), TRUE)))
  p <- sum(vapply(argg, ncol, 1))
  out <- Matrix::sparseMatrix(p, p, x = 1.0)
  idx <- 1
  for (i in 1:length(argg)) {
    block <- argg[[i]]
    blockp <- ncol(block)
    out[idx:(idx+blockp-1), idx:(idx+blockp-1)] <- block
    idx <- idx + blockp
  }
  out
}


# Function to apply methods to a generated dataset
applyMethods <- function(d, pars) {
  
  out <- list()
  
  cat("\n", pars$pid, ": Filter")
  out$filter <- filterSel(
    x = d, 
    p.value = pars$p.value, 
    dir = TRUE
  )
  
  cat("\n", pars$pid, ": CMF")
  out$cmf <- cmf(
    x = d, 
    maxIter = 25, 
    stableLag = 5,
    nStarts = 100, 
    nCores = 1, 
    progressBar = FALSE,
    p.value = pars$p.value
  )
  
  cat("\n", pars$pid, ": HIMA")
  out$hima <- suppressMessages(hima(X = d$x, Y = d$y, M = d[,-c(1,ncol(d))]))
  
  cat("\n", pars$pid, ": Naive Lasso")
  out$nlasso <- nlasso(x = d)
  
  
  return(out)
}


# Filter Method
filterSel <- function(x, p.value = 0.05, dir = TRUE) {
  M <- x[,-c(1,ncol(x))]
  y <- x$y
  x <- x$x
  meds <- names(M)
  
  msel <- numeric(length(meds))
  names(msel) <- meds
  
  n <- length(x)                                  # number of samples
  cpx <- c(crossprod(x))                          # cross product of X
  cpxM <- apply(M, 2, crossprod, y = x)           # cross product of x and M
  alpha <- cpxM/cpx                               # alpha paths
  res_M <- M - tcrossprod(x, alpha)               # residual of m~x+0
  var_M <- apply(res_M, 2, crossprod) / (n - 1)   # rss variance
  var_a <- var_M / cpx                            # variance of alpha
  
  beta <- var_b <- numeric(ncol(M))
  for (i in 1:ncol(M)) {
    m <- M[,i]
    # then the beta path
    if (dir) {
      mm <- cbind(x, m)                                 # model matrix
    } else {
      mm <- cbind(m)
    }
    cpm <- crossprod(mm)                                # cross product of mm
    b <- solve(cpm, crossprod(mm, y))                   # beta
    res_y <- y - mm %*% c(b)                            # residual of y~m+x+0
    var_y <- as.numeric(crossprod(res_y) / (n - 1))     # rss variance
    vb <- diag(var_y * chol2inv(chol(cpm)))             # variance of beta
    beta[i] <- b[2]
    var_b[i] <- vb[2]
  }
  
  stat <- alpha * beta # product of coefficients
  se <- sqrt(alpha^2 * var_b + beta^2 * var_a) # - var_a * var_b
  
  return(stat/se)
}


# Helper for easily generating lavaan model syntax with many variables
expandLavMod <- function(model.syntax = "", verbose = FALSE) {
  revLine <- function(line) {
    sapply(lapply(strsplit(line, NULL), rev), paste, collapse = "")
  }
  if (length(model.syntax) == 0) {
    stop("lavaan ERROR: empty model syntax")
  }
  model.syntax <- gsub("[#!].*(?=\\n)", "", model.syntax, perl = TRUE)
  model.syntax <- gsub(";", "\\n", model.syntax, fixed = TRUE)
  model.syntax <- gsub("[ \\t]+", "", model.syntax, perl = TRUE)
  model.syntax <- gsub("\\n{2,}", "\\n", model.syntax, perl = TRUE)
  model <- unlist(strsplit(model.syntax, "\\n"))
  ops <- "[(\\=\\~)(\\<\\~)(\\~\\*\\~)(\\~\\~)(\\=\\=)(\\:\\=)\\:\\~\\<\\>\\|\\%]"
  expModel <- list()
  for (i in 1:length(model)) {
    line <- model[i]
    if (grepl("%-%", line, fixed = TRUE)) {
      # if expander before operator
      expBefOp <- paste0("\\%\\-\\%.*(?=", ops, ")")
      if (grepl(expBefOp, line, perl = TRUE)) {
        # remove all after operator
        toExpand <- revLine(sub(paste0("^.*?", ops), "", revLine(line), perl = TRUE))
        
        # create n functions
        fromTo <- unlist(strsplit(toExpand, "%-%", fixed = TRUE))
        numbers <- sub(".*?(?=[0-9])", "", fromTo, perl = TRUE)
        name <- sub("(?=[0-9]).*", "", fromTo, perl = TRUE)
        if (name[1] != name[2]) {
          stop("lavaan ERROR: expander failed (varnames not equal)")
        }
        expanded <- lapply(paste0(name[1], numbers[1]:numbers[2]), function(x) sub(toExpand, x, line))
        
      } else {
        # remove all before operator
        toExpand <- sub(paste0("^.*?", ops), "", line, perl = TRUE)
        
        
        # create one additive function
        fromTo <- unlist(strsplit(toExpand, "%-%", fixed = TRUE))
        numbers <- sub(".*?(?=[0-9])", "", fromTo, perl = TRUE)
        name <- sub("(?=[0-9]).*", "", fromTo, perl = TRUE)
        if (name[1] != name[2]) {
          stop("lavaan ERROR: expander failed (varnames not equal)")
        }
        after <- paste(paste0(name[1], numbers[1]:numbers[2]), collapse = " + ")
        
        expanded <- sub(toExpand, after, line)
      }
      
      expModel <- c(expModel, expanded)
      
    } else {
      expModel <- c(expModel, line)
    }
  }
  
  out <- paste(expModel, collapse = ";")
  if (verbose) cat(gsub(";","\n",out), "\n")
  return(invisible(out))
}


# RegSEM Method
xmed <- function(x, dir = FALSE, fitLav = TRUE, cv = TRUE, ...) {
  nmed <- ncol(x) - 2
  
  # create the model: x -> M -> y with expander function loaded before
  lavMod <- paste0("M.1 %-% M.", nmed, " ~ x
                   y ~ M.1 %-% M.", nmed)
  expMod <- expandLavMod(lavMod, FALSE)
  
  # add direct effect
  if (dir) expmod <- paste(expmod, "+ x")
  
  # create lavaan object (optionally, do not fit)
  lav <- suppressWarnings(
    lavaan(expMod, data = scale(x), auto.var = TRUE, do.fit = fitLav, 
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
                        as.numeric(regsemfit$final_pars[(nmed + 1):(2*nmed)])))
    
  } else {
    # fit with a single lambda to be specified in the "..." argument
    regsemfit <- multi_optim(lav, pars_pen = 1:(2*nmed), alpha = 0, ...)
    
    return(as.numeric(as.numeric(regsemfit$out$pars[1:nmed]) *
                        as.numeric(regsemfit$out$pars[(nmed + 1):(2*nmed)])))
  }
}


# Naive lasso
nlasso <- function(x) {
  # inverse regression for the first part
  xm <- glmnet(as.matrix(x[,-c(1, 1002)]), x$x, 
               lambda = cv.glmnet(as.matrix(x[,-c(1, 1002)]), x$x)$lambda.1se)
  
  # normal lasso for the second part
  my <- glmnet(as.matrix(x[,-c(1, 1002)]), x$y, 
               lambda = cv.glmnet(as.matrix(x[,-c(1, 1002)]), x$y)$lambda.1se)
  
  # return the multiplied result
  return((coef(xm) * coef(my))[-1])
}
