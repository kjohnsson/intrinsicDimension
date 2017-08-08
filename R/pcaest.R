# PCA dimension estimation methods given eigenvalues of covariance matrix

################################################################################

pcaLocalDimEst <- function(data, ver, alphaFO = .05, alphaFan = 10, betaFan = .8,
                 PFan = .95, ngap = 5, maxdim = min(dim(data)),
                 verbose = TRUE) {
  lambda <- prcomp(data)$sdev^2
  if (ver == 'FO') return(FO(lambda, alphaFO))
  if (ver == 'fan') return(fan(lambda, alphaFan, betaFan, PFan))
  if (ver == 'maxgap') return(maxgap(lambda))
  if (ver == 'cal') return(cal(lambda, ngap, dim(data)[1], maxdim, verbose))
  stop('Not a valid version.')
}

################################################################################

FO <- function(lambda, alpha) {
  de <- sum(lambda > alpha*lambda[1])
  n <- length(lambda)
  gaps <- lambda[1:(n-1)]/lambda[2:n]
  return(DimEst(de, gap.size = gaps[de]))
}

################################################################################  

fan <- function(lambda, alpha = 10, beta = .8, P = .95) {
  
  n <- length(lambda)
  r <- which(cumsum(lambda)/sum(lambda) > P)[1]
  sigma <- mean(lambda[r:n])
  lambda <- lambda - sigma
  gaps <- lambda[1:(n-1)]/lambda[2:n]
  de <- min(c(which(gaps > alpha), 
               which(cumsum(lambda)/sum(lambda) > beta)))
  return(DimEst(de, gap.size = gaps[de]))
  
}

################################################################################

maxgap <- function(lambda) {
  n <- length(lambda)
  gaps <- lambda[1:(n-1)]/lambda[2:n]
  de <- which.max(gaps)
  return(DimEst(de, gap.size = gaps[de]))
}

################################################################################

cal <- function(lambda, ngap = 5, Ns, maxdim = min(length(lambda), Ns), verbose = TRUE) {
  n <- length(lambda)
  if (Ns < n) {
    n <- Ns
    lambda <- lambda[1:n]
  }
  gaps <- lambda[1:(n-1)]/lambda[2:n]
  D <- order(gaps[1:maxdim], decreasing = TRUE)[1:ngap]
  lik <- rep(NA,ngap)
  for (j in 1:length(D)) {
    d <- D[j]
    sigma.noi <- mean(lambda[(d+1):n])
    sigma.data <- mean(lambda[1:d])
    sigma.ball <- 1/(d+2)
    R <- sqrt(sigma.data/sigma.ball)
    sd.noi <- sqrt(sigma.noi)/R
    ntest <- 100
    if (verbose) {
      cat('Computing likelihood for d =', d, '\n')
      cat('R =', R, '\n')
      cat('sd.noi =', sd.noi, '\n')
    }
    tryCatch({
      lambdamat <- replicate(ntest, prcomp(cutHyperPlane(Ns, d, n, sd.noi))$sdev^2)
      gap <- lambdamat[d,]/lambdamat[d+1,]
      last <- lambdamat[n,]
      lik[j] <- dt((gaps[d] - mean(gap))/sd(gap), ntest-1)*
        dt((lambda[n] - mean(last))/sd(last), ntest-1)
      if (verbose) cat('Likelihood:', lik[j], '\n')
    }, error = function(ex) {
      cat(d, 'removed as possible dimension since reference data could not be constructed \n')
    })
  }
  ind = which.max(lik)
  return(DimEst(D[ind], likelihood = lik(ind)))  
}

################################################################################
