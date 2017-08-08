# Defines functions: maxLikGlobalDimEst, maxLikPointwiseDimEst, maxLikLocalDimEst

################################################################################

maxLikGlobalDimEst <- function(data, k, dnoise = NULL, sigma = 0, n = NULL,
                    integral.approximation = 'Haro', unbiased = FALSE,
                    neighborhood.based = TRUE,
                    neighborhood.aggregation = 'maximum.likelihood',
                    iterations = 5, K = 5) {
  # 'k' is the number of neighbors used for each dimension estimation.
  # 'dnoise' is a vector valued function giving the transition density.
  # 'sigma' is the estimated standard deviation for the noise.
  # 'n' is the dimension of the noise (at least dim(data)[2])
  # integral.approximation can take values 'Haro', 'guaranteed.convergence', 'iteration'
  # neighborhood.based means that estimation is made for each neighborhood,
  # otherwise estimation is based on distances in entire data set.
  # 'K' is number of neighbors per data point that is considered, only used for 
  # neighborhood.based = FALSE

  N <- dim(data)[1]
    
  if (neighborhood.based) {
    mi <- maxLikPointwiseDimEst(data, k, dnoise, sigma, n,
                       integral.approximation = integral.approximation,
                       unbiased = unbiased)$dim.est
    de <- switch(neighborhood.aggregation,
                 maximum.likelihood = 1/mean(1/mi),
                 mean = mean(mi),
                 robust = median(mi)
           )
    return(DimEst(de))
  }

  dist <- yaIkNN(data, 1:N, K)[, (K+1):(2*K)]
  dist <- dist[!duplicated(dist[1:(K*N)])]
  Rs <- dist[rank(dist, ties.method = 'random') <= k]
  de <- maxLikDimEstFromR(Rs, dnoise, sqrt(2)*sigma, n, integral.approximation, unbiased,
            iterations) # Since distances between points are used, noise is 
                        # added at both ends, i.e. variance is doubled. 

  return(DimEst(de, likelihood = NA)) 
  
}


maxLikPointwiseDimEst <- function(data, k, dnoise = NULL, sigma = 0, n = NULL, indices = NULL,
                         integral.approximation = 'Haro', unbiased = FALSE,
                         iterations = 5) {
  ## estimates dimension around each point in data[indices, ]
  #
  # 'indices' give the indexes for which local dimension estimation should
  # be performed.
  # 'k' is the number of neighbors used for each local dimension estimation.
  # 'dnoise' is a vector valued function giving the transition density.
  # 'sigma' is the estimated standard deviation for the noise.
  # 'n' is the dimension of the noise (at least dim(data)[2])
  
  if (is.null(indices))
    indices <- 1:dim(data)[1]
  
  N <- length(indices)
  nbh.dist <- yaIkNN(data, indices, k)[ , (k + 1):(2*k)]
  de <- rep(NA, N) # This vector will hold local dimension estimates
  
  for (i in 1:N) {
    Rs <- nbh.dist[i, ]
    de[i] <- maxLikDimEstFromR(Rs, dnoise, sigma, n, integral.approximation, unbiased, iterations)
  }

  return(DimEstPointwise(de))
}


maxLikLocalDimEst <- function(data, dnoise = NULL, sigma = 0, n = NULL,
                   integral.approximation = 'Haro',
                   unbiased = FALSE, iterations = 5) {
  # assuming data set is local
  center <- apply(data, 2, mean)
  cent.data <- t(t(data) - center)
  Rs <- sort(lens(cent.data))
  de <- maxLikDimEstFromR(Rs, dnoise, sigma, n, integral.approximation, unbiased, iterations)
  return(DimEst(de))
}

################################################################################

maxLikDimEstFromR <- function(Rs, dnoise, sigma, n, integral.approximation = 'Haro',
                unbiased = FALSE, iterations = 5) {

  if (!(integral.approximation %in% c('Haro', 'guaranteed.convergence', 'iteration')))
    stop('Wrong integral approximation')


  if (!(typeof(dnoise) == 'closure') && !is.null(dnoise)) 
    dnoise <- match.fun(dnoise)
  dnoise_orig <- dnoise
  if (!(integral.approximation == 'Haro') && !is.null(dnoise))
    dnoise <- function(r, s, sigma, k) r*dnoise_orig(r, s, sigma, k)
  
  de <- maxLikDimEstFromR_haro_approx(Rs, dnoise, sigma, n, unbiased)
  if (integral.approximation == 'iteration')
    de <- maxLikDimEstFromRIterative(Rs, dnoise_orig, sigma, n, de, unbiased)

  return(de)
}

################################################################################

maxLikDimEstFromR_haro_approx <- function(Rs, dnoise, sigma, n, unbiased = FALSE) {
  # if dnoise is the noise function this is the approximation used in Haro.
  # for 'guaranteed.convergence' dnoise should be r times the noise function
  # with 'unbiased' option, estimator is unbiased if no noise or boundary

  k <- length(Rs)
  kfac <- if (unbiased) k-2 else k-1
  
  Rk <- max(Rs)
  if (is.null(dnoise)) return(kfac/(sum(log(Rk/Rs))))
  
  Rpr <- Rk + 100*sigma
  
  numerator <- rep(NA, k - 1)
  denominator <- rep(NA, k - 1)
  
  for (j in 1:(k-1)) {
    
    Rj <- Rs[j]
    tc <- tryCatch({
      numInt <- integrate(function(x) {
        dnoise(x, Rj, sigma, n) * log(Rk/x)
      }, 0, Rpr, rel.tol = 1e-2)
    }, error = function(ex) {
      print(ex)
      return(NA)
    }, finally = {})
    if (is.na(tc)[1]) return(NA)
    numerator[j] <- numInt$value
    
    tc <- tryCatch({
      denomInt <- integrate(function(x) {
        dnoise(x, Rj, sigma, n)
      }, 0, Rpr, rel.tol = 1e-2)
    }, error = function(ex) {
      print(ex)
      return(NA)
    }, finally = {})
    if (is.na(tc)[1]) return(NA)
    denominator[j] <- denomInt$value
  }

  return(kfac/sum(numerator/denominator)) 
  
}

################################################################################

maxLikDimEstFromRIterative <- function(Rs, dnoise, sigma, n, init = 5,
                  unbiased = FALSE, iterations = 5, verbose = FALSE) {
  m <- init
  if (verbose)
    cat("Start iteration, intial value:", m, "\n")
  for (i in 1:iterations) {
    m <- maxLikDimEstFromRIterative_inner(Rs, dnoise, sigma, n, m, unbiased)
    if (verbose)
      cat("Iteration", i, ":", m, "\n")
  }
  if (verbose)
    cat("\n")
  return(m)

}

maxLikDimEstFromRIterative_inner <- function(Rs, dnoise, sigma, n, m, unbiased) {

  k <- length(Rs)  
  kfac <- if (unbiased) k-2 else k-1
  
  Rk <- max(Rs)
  if (is.null(dnoise)) return(kfac/(sum(log(Rk/Rs))))
  Rpr <- Rk + 100*sigma

  numerator <- rep(NA, k - 1)
  denominator <- rep(NA, k - 1)

  for (j in 1:(k-1)) {
    Rj <- Rs[j]
    numInt <- integrate(function(x) {
                          m <- max(m, 1)
                          x^(m-1)*dnoise(x, Rj, sigma, n) * log(Rk/x)
                        }, 0, Rpr, rel.tol = 1e-2)
    numerator[j] <- numInt$value

    denomInt <- integrate(function(x) {
                            m <- max(m, 1)
                            x^(m-1)*dnoise(x, Rj, sigma, n)
                          }, 0, Rpr, rel.tol = 1e-2)
    denominator[j] <- denomInt$value
  }

  
  return(kfac/sum(numerator/denominator))

}

################################################################################

yaIkNN <- function(data, indexes, k) {
  # ann returns indexes and SQUARED Euclidean distances
  kNN.i.d <- ann(data, data[indexes, , drop = FALSE], k + 1, 
                    verbose = FALSE)$knnIndexDist[ , -c(1, k + 2), drop = FALSE]
  kNN.d <- sqrt(kNN.i.d[ , (k + 1):(2*k)])
  kNN.i.d[ , (k + 1):(2*k)] <- kNN.d
  return(kNN.i.d)
}

################################################################################

