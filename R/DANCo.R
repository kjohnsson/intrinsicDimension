#library(yaImpute)

dancoDimEst <- function(data, k, D, ver = 'DANCo', calibration.data = NULL) {
  cal <- calibration.data
  N <- dim(data)[1]
  if (!is.null(cal)) {
    if (cal$k != k)
      stop(sprintf("Neighborhood parameter k = %s does not agree with neighborhood parameter of calibration data, cal$k = %s",
                   k, cal$k))
    if (cal$N != N)
      stop(sprintf("Number of data points N = %s does not agree with number of data points of calibration data, cal$N = %s",
                   N, cal$N))
  }
  
  if (ver != 'DANCo') return(MIND_MLx(data, k, D, ver))
  
  nocal <- dancoDimEstNoCalibration(data, k, D, ver)
  if (any(is.na(nocal))) return(DimEst(NA, kl.divergence = NA, calibration.data = cal))

  if (is.null(cal)) {
    cal <- DancoCalibrationData(k, N)
  }

  if (cal$maxdim < D) {
    message(sprintf("Computing DANCo calibration data for N = %d, k = %d for dimensions %d to %d",
                    N, k, cal$maxdim+1, D))
  }
  while (cal$maxdim < D) cal <- increaseMaxDimByOne(cal)

  kl <- rep(NA, D)  
  for (d in 1:D) {
    kl[d] <- KL(nocal, cal$calibration.data[[d]], k) 
  }
  de <- which.min(kl)

  if(length(de) == 0) return(DimEst(NA, kl = NA, cal=cal))
  return(DimEst(de, kl.divergence = kl[de], calibration.data = cal))
}

################################################################################

MIND_MLx <- function(data, k, D, ver) {
  nbh.data <- yaIkNN(data, 1:dim(data)[1], k+1)
  rhos <- nbh.data[, k+2]/nbh.data[, 2*k+2]
  d.MIND_MLk <- MIND_MLk(rhos, k, D)
  if (ver == 'MIND_MLk') return(DimEst(d.MIND_MLk))
  d.MIND_MLi <- MIND_MLi(rhos, k, D, d.MIND_MLk)
  if (ver == 'MIND_MLi') return(DimEst(d.MIND_MLi))
  stop(sprintf("Unknown version: %s", ver))
}


dancoDimEstNoCalibration <- function(data, k, D, ver) {
  nbh.data <- yaIkNN(data, 1:dim(data)[1], k+1)
  rhos <- nbh.data[, k+2]/nbh.data[, 2*k+2]
  d.MIND_MLk <- MIND_MLk(rhos, k, D)
  d.MIND_MLi <- MIND_MLi(rhos, k, D, d.MIND_MLk)

  thetas <- angles(data, nbh.data[, 1:k])
  ml.vm <- apply(thetas, 1, ML_VM)
  mu_nu <- mean(ml.vm['nu', ])
  mu_tau <- mean(ml.vm['tau', ])
  if(dim(data)[2] == 1) mu_tau <- 1
  return(c(dhat = d.MIND_MLi, mu_nu = mu_nu, mu_tau = mu_tau))
}

KL <- function(nocal, caldat, k) {
  kld <- KLd(nocal['dhat'], caldat['dhat'], k)
  klnutau <- KLnutau(nocal['mu_nu'], caldat['mu_nu'],
                     nocal['mu_tau'], caldat['mu_tau'])
  #print(klnutau)
  return(c(kl = kld + klnutau))
}

KLd <- function(dhat, dcal, k) {
  H_k <- sum(1/(1:k))
  quo <- dcal/dhat
  a <- (-1)^(0:k)*choose(k, 0:k)*digamma(1 + (0:k)/quo)
  return(H_k*quo - log(quo) - (k-1)*sum(a))
}

KLnutau <- function(nu1, nu2, tau1, tau2) {
  return(log(besselI(tau2, 0)/besselI(tau1, 0)) + 
    besselI(tau1, 1)/besselI(tau1, 0)*
    (tau1 - tau2*cos(nu1-nu2)))
}

################################################################################

MIND_MLk <- function(rhos, k, D) {
  N <- length(rhos)
  d.lik <- rep(NA, D)
  for (d in 1:D) d.lik[d] <- lld(d, rhos, k, N)
  return(which.max(d.lik))
}

MIND_MLi <- function(rhos, k, D, dinit) {
  res <- optim(dinit, nlld, nlld.gr, rhos = rhos, k = k, N = length(rhos),
               method = 'L-BFGS-B', lower = 0, upper = D)
  #if(!is.null(res$message)) print(res$message)
  return(res$par)  
}

angles <- function(data, nbs) {
  N <- dim(data)[1]
  k <- dim(nbs)[2]
  thetas <- matrix(nrow = N, ncol = choose(k, 2))
  for (i in 1:N) {
    nb.data <- data[nbs[i, ], , drop = FALSE]
    thetas[i, ] <- loc.angles(data[i, ], nb.data)    
  }
  return(thetas)
}

ML_VM <- function(thetas) {
  sinth <- sin(thetas)
  costh <- cos(thetas)
  nu <- atan(sum(sinth)/sum(costh))
  eta <- sqrt(mean(costh)^2 + mean(sinth)^2)
  tau <- Ainv(eta)
  return(c(nu = nu, tau = tau))
}

################################################################################

nlld <- function(d, rhos, k, N) return(-lld(d, rhos, k, N))

lld <- function(d, rhos, k, N) {
  if (d == 0) return(-1e30)
  N*log(k*d) + (d-1)*sum(log(rhos)) + (k-1)*sum(log(1-rhos^d))
}

nlld.gr <- function(d, rhos, k, N) {
  if (d == 0) return(-1e30)
  -(N/d + sum(log(rhos) - (k-1)*rhos^d*log(rhos)/(1 - rhos^d)))
}

Ainv <- function(eta) {
  if (eta < .53) return(2*eta + eta^3 + 5*eta^5/6)
  if (eta < .85) return(-.4 + 1.39*eta + .43/(1-eta))
  return(1/(eta^3-4*eta^2+3*eta))
}

loc.angles <- function(pt, nbs) {
  vec <- t(apply(nbs, 1, function(nb) { nb - pt }))
  if (length(pt) == 1) vec <- t(vec)
  vec.len <- lens(vec)
  combs <- combn(dim(nbs)[1], 2)
  sc.prod <- apply(vec[combs[1, ], , drop = FALSE]*vec[combs[2, ], , drop = FALSE], 1, sum)
  #if (length(pt) == 1) {
    #print(sc.prod)
    #print((vec.len[combs[1, ]]*vec.len[combs[2, ]]))
  #}
  cos.th <- sc.prod/(vec.len[combs[1, ]]*vec.len[combs[2, ]])
  if (any(abs(cos.th) > 1)) print(cos.th[abs(cos.th) > 1])
  return(acos(cos.th))
}

