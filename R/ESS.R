################################################################################
# Functions:
#
#   computeEss
#   essLocalDimEst
#   essReference
#
################################################################################

essLocalDimEst <- function(data, ver = 'a', d = 1) {
  essval <- computeEss(data, verbose = FALSE, ver, d)
  if (is.na(essval)) return(c(de = NA, ess = NA))
  mindim <- 1
  maxdim <- 20
  dimvals <- essReference(ver, d, maxdim, mindim)
  while ((ver == 'a' && essval > dimvals[maxdim]) ||
      (ver == 'b' && essval < dimvals[maxdim])) {
    mindim <- maxdim + 1
    maxdim <- 2*(maxdim-1)
    dimvals <- c(dimvals, essReference(ver, d, maxdim, mindim))
  }
  if (ver == 'a') {
    i <- findInterval(essval, dimvals[mindim:maxdim])
  } else {
    i <- length(maxdim:mindim) - findInterval(essval, dimvals[maxdim:mindim])
  }
  de_integer <- mindim+i-1
  de_fractional <- (essval-dimvals[de_integer])/(dimvals[de_integer+1]-dimvals[de_integer])
  de <- de_integer + de_fractional
  return(DimEst(de, ess = essval))
}

################################################################################

computeEss <- function(data, verbose = FALSE, ver, d = 1) {

  p <- d + 1

  n <- dim(data)[2]
  if (p > n) {
    if (ver == 'a') return(0)
    if (ver == 'b') return(1)
    stop('Not a valid version')
  }
  
  vectors <- vecToCs_onedir(data, 1)  
  if (verbose) {
    message('Number of vectors:', dim(vectors)[1], '\n')
  }
  groups <- indnComb(dim(vectors)[1], p)
  if (dim(groups)[1] > 5000) groups <- groups[sample(1:dim(groups)[1], 5000), ]
  if (verbose) {
    message('Number of simple elements:', dim(groups)[1], '\n')
  }
  Allist <- apply(groups, 1, function(row) { list(vectors[row, ]) })
  Alist <- lapply(Allist, function(el) { el[[1]] })
  
  # Compute weights for each simple element
  weight <- sapply(Alist, function(vecgr) { prod(lens(vecgr)) } )
  if (ver == 'a') {
    # Compute the volumes of the simple elements
    vol <- sapply(Alist, function(vecgr) { sqrt(det(vecgr %*% t(vecgr))) } )
    return(sum(vol)/sum(weight))
  }
  if (ver == 'b') {
    if (d == 1) {
      # Compute the projection of one vector onto one other
      proj <- sapply(Alist, function(vecgr) { abs(sum(vecgr[1, ] * vecgr[2, ])) } )
      return(sum(proj)/sum(weight))
    }
    stop('For ver == "b", d > 1 is not supported.')
  }
  stop('Not a valid version')
}

################################################################################

essReference <- function(ver, d, maxdim, mindim=1) {
  
  if (maxdim <= d + 2) {
    stop(paste0("maxdim (", maxdim, ") must be larger than d + 2 (", d + 2, ")"))
  }
  if (ver == 'a') {
    ## dim.val(n) = factor1(n)^d * factor2(n)
    # factor1(n) = gamma(n/2)/gamma((n+1)/2)
    # factor2(n) = gamma(n/2)/gamma((n-d)/2)
    
    ## compute factor1
    # factor1(n) = gamma(n/2)/gamma((n+1)/2)
    # [using the rule gamma(n+1) = n * gamma(n)] repeatedly
    # = gamma(1/2)/gamma(2/2) * prod{j \in J1} j/(j+1) if n is odd
    # = gamma(2/2)/gamma(3/2) * prod(j \in J2) j/(j+1) if n is even
    # where J1 = seq(1, n-2, 2), J2 = seq(2, n-2, 2)
    J1 <- seq(1, maxdim, 2)
    J2 <- seq(2, maxdim, 2)
    factor1.J1 <- gamma(1/2)/gamma(2/2) * c(1, cumprod(J1/(J1+1))[-length(J1)])
    factor1.J2 <- gamma(2/2)/gamma(3/2) * c(1, cumprod(J2/(J2+1))[-length(J2)])
    factor1 <- rep(NA, maxdim)
    factor1[J1] <- factor1.J1
    factor1[J2] <- factor1.J2
    
    ## compute factor2
    # factor2(n) = gamma(n/2)/gamma((n-d)/2)
    # = gamma((d+1)/2)/gamma(1/2) * prod{k \in K1} k/(k-d) if n-d is odd
    # = gamma((d+2)/2)/gamma(2/2) * prod(k \in K2) k/(k-d) if n-d is even
    # where K1 = seq(d+1, n-2, 2), K2 = seq(d+2, n-2, 2)
    # if n > d+2, otherwise 0.
    K1 <- seq(d+1, maxdim, 2)
    K2 <- seq(d+2, maxdim, 2)
    factor2.K1 <- gamma((d+1)/2)/gamma(1/2) * c(1, cumprod(K1/(K1-d))[-length(K1)])
    factor2.K2 <- gamma((d+2)/2)/gamma(2/2) * c(1, cumprod(K2/(K2-d))[-length(K2)])
    factor2 <- rep(0, maxdim)
    factor2[K1] <- factor2.K1
    factor2[K2] <- factor2.K2
    # compute dim.val
    dim.val <- factor1^d * factor2
    dim.val <- dim.val[mindim:maxdim]
    return(dim.val)
  } 
  if (ver == 'b') {
    if (d == 1) {
      # dim.val(n) = 2*pi^(-1/2)/n *gamma((n+1)/2)/gamma((n+2)/2)
      # = gamma(2/2)/gamma(3/2) * prod{j \in J1} (j+1)/(j+2) * 2/sqrt(pi)/n if n is odd
      # = gamma(3/2)/gamma(4/2) * prod(j \in J2) (j+1)/(j+2) * 2/sqrt(pi)/n if n is even
      # where J1 = seq(1, n-2, 2), J2 = seq(2, n-2, 2)
      J1 <- seq(1, maxdim, 2)
      J2 <- seq(2, maxdim, 2)
      dim.val.J1 <- gamma(3/2)/gamma(2/2) * c(1, cumprod((J1+2)/(J1+1))[-length(J1)])
      dim.val.J2 <- gamma(4/2)/gamma(3/2) * c(1, cumprod((J2+2)/(J2+1))[-length(J2)])
      dim.val <- rep(NA, maxdim)
      dim.val[J1] <- dim.val.J1
      dim.val[J2] <- dim.val.J2
      # n <- mindim:maxdim
      # return(gamma((n+2)/2)/gamma((n+1)/2) * 2/sqrt(pi)/n)
      return(dim.val[mindim:maxdim] * 2/sqrt(pi)/(mindim:maxdim))
    }
    stop('For ver == "b", d > 1 is not supported.')
  }
  stop('Not a valid version')
  
}

################################################################################