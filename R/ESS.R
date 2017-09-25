################################################################################
# Functions:
#
#   computeEss
#   essLocalDimEst
#   essReference
#   vnball
#   vnsphere
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
  
  if (ver == 'a') {
    if (d == 1) {
      n <- max(mindim, 2):maxdim
      dim.val <- gamma(n/2)^2/(gamma((n-1)/2)*gamma((n+1)/2))
    } else {
      n <- max(mindim, (d+1)):maxdim
      dim.val <- (gamma(n/2)/gamma((n+1)/2))^d*gamma(n/2)/gamma((n-d)/2)
    }
    if (mindim <= d) {
      dim.val <- c(rep(0, d-mindim+1), dim.val)
    }
    return(dim.val)
  } 
  if (ver == 'b') {
    if (d == 1) {
      n <- mindim:maxdim
      return(2*vnball(n-1)/vnsphere(n-1))
    }
    stop('For ver == "b", d > 1 is not supported.')
  }
  stop('Not a valid version')

}

################################################################################

vnball <- function(n) return(pi^(n/2)/gamma(n/2+1))
vnsphere <- function(n) return((n+1)*vnball(n+1))

################################################################################
