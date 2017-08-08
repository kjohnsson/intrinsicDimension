lens <- function(vectors) 
  sqrt(apply(vectors^2, 1, sum))

rowSubtract <- function(matrix, row)
  t(t(matrix) - row)

pcaOtpmPointwiseDimEst <- function(data, N, alpha = .05) {
  # N is the number of nodes in the graph

  km <- kmeans(data, N)
  pt <- km$centers
  pt_bm <- km$cluster
  pt_sm <- rep(NA, length(pt_bm))

  for (k in 1:dim(data)[1]) {
    pt_sm[k] <- which.min(lens(rowSubtract(pt[-pt_bm[k],  ], data[k, ])))
    if (pt_sm[k] >= pt_bm[k]) pt_sm[k] <- pt_sm[k] + 1
  }

  de.c <- rep(NA, N)
  nbr.nb.c <- rep(NA,N)
  for (k in 1:N) {
    nb <- unique(c(pt_sm[pt_bm == k], pt_bm[pt_sm == k]))
    nbr.nb.c[k] <- length(nb)
    loc.dat <- rowSubtract(pt[nb, ], pt[k, ])
    de.c[k] <- pcaLocalDimEst(loc.dat, ver = 'FO', alphaFO = alpha)$dim.est
  }
  de <- de.c[pt_bm]
  nbr.nb <- nbr.nb.c[pt_bm]
  return(DimEstPointwise(de, list("nbr.neighbors" = nbr.nb)))

}