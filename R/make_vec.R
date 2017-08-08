################################################################################
# Defines the functions:
#
# len <- function(vector)
#   Computes the Euclidean length of 'vector'.
#
# orth <- function(A)
#   Constructs orthogonal matrix from input matrix via Gram-Schmidt on columns.
#
# normProj <- function(vec1, vec2)
#   Computes the projection of vec1/|vec1| onto vec2
#
# rowSubtract <- function(matrix, row)
#   Subtract 'row' from each row in 'matrix'. 'row' should be a vector
#
# indComb <- function(NN)
#   Returns all possible (unique) choices of two indexes smaller 
#   than or equal to NN.
#
# indnComb <- function(NN, n)
#   Returns all possible (unique) choices of n indexes smaller than
#   or equal to NN.
#
# addOpposite <- function(vecOneDir)
#   Add vectors in opposite direction so that vectors are paired as described
#   for the function vectorsBetween.
#
# vectorsBetween <- function(points)
#   Returns all vectors between the points in 'points'. The vectors are paired
#   so that the second vector is the opposite of the first, the fourth is the
#   opposite of the third, and so on.
#
# vectorsToCenter <- function(points, add.mids, weight.mids = 1,
#                             mids.max.dist = Inf))
#   Returns all vectors between points in 'points' and the mass center of the
#   points. The vectors are paired in the same way as for vectorsBetween.
#   Adds vectors to and from midpoints between two points closer to each other
#   than mids.max.dist if add.mids = TRUE.
#
# vectorsToCenter2 <- function(points, weight.mids = 1,
#                              max.proj = 1))
#   Returns all vectors between points in 'points' and the mass center of the
#   points. The vectors are paired in the same way as for vectorsBetween.
#   If add.mids = TRUE, it adds vectors to and from midpoints between any two 
#   points for which the vectors to the points have smaller normalized 
#   projection on each other than max.proj.
#
# vectorsToCenter3 <- function(points, weight.quarts = 1, max.proj = 1)
#   Returns all vectors between points in 'points' and the mass center of the
#   points. For each pair of points for which vectors to them (from the center)
#   have smaller normalized projection on each other than max.proj, three pairs
#   of vectors are added, going to three points evenly spaces between the two
#   orignial points.
#
# vectorsToCenters <- function(points, n.group)
#   For each possible group with n.group points the mass center is computed and 
#   the vectors to the mass center from the points in the group, as well as the
#   the vectors from the mass center to the points in ghe group, are returned
#   for each group.
#
# condense <- function(ind, values, obsolete)
#   The entries in the matrix 'ind' that also are in 'obsolete' and 
#   the corresponding entries in the matrix 'values' are taken away
#   and the remaining entries are shifted rowwise (NA's are introduced in
#   the end of the rows where needed).
#
# depth <- function(index, data)
#   Depth of a point in a data set as defined in Carter et al (2010).
#
################################################################################

len <- function(vector) {
  sqrt(sum(vector^2))
}

################################################################################

lens <- function(vectors) {
  sqrt(apply(vectors^2, 1, sum))
}

################################################################################

normProj <- function(vec1, vec2) {
  sum(vec1 * vec2)/len(vec1)/len(vec2)
}

################################################################################

orth <- function(A) { # Constructs orthogonal matrix from input matrix via 
                      # Gram-Schmidt
  A[, 1] <- A[, 1]/norm(A[, 1, drop = FALSE], type = 'F')
  
  for (k in 2:ncol(A)) {
    w <- t(A[, k, drop = FALSE]) %*% A[, 1:(k-1), drop = FALSE]
    for (j in 1:(k-1)) A[, k] <- A[, k] - w[j]*A[, j]
    A[, k] <- A[, k]/norm(A[, k, drop = FALSE], type = 'F')
  }
  return(A)
  
}

################################################################################

rowSubtract <- function(matrix, row) {

  t(t(matrix) - row)

}

################################################################################

indComb <- function(NN) {

  pt1 <- rep(1:NN, times = NN)
  pt2 <- rep(1:NN, each = NN)

  un <- pt1 > pt2

  pt1 <- pt1[un]
  pt2 <- pt2[un]

  data.frame(pt1 = pt2, pt2 = pt1)

}

################################################################################

indnComb <- function(NN, n) {

  if (n == 1) {
    return(matrix(1:NN, ncol = 1))
  }

  prev <- indnComb(NN, n-1)
  lastind <- prev[ , dim(prev)[2]]
  ind.cf1 <- rep(lastind, each = NN)
  ind.cf2 <- rep(1:NN, times = length(lastind))
  new.ind <- which(ind.cf1 < ind.cf2)
  new.ind1 <- ((new.ind - 1) %/% NN) + 1
  new.ind2 <- new.ind %% NN
  new.ind2[new.ind2 == 0] <- NN
  return(cbind(prev[new.ind1, , drop = FALSE], (1:NN)[new.ind2, drop = FALSE]))

}

################################################################################

addOpposite <- function(vecOneDir) { 

  vectors <- matrix(nrow = 2*nrow(vecOneDir), ncol = ncol(vecOneDir))
  odd <- rep(c(TRUE, FALSE), times = nrow(vecOneDir))
  vectors[odd, ] <- vecOneDir
  vectors[!odd, ] <- -vecOneDir
  return(vectors)

}

################################################################################

vecBw_onedir <- function(points) {

  NN <- dim(points)[1]
  i.c <- indComb(NN)
  return(points[i.c$pt1, , drop = FALSE] - points[i.c$pt2, , drop = FALSE])

}

################################################################################

vectorsBetween <- function(points) {

  vecOneDir <- vecBw_onedir(points)
  addOpposite(vecOneDir)

}

################################################################################

vecToC_onedir <- function(points, add.mids = FALSE, weight.mids = 1,
                          mids.max.dist = Inf) {

  # Mean center data
  center <- apply(points, 2, mean)
  vecOneDir <- rowSubtract(points, center)

  if (add.mids) { # Add midpoints

    i.c <- indComb(dim(vecOneDir)[1])
    mids <- (vecOneDir[i.c$pt1, ] + vecOneDir[i.c$pt2, ])/2
    dist <- lens(vecOneDir[i.c$pt1, ] - vecOneDir[i.c$pt2, ])
    mids <- mids[dist <= mids.max.dist, ] # Remove midpoints for very distant 
                                          # points
    vecOneDir <- rbind(vecOneDir, weight.mids * mids)

  }

  return(vecOneDir)

}
################################################################################

vectorsToCenter <- function(points, add.mids = FALSE, weight.mids = 1,
                            mids.max.dist = Inf) {

  vecOneDir <- vecToC_onedir(points, add.mids, weight.mids = 1, mids.max.dist)
  addOpposite(vecOneDir)

}

################################################################################

vectorsToCenter2 <- function(points, weight.mids = 1,
                             max.proj = 1) {

  # Mean center data
  center <- apply(points, 2, mean)
  vecOneDir <- rowSubtract(points, center)

  i.c <- indComb(dim(vecOneDir)[1])
  mids <- (vecOneDir[i.c$pt1, ] + vecOneDir[i.c$pt2, ])/2
  projections <- sapply(1:dim(vecOneDir)[1], 
                    function(ind, vecs1, vecs2) { 
                      normProj(vecOneDir[i.c$pt1, ], vecOneDir[i.c$pt2, ])
                    }, vecOneDir[i.c[[1]], ], vecOneDir[i.c[[2]], ])
  mids <- mids[abs(projections) <= max.proj, ] # Remove midpoints for vectors 
                                               # with very similar direction
  vecOneDir <- rbind(vecOneDir, weight.mids * mids)

  # Add opposite vectors

  addOpposite(vecOneDir)


}
################################################################################

vectorsToCenter3 <- function(points, weight.quarts = 1,
                             max.proj = 1) {

  # Mean center data
  center <- apply(points, 2, mean)
  vecOneDir <- rowSubtract(points, center)

  i.c <- indComb(dim(vecOneDir)[1])
  mids <- (vecOneDir[i.c$pt1, ] + vecOneDir[i.c$pt2, ])/2
  q1 <- 1/4*vecOneDir[i.c$pt1, ] + 3/4*vecOneDir[i.c$pt2, ]
  q2 <- 3/4*vecOneDir[i.c$pt1, ] + 1/4*vecOneDir[i.c$pt2, ]
  projections <- sapply(1:dim(vecOneDir)[1], 
                    function(ind, vecs1, vecs2) { 
                      normProj(vecOneDir[i.c$pt1, ], vecOneDir[i.c$pt2, ])
                    }, vecOneDir[i.c[[1]], ], vecOneDir[i.c[[2]], ])
  mids <- mids[abs(projections) <= max.proj, ] # Remove midpoints for vectors 
                                               # with very similar direction
  q1 <- q1[abs(projections) <= max.proj, ]
  q2 <- q2[abs(projections) <= max.proj, ]
  quarts <- rbind(mids, q1, q2)
  vecOneDir <- rbind(vecOneDir, weight.quarts * quarts)

  # Add opposite vectors

  addOpposite(vecOneDir)


}

################################################################################

vecToCs_onedir <- function(points, n.group) {

  if (n.group == 1) return(vecToC_onedir(points))

  NN <- dim(points)[1]
  ind.groups <- indnComb(NN, n.group)
  point.groups <- points[t(ind.groups), , drop = FALSE]
  group.centers <- t(apply(ind.groups, 1, 
                           function(ind.group) {
                             pts <- points[ind.group, ]
                             return(apply(pts, 2, mean))
                           } ))
  centers <- group.centers[rep(1:dim(group.centers)[1], each = n.group), ,
                           drop = FALSE]
  return(point.groups - centers)

}

################################################################################

vectorsToCenters <- function(points, n.group) {

  vec.one.dir <- vecToCs_onedir(points, n.group)
  addOpposite(vec.one.dir)

}

################################################################################

allsimpl <- function(points, sdim) {

  NN <- dim(points)[1]
  exdim <- dim(points)[2]
  ind.groups <- indnComb(NN, sdim + 1)
  print(dim(ind.groups))
  point.groups <- points[t(ind.groups), , drop = FALSE]
  basept <- seq(1, dim(point.groups)[1], by = sdim + 1)
  vecs <- point.groups[-basept, ] - point.groups[rep(basept, each = sdim), ]
  #point.groups.array <- aperm(array(point.groups, 
  #                            dim = c(sdim + 1, dim(ind.groups)[1], exdim)),
  #                            c(2, 1, 3))                                
  #print(point.groups)
  #print(point.groups.array)
  #vecs <- point.groups.array[, 1:sdim, ] - # fix this.... point.groups.array[, sdim + 1, ]
  #centers <- group.centers[rep(1:dim(group.centers)[1], each = n.group), ,
  #                         drop = FALSE]
  #vec.one.dir <- point.groups - centers
  #addOpposite(vec.one.dir)

}

################################################################################

condense <- function(ind, values, obsolete) {

  ind[obsolete] <- NA
  values[obsolete] <- NA

  for.ord <- ind

  if(is.logical(obsolete[1])) {
    for.ord[!obsolete] <- 0
  } else {
    for.ord[-obsolete] <- 0
  }

  ord <- t(apply(for.ord, 1, order))

  for (i in 1:dim(ind)[1]) {

    ind[i, ] <- ind[i, ord[i, ]]
    values[i, ] <- values[i, ord[i, ]]

  }
  list(ind, values)
}

################################################################################

depth <- function(index, data) {

  N <- nrow(data)
  d <- ncol(data)

  point <- matrix(data[index, ], nrow = N, ncol = d, byrow = TRUE)

  vectors <- data - point

  vec.len <- lens(vectors)

  zero.dist <- vec.len == 0

  vec.len[zero.dist] <- 1 # Avoid division by zero

  norm.vecs <- vectors/vec.len # Columnwise division

  1 - max(0, (len(apply(norm.vecs, 2, sum)) - sum(zero.dist))/N)

}


################################################################################

################################################################################

#points <- rbind(c(0, 1), c(0, 1), c(1, 1), c(0, 1))
#n.group <- 3
#plot(vectorsToCenters(points, n.group))
#points <- locator(4)
#points <- cbind(points$x, points$y)
#plot(points)
#plot(vectorsToCenters(points, 4))
#plot(vectorsToCenter(points))
#plot(vectorsToCenters(points, 2))
#plot(vectorsBetween(points))
#points <- matrix(1:20, nrow = 5, ncol = 4)
#allsimpl(points, 2)
