################################################################################

len <- function(vector) sqrt(sum(vector^2))

################################################################################

knnDimEst <- function(data, k, ps, M, gamma = 2) {

  n <- nrow(data)
  Q <- length(ps)

  # Compute the distance between any two points in the data set

  dist <- matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:(i - 1)) {

      dist[i, j] = len(data[i, ] - data[j, ])

    }
  }

  dist <- dist + t(dist)

  # Compute weighted graph length for each sample

  L <- matrix(0, nrow = Q, ncol = M)

  for (i in 1:Q) {
    for (j in 1:M) {

      samp.ind <- sample(1:n, ps[i])
      
      for (l in samp.ind) {

        L[i, j] <- L[i, j] + sum(sort(dist[l, samp.ind])[2:(k + 1)]^gamma)
          # Add the weighted sum of the distances to the k nearest neighbors.
          # We should not include the sample itself, to which the distance is 
          # zero. 

      }
    }
  }

  # Least squares solution for m

  d <- ncol(data)
  epsilon <- rep(NA, d)

  for (m in 1:d) {

    alpha = (m - gamma)/m
    ps.alpha = ps^alpha
    hat.c = sum(ps.alpha * apply(L, 1, sum)) / (sum(ps.alpha^2)*M)
    epsilon[m] = sum((L - matrix(hat.c*ps.alpha, nrow = Q, ncol = M))^2)
      # matrix(vec, nrow = length(vec), ncol = b) is a matrix with b
      # identical columns equal to vec
      # sum(matr) is the sum of all elements in the matrix matr

  }

  dest <- which.min(epsilon) # Missing values are discarded
  return(DimEst(dest, residual = epsilon[dest]))

}

################################################################################
