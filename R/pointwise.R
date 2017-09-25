asPointwiseEstimator <- function(estimator, neighborhood.size, indices=NULL, eps=0.0) {
  
  pointwiseEstimator <- function(data, ...) {
    if (is.null(indices)) indices <- 1:dim(data)[1]
    local.data.sets <- neighborhoods(data, neighborhood.size, indices, eps)
    dim.est.res <- lapply(local.data.sets, FUN=function(loc.data) estimator(loc.data, ...))
    dim.ests <- sapply(dim.est.res, function(res) res$dim.est)
    additional.data.names <- names(dim.est.res[[1]])[-1]
    pointwise.additional.data <- list()
    for (name in additional.data.names) {
      pointwise.additional.data[[name]] <- sapply(dim.est.res, function(res) res[[name]])
    }
    return(DimEstPointwise(dim.ests, pointwise.additional.data))
  }
  
  return(pointwiseEstimator)
}

################################################################################

yaIkNN <- function(data, indexes, k, eps=0.0) {
  # ann returns indexes and SQUARED Euclidean distances
  # (not in accordance with documentation)
  # This function corrects for this and returns Euclidean distances
  # The function returns indices and distances to a neighborhood around each
  # data point in indexes, 
  kNN.index.dist <- ann(data, data[indexes, , drop = FALSE], k + 1, verbose = FALSE, 
                        eps=eps)$knnIndexDist[ , -c(1, k + 2), drop = FALSE]
  kNN.dist <- sqrt(kNN.index.dist[ , (k + 1):(2*k)])
  kNN.index.dist[ , (k + 1):(2*k)] <- kNN.dist
  return(kNN.index.dist)
}

################################################################################

neighborhoods <- function(data, k, indices, eps=0.0) {
  nbh.indices.and.dist <- yaIkNN(data, indices, k, eps)
  nbh.indices <- nbh.indices.and.dist[ , 1:k]
  nbh <- lapply(split(nbh.indices, row(nbh.indices)), FUN=function(nb.ind) data[nb.ind, ])
  return(nbh)
}

