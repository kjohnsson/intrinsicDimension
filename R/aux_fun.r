################################################################################

lens <- function(vectors) {
  sqrt(apply(vectors^2, 1, sum))
}

################################################################################

embed <- function(data, n) {

  add.dim = n - ncol(data)
  if (add.dim == 0) return(data)
  if (add.dim < 0)
    stop('Manifold dimension greater than embedding dim')
  cbind(data, matrix(0, ncol = add.dim, nrow = nrow(data)))

}

################################################################################

addNoise <- function(data, n = ncol(data), sd) {
  if (n < ncol(data))
    stop("Dimension of noise must be larger than or equal to dimension of data")
  data <- embed(data, n)
  return(data + matrix(rnorm(prod(dim(data)), sd = sd),
                        nrow = nrow(data), ncol = ncol(data)))
}

################################################################################
