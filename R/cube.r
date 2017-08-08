################################################################################

hyperCube <- function(Ns, n, side = 1) {

	matrix(side*runif(Ns * n), ncol = n)

}

################################################################################

hyperCubeFaces <- function(Ns, n) hyperCubeEdges(Ns, n-1, n)

################################################################################

hyperCubeEdges <- function(Ns, d, n) {

  coord <- matrix(runif(Ns*n), nrow = Ns, ncol = n)

  for (k in 1:Ns) {

    which.singular <- sample(1:n, n - d)
    sing.vals <- sample(c(0, 1), n - d, replace = TRUE)
    coord[k, which.singular] <- sing.vals

  }

  return(coord)

}

################################################################################
