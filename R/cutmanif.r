################################################################################

cutHyperPlane <- function(Ns, d, n, sd) {
  # Cut-off radius is 1
  if (n < d) stop("n has to be at least d")
  N <- round(Ns*20*n/4 * max(sd, 0.1))
  k <- 0
  data <- NULL
  if (n*sd < 1e-4) return(hyperBall(Ns, d, n))
  repeat {
    k <- k + 1
    if (k > 50) {
      warning('Not enough points after max nbr of iterations')
      return(NULL)
    }
    tryCatch({
      new.data <- hyperBall(N, d)
      #new.data <- 2*(hyperCube(N, d) - 0.5)
      new.data <- cbind(new.data, 
                        addNoise(matrix(0, nrow = N, ncol = n - d), sd = sd))
    }, error = function(ex) {
      N <- N/2
      new.data <- hyperBall(N, d)
      #new.data <- 2*(hyperCube(N, d) - 0.5)
      new.data <- cbind(new.data, 
                        addNoise(matrix(0, nrow = N, ncol = n - d), sd = sd))
    })
    data <- rbind(data, new.data[lens(new.data) < 1, , drop = FALSE])
    if (dim(data)[1] >= Ns) break
    if (k <= 5) N <- N*2
  }
  return(data[1:Ns, ])
  
}

################################################################################

cutHyperSphere <- function(Ns, rat, d, n, sd) {

  if (n < d+1) stop("n has to be at least d+1")
  N <- Ns*40*n * max(sd, 0.1)
  k <- 0
  data <- NULL
  repeat {
    k <- k + 1
    if (k > 10) {
      warning('Not enough points after max nbr of iterations, no data returned')
      return(NULL)
    }
    new.data <- addNoise(hyperSphere(N, d), n, sd)
    
    if (rat <= sqrt(2-sqrt(2)))  {
      # Use symmetry to increase data in 1st quadrant
      pos.data <- abs(new.data[, 1:(d+1)])
      if (n > d + 1) new.data <- cbind(pos.data, new.data[,(d+2):n])
      else new.data <- pos.data
    }
    
    mid <- c(rep(1/sqrt(d+1), d+1), rep(0, n - d - 1))
    new.data <- t(t(new.data) - mid)
    new.data <- new.data[lens(new.data) < rat, ]
    data <- rbind(data, new.data)
    if (dim(data)[1] >= Ns) break
  }
  return(data[1:Ns, ])

}

################################################################################