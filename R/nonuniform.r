################################################################################

twinPeaks <- function(Ns, h = 1) hyperTwinPeaks(Ns, 2, h)

################################################################################

hyperTwinPeaks <- function(Ns, n, h = 1) {

  base.coord <- matrix(runif(Ns*n), nrow = Ns, ncol = n)
  height <- h*apply(sin(2*pi*base.coord), 1, prod)
  return(cbind(base.coord, height))

}

################################################################################

oblongNormal <- function(Ns, n) {

  cbind(matrix(rnorm(Ns*ceiling(n/2)), nrow = Ns, ncol = ceiling(n/2)),
        matrix(rnorm(Ns*floor(n/2), sd = sqrt(0.25)), nrow = Ns,
               ncol = floor(n/2)))

}

################################################################################

mHeinManifold <- function(Ns) {

  alpha <- matrix(runif(Ns*12), ncol = 12)
  x <- matrix(nrow = Ns, ncol = 72)
  x[, 2*(1:11)-1] <- alpha[, 2:12]*cos(2*pi*alpha[, 1:11])
  x[, 2*(1:11)] <- alpha[, 2:12]*sin(2*pi*alpha[, 1:11])
  x[, 23] <- alpha[, 1]*cos(2*pi*alpha[, 12])
  x[, 24] <- alpha[, 1]*sin(2*pi*alpha[, 12])
  x[, 1:24 + 24] <- x[, 1:24]
  x[, 1:24 + 48] <- x[, 1:24]
  return(x) 

}

################################################################################

m14Manifold <- function(Ns) M_rozza(Ns, 18)
m15Manifold <- function(Ns) M_rozza(Ns, 24)
  
################################################################################

M_rozza <- function(Ns, d) {

  X <- matrix(runif(Ns*d), ncol = d)
  D <- X*sin(cos(2*pi*X))
  DD <- X*cos(sin(2*pi*X))
  DDD <- cbind(D, DD)
  return(cbind(DDD, DDD))
  
}


################################################################################