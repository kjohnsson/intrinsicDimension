################################################################################

swissRoll <- function(Ns, a = 1, b = 2, nturn = 1.5, h = 4)
  swissRoll3Sph(Ns, Nsph = 0, a, b, nturn, h)[, 1:3]

################################################################################

swissRoll3Sph <- function(Ns, Nsph, a = 1, b = 2, nturn = 1.5, h = 4) {

  if (Ns > 0) { # Sample from Swiss roll
  
    omega <- 2*pi*nturn

    # Swiss roll profile parametrization: 0 < r < 1
    #   x(r) = (a + b * r) * cos(omega * r)
    #   y(r) = (a + b * r) * sin(omega * r)
    # Length increment:
    #   dl(r) = sqrt( (dx/dr)^2 + (dy/dr)^2 )
    #         = sqrt( b^2 + omega^2 * (a + b * r)^2 )
    # To get uniform sampling on Swiss roll profile we need to sample r from
    # a pdf proprotional to dl(s), 0 < s < 1.

    # Rejection sampling from pdf proportional to dl(r)
    dl <- function(r) sqrt(b^2 + omega^2 * (a + b * r)^2)
    ok <- 0
    while(sum(ok) < Ns) {
      r.samp <- runif(3 * Ns)
      ok <-(runif(3 * Ns) < dl(r.samp)/dl(1))
    }
    r.samp <- r.samp[ok][1:Ns]

    # Map to Swiss roll profile
    x <- (a + b * r.samp) * cos(omega * r.samp)
    y <- (a + b * r.samp) * sin(omega * r.samp)

    # Construct coordinates for other dimensions
    z <- runif(Ns, -h, h)
    w <- rep(0, Ns)

  } else {

    x <- y <- z <- w <- NULL

  }

  if (Nsph > 0) { # Sample from 3-sphere

    sph <- hyperSphere(Nsph, 3) * a

    x <- c(x, sph[ , 1])
    y <- c(y, sph[ , 2])
    z <- c(z, sph[ , 3])
    w <- c(w, sph[ , 4])


  }

  cbind(x, y, z, w)

}

################################################################################

cornerPlane = function(Ns, theta = pi/4) {

  x <- runif(Ns, -1, 1)
  y <- rep(0, Ns)
  z <- runif(Ns, -1, 1)

  rigid <- (x <= 0)
  bend <- (x > 0)

  x.ri <- x[rigid]
  y.ri <- y[rigid]

  x.be.bef <- x[bend]
  y.be.bef <- y[bend]

  x.be <- x.be.bef * cos(pi - theta)
  y.be <- x.be.bef * sin(pi - theta)

  x <- c(x.ri,x.be)
  y <- c(y.ri,y.be)

  cbind(x, y, z)

}

################################################################################
