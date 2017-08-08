################################################################################

hyperBall <- function(Ns, d, n = d, sd = 0) {

  # Sample on unit sphere
  sph <- hyperSphere(Ns, d - 1)
  # Sample radius and multiply
  ball <- runif(Ns)^(1/d)*sph
  return(addNoise(ball, n, sd = sd))

}

################################################################################

hyperSphere <- function(Ns, d, n = d + 1, sd = 0) {

  sph <- matrix(rnorm(Ns * (d + 1)), nrow = Ns, ncol = d + 1)
  sph <- sph/lens(sph) # divides each column in sph with lens(sph)
  return(addNoise(sph, n, sd = sd))

}

################################################################################

isotropicNormal <- function(Ns, d, n = d, sd = 0) {
  nor <- matrix(rnorm(Ns * d), ncol = d)
  return(addNoise(nor, n, sd = sd))
}

################################################################################
