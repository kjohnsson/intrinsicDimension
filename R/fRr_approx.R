# Functions to use for approximation of f(R_i|r) as a weighted Gaussian
# truncated at 0.
# Noise dimension: k; Noise sd: sigma

# Approximation (Gaussian)
#--------------------------
# w1 = inthr(R, k, sigma) (for scaling f(R|r) so that it integrates to 1)
# Mean: intrhr(R, k, sigma)/w1
# Sd: sqrt(intr2hr/w1 - mean^2)
# w2 = 1 - pnorm(0, mean = mu, sd = sqrt(tsigma)) Compensation for truncation
# Weight: w1/w2 (The trucated Gaussian should integrate to same as f(R_i|r))

################################################################################

#source('L:/Documents/Forskning/Code/TransPoiss/myKummer.R')
  # Load function: myKummerM

################################################################################

h <- function(R, k, r, sigma) 1/sigma*g(R/sigma, k, r/sigma) # f(R|r)
g <- function(x, k, lambda) 2*x*dchisq(x^2, k, lambda^2) # Density function
                                                # for non-central chi-distr.

################################################################################

inthr <- function(R, k, sigma) intglam(R/sigma, k)

intglam <- function(x, k) {
  y <- x^2/2
  res <- NA
  tryCatch({res <- sqrt(pi)*exp(-y)*y^((k-1)/2)/gamma(k/2)*
                     myKummerM(1/2, k/2, y)},
           error = function(ex) {
             cat('Did not manage to compute kummerM,', 
                    'numerical integration used instead.\n')
           }, finally = {if(is.na(res)) res <- numintglam(x, k)})
  return(res)         
}

intrhr <- function(R, k, sigma) {
  x <- R/sigma
  y <- x^2/2
  res <- NA
  tryCatch({res <- sigma*sqrt(2)*y^((k-1)/2)*exp(-y)/gamma(k/2)*
                     myKummerM(1, k/2, y)},
           error = function(ex) {
             cat('Did not manage to compute kummerM,', 
                    'numerical integration used instead.\n')
           }, finally = {if(is.na(res)) res <- numintrhr(R, k, sigma)})
  return(res)
}

intr2hr <- function(R, k, sigma) {
  x <- R/sigma
  y <- x^2/2
  res <- NA
  tryCatch({res <- sigma^2*sqrt(pi)*y^((k-1)/2)*exp(-y)/gamma(k/2)*
                     myKummerM(3/2, k/2, y)},
           error = function(ex) {
             cat('Did not manage to compute kummerM,', 
                    'numerical integration used instead.\n')
           }, finally = {if(is.na(res)) res <- numintr2hr(R, k, sigma)})
  return(res)

}


################################################################################

upintlim <- function(x, k) sqrt(qchisq(1 - 1e-5, df = k, ncp = x^2))

numintglam <- function(x, k) {
  integrate(function(lam) {
              g(x, k, lam)
            }, 0, upintlim(x, k))$value
}

numintrhr <- function(R, k, sigma) {
  integrate(function(r) {
              r*h(R, k, r, sigma)
            }, 0, upintlim(R/sigma, k))$value
}

numintr2hr <- function(R, k, sigma) {
  integrate(function(r) {
              r^2*h(R, k, r, sigma)
            }, 0, upintlim(R/sigma, k))$value
}

################################################################################


