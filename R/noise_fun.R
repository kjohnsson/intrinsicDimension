################################################################################

dnoiseGaussH <- function(r, s, sigma, k = NULL) {

  if (length(r) > 1 && length(s) > 1) 
    stop('r and s cannot both be vectors')
  dnorm(s, r, sigma)              # f(s|r) in Haro et al. (2008) w/ Gaussian
                                  # transition density
                                  # 'k' is not used, but is input
                                  # for compatibility

}

################################################################################

dnoiseGaussB <- function(r, s, sigma, k) {

  if (length(r) > 1 && length(s) > 1) 
    stop('r and s cannot both be vectors')
  w1 <- inthr(s, k, sigma)
  mu <- intrhr(s, k, sigma)/w1
  tsigma2 <- intr2hr(s, k, sigma)/w1 - mu^2
  w2 <- 1 - pnorm(0, mean = mu, sd = sqrt(tsigma2))
  w <- w1/w2
  gau <- w*dnorm(r, mu, sqrt(tsigma2))
  return(gau*(r > 0))   # Best approximation of f(s|r) by truncated Gaussian
                        # when Gaussian k-dim noise.
                        # NB! The noncentral chi distribution is
                        # f(r|s), not f(s|r).

}

################################################################################

dnoiseNcChi <- function(r, s, sigma, k) {

  if (length(r) > 1 && length(s) > 1) 
    stop('r and s cannot both be vectors')
  lambda <- r/sigma
  2*s/sigma^2*dchisq((s/sigma)^2, k, lambda^2) # 'k' is the number of
                                               # degrees of freedom

}

################################################################################

