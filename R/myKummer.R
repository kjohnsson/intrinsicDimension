################################################################################

# NB! myKummer only works for positive values of x, which is not what we have
# when we want to compute the mean of the non-central Chi distribution. In that
# case use kummerMmod (se below)

# Defines functions: myKummerM, kummerapprox, kummerMmod, kummerM_sum,
# kummerM_sum_pos

################################################################################

myKummerM <- function(a, b, x) {

  if (a <= 0) stop("a must be positive")
  if (b <= a) stop("Must have b > a")
  gamma(b)/gamma(b-a)/gamma(a)*
    sapply(x, function(s) {
                integrate(function(t) {
                            res <- rep(NA, length(t))
                            boundary <- t == 0 | t == 1
                            res[boundary] <- 0
                            res[!boundary] <- 
                              exp(s*t + 
                                  log(t[!boundary])*(a-1) + 
                                  log(1-t[!boundary])*(b-a-1))
                            return(res)
                          }, 0, 1)$value
              })

}

#x <- seq(0, -200, length.out = 200)  # Faster that kummerMmod
#plot(x, myKummerM(1/2, 5/2, x), type = 'l')

################################################################################

#kummerapprox <- function(a, b, x) { # x large and positive
#
#  N <- 10
#  sapply(x, function(z) {
#    exp(z)*z^(a-b)*gamma(b)/gamma(a)*sum(c(1, 
#      sapply(1:N, function(n) {
#                    prod(b-a + 0:(n-1))*prod(1-a + 0:(n-1))/prod(1:n)/z^n
#                  })))})
#
#}

################################################################################

#library(fAsianOptions)
#
#kummerMmod <- function(a, b, x) {
#
#  if(any(!is.real(x))) stop('x must be real-valued')
#  res <- rep(NA, length(x))
#  neg <- x < 0
#  res[neg] <- exp(x[neg])*kummerM(-x[neg], b-a, b)
#  res[!neg] <- kummerM(x[!neg], a, b)
#  return(Re(res))
#
#}

#x <- seq(-200, 0, length.out = 200)
#plot(x, kummerMmod(1, 10, x), type = 'l')  # Compare to Maple
#lines(x, kummerapprox2(1, 10, x), col = 'red')

#system.time(replicate(100,kummerapprox2(1, 10, -x)))
#system.time(replicate(100,kummerM(1, 10, x)))

################################################################################

#kummerM_sum <- function(a, b, x) {
#
#  res <- rep(NA, length(x))
#  neg <- x < 0
#  if(any(neg)) res[neg] <- exp(x[neg])*kummerM_sum_pos(b-a, b, -x[neg])
#  if(any(!neg)) res[!neg] <- kummerM_sum_pos(a, b, x[!neg])
#  return(res)
#
#}
#
#kummerM_sum_pos <- function(a, b, x) {
#
#  if(any(x < 0)) stop('x must be non-negative')
#
#  tol <- 1e-10
#  nmin <- 10
#
#  term <- x*a/b
#  res <- 1 + term
#  n <- 1
#  an <- a
#  bn <- b
#  while(n < nmin | max(abs(term)) > tol) {
#    n <- n + 1
#    an <- an + 1
#    bn <- bn + 1
#    term <- x*term*an/bn/n
#    res <- res + term
#  }
#  return(res)
#
#}

################################################################################