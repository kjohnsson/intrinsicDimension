library(testthat)
library(intrinsicDimension)

#test_check("intrinsicDimension")

################################################################################

vnball <- function(n) return(pi^(n/2)/gamma(n/2+1))
vnsphere <- function(n) return((n+1)*vnball(n+1))

essReferenceNonStable <- function(ver, d, maxdim, mindim=1) {
  # this reference was used in version 1.0 of the package
  if (ver == 'a') {
    if (d == 1) {
      n <- max(mindim, 2):maxdim
      dim.val <- gamma(n/2)^2/(gamma((n-1)/2)*gamma((n+1)/2))
    } else {
      n <- max(mindim, (d+1)):maxdim
      dim.val <- (gamma(n/2)/gamma((n+1)/2))^d*gamma(n/2)/gamma((n-d)/2)
    }
    if (mindim <= d) {
      dim.val <- c(rep(0, d-mindim+1), dim.val)
    }
    return(dim.val)
  } 
  if (ver == 'b') {
    if (d == 1) {
      n <- mindim:maxdim
      return(2*vnball(n-1)/vnsphere(n-1))
    }
    stop('For ver == "b", d > 1 is not supported.')
  }
  stop('Not a valid version')
  
}

################################################################################

# M <- 210
# m <- 190
# plot(essReference('a', 1, M)[m:M])
# points(essReferenceNonStable('a', 1, M)[m:M], col='blue', pch='+')

test_that("ESSa values with stable algorithm equals values with original 
          algorithm", {
  original.stable.to <- 190
  for (d in 1:(original.stable.to-3)) {
    expect_equal(essReferenceNonStable('a', d, original.stable.to),
                 essReference('a', d, original.stable.to), info = sprintf("d = %d", d))
  }
})

# M <- 350
# m <- 340
# plot(essReference('b', 1, M)[m:M])
# points(essReferenceNonStable('b', 1, M)[m:M], col='blue', pch='+')

test_that("ESSb values with stable algorithm equals values with original 
          algorithm", {
            original.stable.to <- 340
            expect_equal(essReferenceNonStable('b', 1, original.stable.to),
                         essReference('b', 1, original.stable.to))
          })