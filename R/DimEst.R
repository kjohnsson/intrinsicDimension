DimEst <- function(dim.est, ...) {
  if (!is.numeric(dim.est) || length(dim.est) > 1) {
    stop("Non-numeric or too long input")
  }
  self <- c(list("dim.est" = dim.est), list(...))
  class(self) <- append(class(self), "DimEst")
  return(self)
}

print.DimEst <- function(x, ...) {
  cat("Dimension estimate:", x$dim.est, "\n")
  if (length(x) > 1) {
    cat("Additional data:", paste(names(x)[-1], collapse=', '), "\n")
  }
}

################

DimEstPointwise <- function(dim.ests, pointwise.data = NULL) {
  
  if (!is.numeric(dim.ests)) {
    stop("Non-numeric input")
  }
  self <- data.frame("dim.est" = dim.ests)
  class(self) <- c("DimEstPointwise", class(self))
  
  if (is.null(pointwise.data)) return(self)
  
  for (data in names(pointwise.data)) {
    self[data] <- pointwise.data[data]
  } 
  return(self)
  
}

print.DimEstPointwise <- function(x, ...) {
  cat("Dimension estimates at", dim(x)[1], "data points.\n",
      "min:", min(x$dim.est), "; max:", max(x), "\n")
  if (dim(x)[2] > 1)
    cat("Additional data:", paste(names(x)[-1], collapse=', '), '\n')
}

