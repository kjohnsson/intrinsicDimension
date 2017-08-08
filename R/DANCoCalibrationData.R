DancoCalibrationData <- function(k, N) {
  me <- list(
    k = k,
    N = N,
    calibration.data = list(),
    maxdim = 0
  )
  class(me) <- append(class(me), 'DANCoCalibrationData')
  return(me)
}

print.DancoCalibrationData <- function(x, ...) {
  cat('DANCo Calibration Data with parameters:\n',
      'k =', x$k, ',', 'N =', x$N, ',',
      'maxdim =', x$maxdim, '\n')
}

increaseMaxDimByOne <- function(dancoCalDat) {
  newdim <- dancoCalDat$maxdim + 1
  MIND_MLx.maxdim <- newdim*2+5
  dancoCalDat$calibration.data[[newdim]] <-
    dancoDimEstNoCalibration(hyperBall(dancoCalDat$N, newdim), dancoCalDat$k, MIND_MLx.maxdim,
                ver = 'DANCo')
  dancoCalDat$maxdim <- newdim
  return(dancoCalDat)
}

