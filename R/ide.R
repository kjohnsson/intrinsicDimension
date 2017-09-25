localIntrinsicDimension <- function(.data, .method, ...) {
  IDEs <- c('essLocalDimEst', 'dancoDimEst', 'pcaLocalDimEst', 'maxLikLocalDimEst', 'kNN')
  if (!(.method %in% IDEs)) {
    stop(paste(.method, 'not local intrinsic dimension estimator. Should be one of',
               paste(IDEs, collapse=',')))
  }
  return(do.call(.method, c(list(.data), list(...))))  
}

globalIntrinsicDimension <- function(.data, .method, ...) {
  IDEs <- c('dancoDimEst', 'maxLikGlobalDimEst', 'kNN')
  if (!(.method %in% IDEs)) {
    stop(paste(.method, 'not global intrinsic dimension estimator. Should be one of',
               paste(IDEs, collapse=',')))
  }
  return(do.call(.method, c(list(.data), list(...))))  
}

pointwiseIntrinsicDimension <- function(.data, .method, ...) {
  IDEs <- c('pcaOtpmPointwiseDimEst', 'maxLikPointwiseDimEst')
  if (!(.method %in% IDEs)) {
    stop(paste(.method, 'not pointwise intrinsic dimension estimator. Should be one of',
               paste(IDEs, collapse=',')))
  }
  return(do.call(.method, c(list(.data), list(...))))    
}

