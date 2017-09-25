\name{intrinsicDimension-package}
\alias{intrinsicDimension}
\docType{package}
\title{
  Intrinsic Dimension Estimation and Data on Manifolds
}
\description{
The intrinsic dimension of a data set is a measure of its complexity. In technical terms it typically means the manifold or Hausdorff (fractal) dimension of the support of the probability distribution generating the data. This package contains functions for estimating intrinsic dimension and generating ground truth data sets with known intrinsic dimension.
}
\details{
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
Data sets that can be accurately described with a few parameters have low intrinsic dimension. It is expected that the performance of many machine learning algorithms is dependent on the intrinsic dimension of the data. Is has also been proposed to use estimates of intrinsic dimension for applications such as network anomaly detection and image analysis.

This package contains implementations of a variety of approaches for intrinsic dimension estimation: modeling distances by for example Maximum Likelihood, approximating hyperplanes using Principal Component Analysis (PCA) and modeling angular information and concentration of measure (ESS and DANCo methods). Ground truth data, i.e. data with known intrinsic dimension, can be generated with a number of functions modeling manifolds. The manifold dimension is the intrinsic dimension.

The package distinguishes between local, global and pointwise estimators of intrinsic dimension. Local estimators estimate dimension of a _local data set_, for example a neighborhood from a larger data set. For this estimate to be accurate the noise and the curvature of the data has to be small relative to the neighborhood diameter. A global estimator takes the entire data set and returns one estimate of intrinsic dimension. Global estimators has the potential to handle higher noise and curvature levels than local estimators, but require that the entire data set has the same intrinsic dimension. Pointwise estimators are essentially local estimators applied neighborhoods around each point in the data set, but sometimes information beyond the neighborhood is used, as in PCA with Optimally Topolgy Perserving Maps. Any local estimator can be converted into a pointwise estimator.

  Functions for estimating intrinsic dimension: \code{\link{localIntrinsicDimension}}, \code{\link{globalIntrinsicDimension}}, \code{\link{pointwiseIntrinsicDimension}}, \code{\link{essLocalDimEst}}, \code{\link{dancoDimEst}}, \code{\link{pcaLocalDimEst}}, \code{\link{pcaOtpmPointwiseDimEst}},
  \code{\link{maxLikGlobalDimEst}}, \code{\link{maxLikLocalDimEst}}, \code{\link{maxLikPointwiseDimEst}}, \code{\link{knnDimEst}}.

  Functions for generating data points from (usually uniform) distributions on 
  manifolds (possibly with noise): \code{\link{hyperBall}}, \code{\link{hyperSphere}}, \code{\link{hyperCube}},
  \code{\link{isotropicNormal}}, \code{\link{hyperCubeFaces}}, \code{\link{hyperCubeEdges}}, \code{\link{cutHyperPlane}}, \code{\link{cutHyperSphere}}, \code{\link{oblongNormal}}, \code{\link{swissRoll}}, \code{\link{swissRoll3Sph}}, \code{\link{twinPeaks}}, \code{\link{hyperTwinPeaks}},
  \code{\link{cornerPlane}}, \code{\link{mHeinManifold}}, \code{\link{m14Manifold}}, \code{\link{m15Manifold}}.

  Functions for applying local estimators to non-local data: \code{\link{asPointwiseEstimator}}, \code{\link{neighborhoods}}
}
\author{
Kerstin Johnsson, Lund University

Maintainer: Kerstin Johnsson <kerstin.johnsson@hotmail.com>
}
\references{
%~~ Literature or other references for background information ~~
Johnsson, K (2016). Structures in high-dimensional data: Intrinsic dimension and cluster analysis. PhD thesis. Lund University. \url{http://portal.research.lu.se/portal/sv/publications/structures-in-highdimensional-data-intrinsic-dimension-and-cluster-analysis(8404f72e-e760-436d-ad7f-1be15af4b3d1).html}
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\keyword{ datagen }
%~~ simple examples of the most important functions ~~
