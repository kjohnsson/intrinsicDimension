\name{intrinsicDimension-package}
\alias{intrinsicDimension}
\docType{package}
\title{
  Intrinsic Dimension Estimation and Data on Manifolds
}
\description{
The intrinsic dimension of a data set typically means the manifold or Hausdorff dimension of the support of the probability distribution generating the data. It can also mean the fractal dimension of the support. See references for details.
}
\details{
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
  Functions for estimating intrinsic dimension: \code{\link{localIntrinsicDimension}}, \code{\link{globalIntrinsicDimension}}, \code{\link{pointwiseIntrinsicDimension}}, \code{\link{essDimEst}}, \code{\link{dancoDimEst}}, \code{\link{pcaLocalDimEst}}, \code{\link{pcaOtpmPointwiseDimEst}},
  \code{\link{maxLikGlobalDimEst}}, \code{\link{maxLikLocalDimEst}}, \code{\link{maxLikPointwiseDimEst}}, \code{\link{knnDimEst}}.

  Functions for generating data points from (usually uniform) distributions on 
  manifolds (possibly with noise): \code{\link{hyperBall}}, \code{\link{hyperSphere}}, \code{\link{hyperCube}},
  \code{\link{isotropicNormal}}, \code{\link{hyperCubeFaces}}, \code{\link{hyperCubeEdges}}, \code{\link{cutHyperPlane}}, \code{\link{cutHyperSphere}}, \code{\link{oblongNormal}}, \code{\link{swissRoll}}, \code{\link{swissRoll3Sph}}, \code{\link{twinPeaks}}, \code{\link{hyperTwinPeaks}},
  \code{\link{cornerPlane}}, \code{\link{mHeinManifold}}, \code{\link{m14Manifold}}, \code{\link{m15Manifold}}.
}
\author{
Kerstin Johnsson, Lund University

Maintainer: Kerstin Johnsson <kerstin.johnsson@hotmail.com>
}
\references{
%~~ Literature or other references for background information ~~
Johnsson, K (2016). Structures in high-dimensional data: Intrinsic dimension and cluster analysis. PhD thesis. Lund University.

Johnsson, K., Soneson, C. and Fontes, M. (2015). Low Bias Local Intrinsic 
  Dimension Estimation from Expected Simplex Skewness. \emph{IEEE Trans. Pattern Anal. 
  Mach. Intell.}, \bold{37}(1), 196-202.

Ceruti, C. et al. (2012). DANCo: Dimensionality from Angle and Norm Concentration.
  \emph{arXiv preprint} 1206.3881.
  
Rozza, A et al. (2012). Novel high intrinsic dimensionality estimators. \emph{Machine learning}
  \bold{89}, 37-65.

Fukunaga, K. and Olsen, D. R. (1971). An algorithm for finding intrinsic dimensionality
  of data. \emph{IEEE Trans. Comput.}, \bold{c-20}(2):176-183.

Fan, M. et al. (2010). Intrinsic dimension estimation of data by principal component 
  analysis. \emph{arXiv preprint} 1002.2050.

Bruske, J. and Sommer, G. (1998) Intrinsic dimensionality estimation with
  optimally topology perserving maps. \emph{IEEE Trans. on Pattern Anal. and Mach.
  Intell.}, \bold{20}(5), 572-575.

Haro, G., Randall, G. and Sapiro, G. (2008) Translated Poisson Mixture Model
  for Stratification Learning. \emph{Int. J. Comput. Vis.}, \bold{80}, 358-374.

Hill, B. M. (1975) A simple general approach to inference about the tail of a distribution. \emph{Ann. Stat.}, \bold{3}(5) 1163-1174.

Carter, K.M., Raich, R. and Hero, A.O. (2010) On local intrinsic dimension 
  estimation and its applications. \emph{IEEE Trans. on Sig. Proc.}, 
  \bold{58}(2), 650-663.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\keyword{ datagen }
\examples{
%~~ simple examples of the most important functions ~~
data <- hyperBall(100, 4, 15, .05)
localIntrinsicDimension(data, .method='essDimEst', ver = 's', d = 1)
}
