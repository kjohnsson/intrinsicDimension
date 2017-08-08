\name{ide}
\alias{localIntrinsicDimension}
\alias{globalIntrinsicDimension}
\alias{pointwiseIntrinsicDimension}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Intrinsic Dimension Estimation
}
\description{
Intrinsic dimension estimation with method given as parameter.
}
\usage{
localIntrinsicDimension(.data, .method, ...)
globalIntrinsicDimension(.data, .method, ...)
pointwiseIntrinsicDimension(.data, .method, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.data}{ Data set for which dimension should be estimated.}
  \item{.method}{ For \code{local.dimension.estimation}, one of \code{'essDimEst'}, 
                 \code{'dancoDimEst'}, \code{'pcaLocalDimEst'},
                 \code{'maxLikLocalDimEst'}, 
                 \code{'knnDimEst'}. For \code{global.dimension.estimation}, one of
                 \code{'dancoDimEst'}, \code{'maxLikGlobalDimEst'}, \code{'knnDimEst'}.
                  For \code{pointwise.dimension.estimation},
                 \code{'pcaOtpmLocalDimEst'} or \code{'maxLikPointwiseDimEst'}.}
  \item{...}{ arguments passed to intrinsic dimension estimator. }
}
\value{
  For \code{localIntrinsicDimension} and \code{globalIntrinsicDimension}, a \code{DimEst} object with the slot \code{dim.est} containing the dimension estimate and possibly additional slots containing additional information about the estimation process.
%  
  For \code{pointwiseIntrinsicDimension}, a \code{DimEstPointwise} object, inhereting \code{data.frame}, with one slot \code{dim.est} containing the dimension estimates and possibly additional slots containing additional information about the estimation process.
}
\details{
  For the \code{localIntrinsicDimension} function, \code{.data} should be a
  local data set, i.e. a piece of a data set that is well approxiated by a
  hyperplane (meaning that the curvature should be low in the local data set).

  The function \code{pointwiseIntrinsicDimension} estimates local dimension
  around each data point in the data set.
}
\references{
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

  Levina, E. and Bickel., P. J. (2005) Maximum likelihood estimation of intrinsic dimension. \emph{Advances in Neural Information Processing Systems 17}, 777-784. MIT Press.

Carter, K.M., Raich, R. and Hero, A.O. (2010) On local intrinsic dimension 
  estimation and its applications. \emph{IEEE Trans. on Sig. Proc.}, 
  \bold{58}(2), 650-663.
}
\author{
Kerstin Johnsson, Lund University
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{essDimEst}}, \code{\link{dancoDimEst}}, \code{\link{pcaLocalDimEst}}, \code{\link{knnDimEst}}
\code{\link{pcaOtpmPointwiseDimEst}}, \code{\link{maxLikGlobalDimEst}}, \code{\link{maxLikLocalDimEst}},
\code{\link{maxLikPointwiseDimEst}}
}
\examples{
data <- hyperBall(100, 4, 15, .05)
localIntrinsicDimension(data, .method='essDimEst', ver = 's', d = 1)
globalIntrinsicDimension(data, 'dancoDimEst', k = 8, D = 20)
pointwiseIntrinsicDimension(data, .method='maxLikPointwiseDimEst', k = 8, dnoise = NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
