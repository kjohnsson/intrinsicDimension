\name{asPointwiseEstimator}
\alias{asPointwiseEstimator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Turn a local estimator into a pointwise estimator
}
\description{
Returns a function that can be used as a pointwise estimator of intrinsic dimension that uses local data sets with a fixed number of data points.
}
\usage{
asPointwiseEstimator(estimator, neighborhood.size, indices=NULL, eps=0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimator}{ A local intrinsic dimension estimator.}
  \item{neighborhood.size}{ The number of neighbors used for each dimension estimate.}
  \item{indices}{ A vector with indices of the points in \code{data} (as sent to the estimator function) that should be used as center for neighborhoods. }
  \item{eps}{ If non-zero, the relative error in distance allowed when finding nearest neighbors. See Details.}
}
\details{
The \code{ann} function of the package \code{yaImpute} is used for finding the \code{k} nearest neighbors. The \code{eps} parameter to \code{neighborhoods} is used in the \code{ann} function.
}
\value{
A function that can be used as a pointwise dimension estimator.
}
\author{
Kerstin Johnsson, Lund University
}
\examples{
data <- swissRoll3Sph(300, 300) 
 # the first 300 data points are on the swiss roll (ID=2) , the last 300 on the 3-sphere (ID=3)
essPointwiseDimEst <- asPointwiseEstimator(essLocalDimEst, neighborhood.size=10,
                                           indices = c(1:10, 301:310))
ess.pw.res <- essPointwiseDimEst(data)
ess.pw.res$dim.est
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
