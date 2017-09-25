\name{essLocalDimEst}
\alias{essLocalDimEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Expected Simplex Skewness Local Dimension Estimation
}
\description{
Local intrinsic dimension estimation with the ESS method
}
\usage{
essLocalDimEst(data, ver, d = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Local data set for which dimension should be estimated.}
  \item{ver}{ Possible values: \code{'a'} and \code{'b'}. See Johnsson et al. (2015).}
  \item{d}{ For \code{ver = 'a'}, any value of \code{d} is possible, for \code{ver = 'b'}, only \code{d = 1} is supported. }
}
\details{
The ESS method assumes that the data is local, i.e. that it is a neighborhood taken from a larger data set,
such that the curvature and the noise within the neighborhood is relatively small. In the ideal case
(no noise, no curvature) this is equivalent to the data being uniformly distributed over a hyper ball.
}
\value{
A \code{DimEst} object with two slots:
\item{dim.est}{The interpolated dimension estimate.}
\item{ess}{The ESS value produced by the algorithm.}
}
\references{
Johnsson, K., Soneson, C., & Fontes, M. (2015). Low Bias Local Intrinsic Dimension Estimation from Expected Simplex Skewness. IEEE Trans. Pattern Anal. Mach. Intell., 37(1), 196-202.
}
\author{
Kerstin Johnsson, Lund University
}
\examples{
data <- hyperBall(100, 4, 15, .05)
essLocalDimEst(data, ver = 'a', d = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
