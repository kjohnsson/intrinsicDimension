\name{essReference}
\alias{essReference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ESS Reference Values
}
\description{
Reference values for the ESS dimension estimation method
}
\usage{
essReference(ver, d, maxdim, mindim=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ver}{ Possible values: \code{'a'} and \code{'b'}. See Johnsson et al. (2015).}
  \item{d}{ For \code{ver = 'a'}, any value of \code{d} is possible, for \code{ver = 'b'}, only \code{d = 1} is supported. }
  \item{maxdim}{ The largest dimension for which reference values should be computed. }
  \item{mindim}{ The smallest dimension for which reference values should be computed. }
}
\details{
The ESS reference values are used by the ESS algorithm (\link{essLocalDimEst}) to compute the final dimension estimate.
}
\value{
A vector of length maxdim-(mindim-1), where each slot represents the reference value.
}
\references{
Johnsson, K., Soneson, C., & Fontes, M. (2015). Low Bias Local Intrinsic Dimension Estimation from Expected Simplex Skewness. IEEE Trans. Pattern Anal. Mach. Intell., 37(1), 196-202.
}
\author{
Kerstin Johnsson, Lund University
}
\examples{
essReference('a', 3, maxdim=500)
essReference('b', 1, maxdim=30, mindim=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
