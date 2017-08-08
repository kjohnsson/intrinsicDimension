\name{mHeinManifold}
\alias{mHeinManifold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 12-dimensional manifold from Hein and Audibert (2005) }
\description{
  Generates a 12-dimensional manifold with extrinsic dimension 72
  (not uniformly sampled).
}
\usage{
mHeinManifold(Ns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ns}{ number of data points. }
}
%\details{
%  The height of the points is computed as \eqn{\prod_1^n \sin(x_i)}{
%  prod(sin(x_i))}, where \eqn{x_1,...,x_n} are the coordinates of the point in
%  the (hyper) plane.
%}
\value{
  A \code{72}-dimensional data set.
}
\references{
  Hein, M. and Audibert, J-Y. (2005) Intrinsic Dimensionality Estimation of
  Submanifolds in R^d. \emph{Proceedings of ICML}, 289-296.
}
\author{
  Kerstin Johnsson, Lund University.
}

\examples{
datap <- mHeinManifold(800)
par(mfrow = c(1, 3))
plot(datap[,1], datap[,3])
plot(datap[,2], datap[,3])
plot(datap[,1], datap[,2])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
