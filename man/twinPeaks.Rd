\name{twinPeaks}
\alias{twinPeaks}
\alias{hyperTwinPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Twin Peaks }
\description{
Generates data points from a two- or higher-dimensional Twin Peaks manifold.
}
\usage{
twinPeaks(Ns, h = 1)
hyperTwinPeaks(Ns, n, h = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ns}{ number of data points. }
  \item{n}{ dimension of the (hyper) plane from which the peaks stand out. For
            \code{twinPeaks} \code{n} is 2. }
  \item{h}{ height of the peaks. }
}
\details{
  The height of the points is computed as \eqn{\prod_1^n \sin(x_i)}{
  prod(sin(x_i))}, where \eqn{x_1,...,x_n} are the coordinates of the point in
  the (hyper) plane.
}
\value{
  A \code{n+1}-dimensional data set, where the last dimension represents the
  height of the points.
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
  Kerstin Johnsson, Lund University.
}

\examples{
datap <- twinPeaks(400)
par(mfrow = c(1, 3))
plot(datap[,1], datap[,3])
plot(datap[,2], datap[,3])
plot(datap[,1], datap[,2])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
