\name{hyperCube}
\alias{hyperCube}
\alias{hyperCubeFaces}
\alias{hyperCubeEdges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hypercube }
\description{
  Generates a sample from a uniform distribution on a hypercube, the faces
  of a hypercube or the ``edges'' of a hyper cube.
}
\usage{
hyperCube(Ns, n, side = 1)
hyperCubeFaces(Ns, n)
hyperCubeEdges(Ns, d, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ns}{ number of data points. }
  \item{d}{ dimension of edges. }
  \item{n}{ dimension of the hypercube. }
  \item{side}{ the length of the side of the hyper cube. }
}
\details{
  The hypercube is \eqn{[0,1]^n}. The edges of dimension \code{d} of the 
  hypercube are the \code{d}-dimensional boundaries of the hypercube. The 
  hypercube faces are the hyper cube edges of dimension \code{n-1}. 
}
\value{
A \code{Ns} by \code{n} matrix.
}
\author{
  Kerstin Johnsson, Lund University.
}

\examples{
datap <- hyperCubeEdges(200, 1, 3)
par(mfrow = c(1, 3))
plot(datap[, 1], datap[, 2])
plot(datap[, 1], datap[, 3])
plot(datap[, 2], datap[, 3])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
