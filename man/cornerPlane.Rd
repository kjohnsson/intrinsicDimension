\name{cornerPlane}
\alias{cornerPlane}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Corner Plane }
\description{
  Generates a sample from a uniform distribution on a bent plane. Half of the
  plane is in the xz-plane and half of the plane is bent over the x-axis, so
  that the resulting surface has an edge along the x-axis.
}
\usage{
cornerPlane(Ns, theta = pi/4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ns}{ number of data points. }
  \item{theta}{ angle at the x-axis. }
}
\value{
A \code{Ns} x 3 matrix with columns x, y and z.
}
\author{
  Kerstin Johnsson, Lund University
}
\examples{
datap <- cornerPlane(400)
par(mfrow = c(1, 2))
plot(datap[,1], datap[,2])
plot(datap[,1], datap[,3])
}
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
