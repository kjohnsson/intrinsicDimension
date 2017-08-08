\name{swissRoll}
\alias{swissRoll}
\alias{swissRoll3Sph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Swiss roll with or without 3-sphere inside
}
\description{
  Generates a sample from a uniform distribution on a Swiss roll-surface, 
  possibly together with a sample from a uniform distribution on a 3-sphere
  inside the Swiss roll.
}
\usage{
swissRoll(Ns, a = 1, b = 2, nturn = 1.5, h = 4)
swissRoll3Sph(Ns, Nsph, a = 1, b = 2, nturn = 1.5, h = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ns}{ number of data points on the Swiss roll. }
  \item{Nsph}{ number of data points on the 3-sphere. }
  \item{a}{ minimal radius of Swiss roll and radius of 3-sphere. }
  \item{b}{ maximal radius of Swiss roll. }
  \item{nturn}{ number of turns of the surface. }
  \item{h}{ height of Swiss roll. }
}
%\details{
%
%}
\value{
  \code{swissRoll} returns three-dimensional data points.
  \code{swissRoll3Sph} returns four-dimensional data points with the Swiss roll
  in the three first dimensions (columns). The \code{Ns} first data points
  lie on the Swiss roll and the \code{Nsph} last data points lie on the 
  3-sphere.
}

\author{
  Kerstin Johnsson, Lund University.
}

\examples{
datap <- swissRoll3Sph(300, 100)
par(mfrow = c(1, 3))
plot(datap[,1], datap[,2])
plot(datap[,1], datap[,3])
plot(datap[,1], datap[,4])
}
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
