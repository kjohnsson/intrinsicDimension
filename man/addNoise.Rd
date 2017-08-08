\name{addNoise}
\alias{addNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Add Noise to Data Set
}
\description{
  Embeds the data in \code{n} dimensions and adds normal isotropic noise to 
  the data set. Hence \code{n} has to be at least equal to the dimension (the
  number of columns) of the data set, otherwise the function terminates with 
  an error.
}
\usage{
addNoise(data, n = ncol(data), sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ data set. Each row corresponds to a data point. }
  \item{n}{ dimension of noise. }
  \item{sd}{ standard deviation of noise. The covariance matrix of the noise
             is \eqn{sd^2 \cdot I}{sd^2*I}.}
}
\value{
Matrix of same size as \code{data}.
}
\author{
  Kerstin Johnsson, Lund University
}

\examples{
datap <- hyperCubeEdges(100, 1, 2)
datap <- addNoise(datap, 3, .05)
par(mfrow = c(1, 2))
plot(datap[, 1], datap[, 2])
plot(datap[, 1], datap[, 3])
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
