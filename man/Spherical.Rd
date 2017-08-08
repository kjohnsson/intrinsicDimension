\name{Spherical}
\alias{hyperBall}
\alias{hyperSphere}
\alias{isotropicNormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Isotropic Distributions With or Without Noise }
\description{
  Generates a sample from isotropic distributions in \code{d} dimensions with
  \code{n}-dimensional noise added to it.
}
\usage{
hyperBall(Ns, d, n = d, sd = 0)
hyperSphere(Ns, d, n = d + 1, sd = 0)
isotropicNormal(Ns, d, n = d, sd = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ns}{ number of points. }
  \item{d}{ intrinsic dimension of the support of the distribution 
            (the manifold.) }
  \item{n}{ dimension of noise. }
  \item{sd}{ standard deviation of noise. }
}
\details{
  \code{hyperBall} draws a sample from a uniform distribution on a hyperball of 
  radius 1.
  \code{hyperSphere} draws a sample from a uniform distribution on a hypersphere
  of radius 1.
  \code{isotropicNormal} draws a sample from a isotropic normal distribution with
  identity covariance matrix.
}

\author{
  Kerstin Johnsson, Lund University
}

\examples{
datap <- hyperSphere(100, 1, 3, sd = .1)
par(mfrow = c(1, 2))
plot(datap[, 1], datap[, 2])
plot(datap[, 1], datap[, 3])    

}
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
