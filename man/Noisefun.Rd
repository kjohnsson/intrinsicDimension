\name{Noisefun}
\alias{dnoiseNcChi}
\alias{dnoiseGaussH}
\title{ Transition Functions Describing Noise }
\description{
  Transition functions \eqn{f(s|r)} describing the shift in lengths of vectors
  when Gaussian noise is added. Given a length \eqn{r}, \eqn{f(s|r)} is the
  probability density for the length after noise is added to one endpoint.
}
\usage{
dnoiseNcChi(r, s, sigma, k)
dnoiseGaussH(r, s, sigma, k)
}
\arguments{
  \item{r}{ length or vector of lengths of original vector. }
  \item{s}{ length or vector of lengths of perturbed vector. }
  \item{sigma}{ noise standard deviation. }
  \item{k}{ dimension of noise. }
}
\details{
  \code{dnoiseNcChi} is the true transition function density when the noise
  is Gaussian, the other transition functions are approximations of this.
  \code{dnoiseGaussH} is the Gaussian approximation used in Haro et al.
  
  If Gaussian noise is added to both endpoints of the vector, \code{sigma}
  should be replaced by 
  \code{sqrt(2)*sigma}.
}
\value{
  Vector of probability densities.
}
\references{
  Haro, G., Randall, G. and Sapiro, G. (2008) Translated Poisson Mixture Model
  for Stratification Learning. \emph{Int. J. Comput. Vis.}, \bold{80}, 358-374.
}
\author{
  Kerstin Johnsson, Lund University
}
\note{
  Only \code{r} or \code{s} can be a vector.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{maxLikPointwiseDimEst}}, \code{\link{maxLikGlobalDimEst}}, \code{\link{maxLikLocalDimEst}}
}
\examples{
# High SNR, high-dim noise
sigma <- 0.05
x <- seq(0, 1.5, length.out = 200)
y <- dnoiseNcChi(x, s = .5, sigma, k = 20)
plot(x, y, type = 'l', main = 'Noise dim = 20')
y2 <- dnoiseGaussH(x, s = .5, sigma, k = 20)
lines(x, y2, lty = 2)

# Low SNR
par(mfrow = c(2, 3))
sigma <- 0.2
x <- seq(0, 1.5, length.out = 200)
y <- dnoiseNcChi(x, s = .5, sigma, k = 4)
plot(x, y, type = 'l', main = 'Noise approximations')
y2 <- dnoiseGaussH(x, s = .5, sigma, k)
lines(x, y2, lty = 2)

# High SNR, low-dim noise
sigma <- 0.05
x <- seq(0, 1.5, length.out = 200)
y <- dnoiseNcChi(x, s = .5, sigma, k = 4)
plot(x, y, type = 'l', main = 'Noise dim = 4')
y2 <- dnoiseGaussH(x, s = .5, sigma, k)
lines(x, y2, lty = 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
