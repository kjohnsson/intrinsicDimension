\name{tp}
\alias{maxLikGlobalDimEst}
\alias{maxLikPointwiseDimEst}
\alias{maxLikLocalDimEst}
\alias{hill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dimension Estimation via Translated Poisson Distributions }
\description{
Estimates the intrinsic dimension of a data set using models of translated
Poisson distributions.
}
\usage{
maxLikGlobalDimEst(data, k, dnoise = NULL, sigma = 0, n = NULL,
        integral.approximation = 'Haro', unbiased = FALSE,
        neighborhood.based = TRUE,
        neighborhood.aggregation = 'maximum.likelihood', iterations = 5, K = 5)
maxLikPointwiseDimEst(data, k, dnoise = NULL, sigma = 0, n = NULL, indices = NULL,
             integral.approximation = 'Haro', unbiased = FALSE, iterations = 5)
maxLikLocalDimEst(data, dnoise = NULL, sigma = 0, n = NULL,
       integral.approximation = 'Haro',
       unbiased = FALSE, iterations = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data set with each row describing a data point.}
  \item{k}{the number of distances that should be used for each dimension estimation.}
  \item{dnoise}{a function or a name of a function giving the translation density. If NULL, no noise is modeled, and the estimator turns into the Hill estimator (see References). Translation densities \code{\link{dnoiseGaussH}}, \code{\link{dnoiseGaussB}} and \code{\link{dnoiseNcChi}} are provided in the package.}
  \item{sigma}{(estimated) standard deviation of the (isotropic) noise.}
  \item{n}{dimension of the noise.}
  \item{indices}{the indices of the data points for which local dimension
    estimation should be made.}
  \item{integral.approximation}{ how to approximate the integrals in eq. (5) in Haro et al. (2008). Possible values: \code{'Haro'}, \code{'guaranteed.convergence'}, \code{'iteration'}. See Details.}
  \item{unbiased}{ if \code{TRUE}, a factor \code{k-2} is used instead of the factor \code{k-1} that was used in Haro et al. (2008). This makes the estimator is unbiased in the case of data without noise or boundary.}
  \item{neighborhood.based}{ if TRUE, dimension estimation is first made for neighborhoods around each data point and final value is aggregated from this. Otherwise dimension estimation is made once, based on distances in entire data set.}
  \item{neighborhood.aggregation}{ if \code{neighborhood.based}, how should dimension estimates from different neighborhoods be combined. Possible values: \code{'maximum.liklihood'} follows Haro et al. (2008) in maximizing likelihood by using the harmonic mean, \code{'mean'} follows Levina and Bickel (2005) and takes the mean, \code{'robust'} takes the median, to remove influence from possible outliers.}
  \item{iterations}{ for \code{integral.approxmation = 'iteration'}, how many iterations should be made.}
  \item{K}{ for \code{neighborhood.based = FALSE}, how many distances for each data point should be considered when looking for the \code{k} shortest distances in the entire data set.}
}
\details{
The estimators are based on the referenced paper by Haro et al. (2008), using the assumption that there is a single manifold. The estimator in the paper is obtained using default parameters and \code{dnoise = dnoiseGaussH}.

  With \code{integral.approximation = 'Haro'} the Taylor expansion approximation of \code{r^(m-1)} that Haro et al. (2008) used are employed. With \code{integral.approximation = 'guaranteed.convergence'}, \code{r} is factored out and kept and \code{r^(m-2)} is approximated with the corresponding Taylor expansion. This guarantees convergence of the integrals. Divergence might be an issue when the noise is not sufficiently small in comparison to the smallest distances. With \code{integral.approximation = 'iteration'}, five iterations is used to determine \code{m}. 
  
  \code{maxLikLocalDimEst} assumes that the data set is local i.e. a piece of a data set cut out by a sphere with a radius such that the data set is well approxiated by a hyperplane (meaning that the curvature should be low in the local data set). See \code{\link{localIntrinsicDimension}}.
}
\value{
For \code{maxLikGlobalDimEst} and \code{maxLikLocalDimEst}, a \code{DimEst} object with one slot:
  \item{dim.est}{ the dimension estimate}
%
For \code{maxLikPointwiseDimEst}, a \code{DimEstPointwise} object, inhereting \code{data.frame}, with one slot:
  \item{dim.est}{ the dimension estimate for each data point. Row \code{i} has the local dimension estimate at point \code{data[indices[i], ]}.}
}
\references{
  Haro, G., Randall, G. and Sapiro, G. (2008) Translated Poisson Mixture Model
  for Stratification Learning. \emph{Int. J. Comput. Vis.}, \bold{80}, 358-374.
  
  Hill, B. M. (1975) A simple general approach to inference about the tail of a distribution. \emph{Ann. Stat.}, \bold{3}(5) 1163-1174.
  
  Levina, E. and Bickel., P. J. (2005) Maximum likelihood estimation of intrinsic dimension. \emph{Advances in Neural Information Processing Systems 17}, 777-784. MIT Press.
}
\author{
Kerstin Johnsson, Lund University.
}

\examples{
data <- hyperBall(1000, d = 7, n = 13, sd = 0.01)
maxLikGlobalDimEst(data, 10, dnoiseGaussH, 0.01, 13)
maxLikGlobalDimEst(data, 10, dnoiseNcChi, 0.01, 13)
maxLikGlobalDimEst(data, 10, dnoiseNcChi, 0.01, 13, neighborhood.aggregation = 'robust')
maxLikGlobalDimEst(data, 10, dnoiseNcChi, 0.01, 13,
        integral.approximation = 'guaranteed.convergence')
maxLikGlobalDimEst(data, 10, dnoiseNcChi, 0.01, 13,
        integral.approximation = 'iteration', unbiased = TRUE)

data <- hyperBall(1000, d = 7, n = 13, sd = 0.01)
maxLikGlobalDimEst(data, 500, dnoiseGaussH, 0.01, 13,
        neighborhood.based = FALSE)
maxLikGlobalDimEst(data, 500, dnoiseNcChi, 0.01, 13,
        neighborhood.based = FALSE)
maxLikGlobalDimEst(data, 500, dnoiseNcChi, 0.01, 13,
        integral.approximation = 'guaranteed.convergence',
        neighborhood.based = FALSE)
maxLikGlobalDimEst(data, 500, dnoiseNcChi, 0.01, 13,
        integral.approximation = 'iteration',
        neighborhood.based = FALSE)
        
data <- hyperBall(100, d = 7, n = 13, sd = 0.01)
maxLikPointwiseDimEst(data, 10, dnoiseGaussH, 0.01, 13, indices=1:10)
maxLikPointwiseDimEst(data, 10, dnoiseNcChi, 0.01, 13, indices=1:10)

data <- cutHyperPlane(50, d = 7, n = 13, sd = 0.01)
maxLikLocalDimEst(data, dnoiseGaussH, 0.1, 3)
maxLikLocalDimEst(data, dnoiseNcChi, 0.1, 3)
maxLikLocalDimEst(data, dnoiseNcChi, 0.1, 3,
       integral.approximation = 'guaranteed.convergence')

}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
