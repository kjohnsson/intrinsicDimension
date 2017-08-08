\name{knnDimEst}
\alias{knnDimEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dimension Estimation from kNN Distances }
\description{
  Estimates the intrinsic dimension of a data set using weighted average kNN
  distances.
}
\usage{
knnDimEst(data, k, ps, M, gamma = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data set with each row describing a data point.}
  \item{k}{number of distances to neighbors used at a time.}
  \item{ps}{vector with sample sizes; each sample size has to be larger than
            k and smaller than \code{nrow(data)}.}
  \item{M}{number of bootstrap samples for each sample size. }
  \item{gamma}{weighting constant.}
}
\details{
  This is a somewhat simplified version of the kNN dimension estimation method
  described by Carter et al. (2010), the difference being that block 
  bootstrapping is not used.
}
\value{ A \code{DimEst} object with slots:
  \item{dim.est}{the intrinsic dimension estimate (integer).}
  \item{residual}{the residual, see Carter et al. (2010).} 
}
\references{
  Carter, K.M., Raich, R. and Hero, A.O. (2010) On local intrinsic dimension 
  estimation and its applications. \emph{IEEE Trans. on Sig. Proc.}, 
  \bold{58}(2), 650-663.
}
\author{
  Kerstin Johnsson, Lund University.
}
\examples{
N <- 50
data <- hyperBall(N, 5)

k <- 2
ps <- seq(max(k + 1, round(N/2)), N - 1, by = 3)
knnDimEst(data, k, ps, M = 10, gamma = 2)

}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
