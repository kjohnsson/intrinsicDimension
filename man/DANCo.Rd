\name{dancoDimEst}
\alias{dancoDimEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Dimension Estimation With the DANCo and MIND Methods
}
\description{
  Intrinsic dimension estimation with the DANCo (Ceruti et al. 2012),
  MIND_MLi and MIND_MLk (Rozza et al. 2012) methods.
}
\usage{
dancoDimEst(data, k, D, ver = "DANCo", calibration.data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ a data set for which the intrinsic dimension is estimated. }
  \item{k}{ neighborhood parameter. }
  \item{D}{ maximal dimension. }
  \item{ver}{ possible values: 'DANCo', 'MIND_MLi', 'MIND_MLk'. }
  \item{calibration.data}{ precomputed calibration data. }
}
\details{
If \code{cal = NULL} or the \code{cal$maxdim < D} new calibration data will be computed as needed.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{DimEst} object with slots:
  \item{dim.est}{the intrinsic dimension estimate.}
  \item{kl.divergence}{the KL divergence between data and reference data for the estimated dimension (if ver == 'DANCo').}
  \item{calibration.data}{calibration data that can be reused when applying DANCo to data sets
  of the same size with the same neighborhood parameter k.}
%% ...
}
\references{
  Ceruti, C. et al. (2012) DANCo: Dimensionality from Angle and Norm Concentration.
  \emph{arXiv preprint} 1206.3881.
  
  Rozza, A et al. (2012) Novel high intrinsic dimensionality estimators. \emph{Machine learning}
  \bold{89}, 37-65.
}
\author{
  Kerstin Johnsson, Lund University
}
\examples{
data <- hyperBall(50, 10)
res <- dancoDimEst(data, 8, 20)
print(res)

## Reusing calibration data
data2 <- hyperBall(50, 5)
dancoDimEst(data2, 8, 20, calibration.data=res$calibration.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
