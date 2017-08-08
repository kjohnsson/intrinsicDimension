\name{pcaOtpmPointwiseDimEst}
\alias{pcaOtpmPointwiseDimEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dimension Estimation With Optimally Topology Perserving Maps }
\description{
  Intrinsic dimension estimation with the metod proposed in Bruske and Sommer
  (1998). A graph called optimally topology perserving map (OTPM) is constructed
  and on this local PCA is made with the Fukunaga-Olsen criterion to determine
  which eigenvalues that are significant.
}
\usage{ pcaOtpmPointwiseDimEst(data, N, alpha = .05) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ a data set for which dimension should be estimated. }
  \item{N}{ the number of the nodes in the OTPM.}
  \item{alpha}{ the significance level for the Fukunaga-Olsen method. }
}

\value{
A \code{DimEstPointwise} object, inhereting \code{data.frame}, with two columns:
\item{dim.est }{The dimension estimate at each point.}
\item{nbr.nb}{The number of neighbouring nodes used for the dimension estimate at each point.}
}
\references{
  Bruske, J. and Sommer, G. (1998) Intrinsic dimensionality estimation with
  optimally topology perserving maps. \emph{IEEE Trans. on Pattern Anal. and Mach.
  Intell.}, \bold{20}(5), 572-575.
}
\author{
  Kerstin Johnsson, Lund University
}

\seealso{
 \code{\link{pcaLocalDimEst}}
}
\examples{
data <- hyperBall(1000, 5)
pcaOtpmPointwiseDimEst(data, 400)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
