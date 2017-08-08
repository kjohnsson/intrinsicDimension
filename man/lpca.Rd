\name{pcaLocalDimEst}
\alias{pcaLocalDimEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Local Dimension Estimation with PCA }
\description{
  Estimates local manifold dimension using the largest singular values of the covariance matrix.
}
\usage{
pcaLocalDimEst(data, ver, alphaFO = .05, alphaFan = 10, betaFan = .8, PFan = .95,
     ngap = 5, maxdim = min(dim(data)), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ a local data set for which dimension should be estimated.}
  \item{ver}{ possible values: \code{'FO'}, \code{'fan'}, \code{'maxgap'}, \code{'cal'}. \code{'cal'} is often very slow. }
  \item{alphaFO}{ only for \code{ver = 'FO'}. An eigenvalue is considered significant if it is larger than alpha times the largest eigenvalue. }
  \item{alphaFan}{ only for \code{ver = 'Fan'}. The alpha parameter (large gap threshold).}
  \item{betaFan}{ only for \code{ver = 'Fan'}. The beta parameter (total covariance threshold).}
  \item{PFan}{ only for \code{ver = 'Fan'}. Total covariance in non-noise.}
  \item{ngap}{ only for \code{ver = 'cal'}. How many of the largest gaps that should be considered. }
  \item{maxdim}{ only for \code{ver = 'cal'}. The maximal manifold dimension of the data. }
  \item{verbose}{ should information about the process be printed out? }
}
\details{
  Version \code{'FO'} is the method by Fukunaga-Olsen, version \code{'fan'} is the method by Fan et al..

  Version \code{'maxgap'} returns the position of the largest relative gap in the sequence of singular values.

  Version \code{'cal'} considers the positions of the \code{ngap} largest relative gaps in the
  sequence of singular values and generates calibration data to determine which one of them is most likely.

  All versions assume that the data is local, i.e. that it is a neighborhood taken from a larger data set,
such that the curvature and the noise within the neighborhood is relatively small. In the ideal case
(no noise, no curvature) this is equivalent to the data being uniformly distributed over a hyper ball.
}
\value{
A \code{DimEst} object with slots:
  \item{dim.est}{the dimension estimate}
  \item{gap.size}{if \code{ver} is not \code{'cal'}, the size of the gap in singular values corresponding to the estimated dimension}
  \item{likelihood}{if \code{ver} is \code{cal}, the likelihood of the estimated dimension.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Fukunaga, K. and Olsen, D. R. (1971). An algorithm for finding intrinsic dimensionality
  of data. \emph{IEEE Trans. Comput.}, \bold{c-20}(2):176-183.

Fan, M. et al. (2010). Intrinsic dimension estimation of data by principal component 
  analysis. \emph{arXiv preprint} 1002.2050.
}
\author{
  Kerstin Johnsson, Lund University
}
\seealso{
 \code{\link{pcaOtpmPointwiseDimEst}}
}
\examples{
data <- cutHyperPlane(100, 4, 10, .05)
pcaLocalDimEst(data, 'fan')
pcaLocalDimEst(data, 'FO')
pcaLocalDimEst(data, 'maxgap')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
