\name{neighborhoods}
\alias{neighborhoods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Obtaining neighborhoods (local data) from a data set
}
\description{
Get a list of neighborhoods, each containing the \code{k} nearest neighbors (not including itself) to a point in the data set. 
}
\usage{
neighborhoods(data, k, indices, eps=0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data set.}
  \item{k}{ The number of neighbors in each neighborhood.}
  \item{indices}{ A vector with indices of the points in \code{data} that should be used as center for neighborhoods. }
  \item{eps}{ If non-zero, the relative error in distance allowed when finding nearest neighbors. See Details.}
}
\details{
The \code{ann} function of the package \code{yaImpute} is used for finding the \code{k} nearest neighbors. The \code{eps} parameter to \code{neighborhoods} is used in the \code{ann} function.
}
\value{
A list of neighborhoods where each item corresponds to one index in \code{indices} and each item contains a data set with \code{k} data points.
}
\author{
Kerstin Johnsson, Lund University
}
\examples{
data <- swissRoll3Sph(300, 300)
neighborhoods(data, 10, 1:10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
