\name{oblongNormal}
\alias{oblongNormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Oblong Normal Distribution }
\description{
  Generates a sample from a certain anisotropic normal distribution centered
  around the origin.
}
\usage{ oblongNormal(Ns, n) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ns}{ number of data points. }
  \item{n}{ dimension of the distribution (and the data points). }
}
\details{
  In the first half of the dimensions (rounded down if \code{n} is odd)
  the standard deviation is 1 and in the rest the standard deviation is 0.25 . 
}
\value{
A \code{Ns} by \code{n} matrix.
}
\author{
  Kerstin Johnsson, Lund University
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
datap <- oblongNormal(100, 10)
par(mfrow = c(1, 2))
plot(datap[, 1], datap[, 2])
plot(datap[, 1], datap[, 6])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
