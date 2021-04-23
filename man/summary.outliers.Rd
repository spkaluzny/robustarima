\name{summary.outliers}
\alias{summary.outliers}
\title{
  Summary Method for outliers Objects
}
\description{
Returns an object of class \code{"summary.outliers"}
containing a summary list
for an \code{"outliers"} object.
A null value will be returned if printing is invoked.
}
\usage{
\method{summary}{outliers}(object, \dots)
}
\arguments{
\item{object}{
an object of class \code{"outliers"}.
}
\item{\dots}{
extra arguments passed to or from other methods.
The summary method here ignore these arguments.
}

}
\value{
a list is returned with the following components:

\item{nout}{
number of outliers (and level shifts) detected.
}
\item{outliers.table}{
the index (or \code{"timeDate"} information if the
argument \code{data} passed to
\code{arima.rob} is of
class \code{"timeSeries"}), type,
impact and t value for each of the detected outlier
(or level shift).
}
\item{sigma0}{
the estimate of the innovation scale before correcting outliers.
}
\item{sigma}{
the estimate of the innovation scale after correcting outliers.
}
}
\details{
This function is a method for the generic function
\code{summary} for class
\code{"outliers"}. It can be invoked by
calling \code{summary} for an object of the
appropriate class, or directly by calling
\code{summary.outliers} regardless of
the class of the object.
}
\seealso{
\code{\link{outliers}},
\code{\link{outliers.object}},
\code{\link{summary}}.
}
\examples{
frip.rr <- arima.rob(log(frip.dat) ~ 1, p=2, d=1)
frip.outliers <- outliers(frip.rr)
summary(frip.outliers)
}
\keyword{ts}
% docclass is function
% Converted by mySd2Rd version 47442.
