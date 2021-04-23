\name{outliers}
\alias{outliers}
\title{
  Outliers Extraction for an <code>arima.rob</code> Object
}
\description{
Returns an object of class \code{"outliers"}.
}
\usage{
outliers(object, iter=NULL)
}
\arguments{
\item{object}{
an object of class \code{"arima.rob"}.
}
\item{iter}{
a number specifying from which iteration to extract the detected outliers,
if the \code{iter} argument passed to the
function \code{arima.rob} that
produced \code{object}
is non-zero. The default is set to \code{NULL}.
}
}
\value{
an object of class \code{"outliers"}.
If \code{iter} is \code{NULL},
the object contains
all the detected outliers (and level shifts).
If \code{iter} is not
\code{NULL}, the
object contains the outliers (and level shifts) detected in iteration
\code{iter}.
See \code{outliers.object} for
components of the returned object.
}
\seealso{
\code{\link{outliers.object}},
\code{\link{arima.rob}},
\code{\link{arima.rob.object}}.
}
\examples{
frip.rr <- arima.rob(log(frip.dat) ~ 1, p=2, d=1, iter=2)
frip.outliers.all <- outliers(frip.rr)
frip.outliers.2 <- outliers(frip.rr, iter=2)
}
\keyword{ts}
% docclass is function
% Converted by mySd2Rd version 47442.
