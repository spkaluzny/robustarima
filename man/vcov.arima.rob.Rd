\name{vcov.arima.rob}
\alias{vcov.arima.rob}
\title{
  Use <code>vcov()</code> on an <code>arima.rob</code> Model Object
}
\description{
This is a method for the function \code{vcov()} for objects
inheriting from class \code{"arima.rob"}.
See \code{vcov} for the general behavior of this
function and for the interpretation of \code{object}.
}
\usage{
\method{vcov}{arima.rob}(object, nterm=1000, \dots)
}
\arguments{
\item{object}{
an object inheriting from class \code{"arima.rob"}.
}
\item{nterm}{
an interger specifying the number of terms.
}
\item{\dots}{
extra arguments passed to or from other methods.
The method here ignore these arguments.
}
}
\value{
the variance-covariance matrix of fitted arima.rob model returned
from a call \code{arima.rob}.
}
\seealso{
\code{\link{arima.rob.object}},
\code{\link{vcov}}.
}
\keyword{methods}
\keyword{models}
% docclass is function
% Converted by mySd2Rd version 47442.
