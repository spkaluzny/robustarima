\name{predict.arima.rob}
\alias{predict.arima.rob}
\title{
  Use predict on an arima.rob Object
}
\description{
Predicts from a fitted \code{"arima.rob"} object.
}
\usage{
\method{predict}{arima.rob}(object, n.predict=1, newdata=NULL,
                  olddata=NULL, se.fit=FALSE, \dots)
}
\arguments{
\item{object}{
an object of class \code{"arima.rob"}.
}
\item{n.predict}{
the number of predictions to be returned.
}
\item{newdata}{
a data frame containing the future values of exogenous variables, if
any, at which predictions are required.
}
\item{olddata}{
a data frame containing the original data used to fit
\code{object}. This is only required if \code{tslag}
is used to create distributed lags of exogenous variables in the
original call that generated \code{object}.
}
\item{se.fit}{
logical flag: if \code{TRUE},
pointwise standard errors are computed along with the predictions.
}
\item{\dots}{
extra arguments passed to or from other methods.
The predict method here ignore these arguments.
}
}
\value{
a list containing the following components:
\item{values}{
the predicted values.
}
\item{std.err}{
the pointwise standard errors for the predictions (if \code{se.fit=TRUE}).
% an object of class \code{"forecast"}. See \code{forecast.object} for details.
}
}
\details{
This function is a method for the generic function \code{predict} for class
\code{"arima.rob"}.
It can be invoked by calling \code{predict} for an object of the
appropriate class,
or directly by calling \code{predict.arima.rob} regardless of the
class of the object.
}
\seealso{
\code{\link{arima.rob}}
}
\examples{
import.rr <- arima.rob(import ~ taxes-1, data=import.dat, p=2, d=1)
import.hat <- predict(import.rr, 5, newdata=newtaxes.dat, se=TRUE)
}
\keyword{ts}
% docclass is function
% Converted by mySd2Rd version 47442.
