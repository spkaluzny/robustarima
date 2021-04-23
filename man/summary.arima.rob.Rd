\name{summary.arima.rob}
\alias{summary.arima.rob}
\title{
  Summary Method for arima.rob Objects
}
\description{
Returns a summary list for an \code{"arima.rob"} object.
}
\usage{
\method{summary}{arima.rob}(object, correlation = FALSE, \dots)
}
\arguments{
\item{object}{
an object of class \code{"arima.rob"}.
}
\item{correlation}{
a logical flag: if \code{TRUE}, correlation matrices of
regression coefficients and ARIMA coefficients are also produced.
The default is \code{FALSE}.
}
\item{\dots}{
extra arguments passed to or from other methods.
The summary method here ignore these arguments.
}
}
\value{
an object of class \code{"summary.arima.rob"} which
must contain the following components:

\item{ARIMA.model}{
the same list as the \code{model} component of \code{object}.
See \code{arima.rob.object} for details.
}
\item{reg.coef}{
a matrix with four columns, containing the regression coefficients, their
standard errors, the t-statistics and the corresponding p-values.
}
\item{regcoef.cov}{
the estimated covariance matrix for the regression coefficients.
}
\item{regcoef.corr}{
the estimated correlation matrix for the regression coefficients.
This is only present if \code{corr=T}.
}
\item{AR.coef}{
a matrix with four columns, containing the AR coefficients, their standard
errors, the t-statistics and the p-values.
}
\item{MA.coef}{
a matrix with four columns, containing the MA coefficients, their standard
errors, the t-statistics and the p-values.
}
\item{sMA.coef}{
an array which contains the seasonal moving average parameter, its standard
error, the t-statistic and the p-value.
}
\item{ARIMA.cov}{
the estimated covariance matrix of the ARMA coefficients.
}
\item{ARIMA.corr}{
the estimated correlation matrix of the ARMA coefficients.
This is only present if \code{corr=T}.
}
\item{n}{
the length of the time series.
}
\item{df}{
the number of degrees of freedom for the model.
}
\item{sigma}{
the estimate of the innovations scale.
}
\item{call}{
the image of the original call to \code{arima.rob}.
}
\item{outliers}{
an object of class \code{"summary.outliers"}.
}
}
\details{
This function is a method for the generic function
\code{summary} for class
\code{"arima.rob"}. It can be invoked by
calling \code{summary} for an object of the
appropriate class, or directly by calling
\code{summary.arima.rob} regardless of the
class of the object.
}
\seealso{
\code{\link{arima.rob}},
\code{\link{arima.rob.object}},
\code{\link{summary}}.
}
\examples{
frip.rr <- arima.rob(log(frip.dat) ~ 1, p=2, d=1)
summary(frip.rr)
}
\keyword{regression}
\keyword{robust}
\keyword{ts}
% docclass is function
% Converted by mySd2Rd version 47442.
