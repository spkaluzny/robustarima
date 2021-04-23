\name{arima.rob.object}
\alias{arima.rob.object}
\title{
  Robust REGARIMA Model and Outliers Detection Objects
}
\description{
These are objects of class \code{"arima.rob"} which represent
the robust fit of a
regression model with ARIMA errors. It also contains information about the
detected outliers.
}
\section{Generation}{
This class of objects is returned from the
\code{arima.rob} function.
}
\section{Methods}{
\code{coef}, \code{formula},
\code{outliers},
\code{predict}, \code{print}, \code{summary}.
}
\arguments{
The following components must be included in a legitimate
\code{"arima.rob"} object:

\item{x}{
the model matrix.
}
\item{y}{
the response variable.
}
\item{model}{
a list with the following named components: \code{"freq"}
which is the frequency of the original data, \code{"sfreq"}
which is the seasonal frequency of the original data,
\code{"d"}
which is the number of regular differences,
\code{"sd"}
which is the number of seasonal differences,
\code{"ar"} which is the estimated AR coefficients,
\code{"ma"} which is the estimated MA coefficients,
\code{"sma"} which is the seasonal MA coefficient if
estimated.
}
\item{regcoef}{
the estimates of regression coefficients.
}
\item{regcoef.cov}{
the estimated covariance matrix of the regression coefficients.
}
\item{innov}{
the estimated innovations.
}
\item{innov.acf}{
a series whose autocorrelations or partial autocorrelations are the robust
estimates of the innovation autocorrelations or partial autocorrelations.
}
\item{regresid}{
the estimated regression residuals cleaned of additive outliers by the
robust filter.
}
\item{regresid.acf}{
a series whose autocorrelations or partial autocorrelations are the robust
estimates of the autocorrelations or partial autocorrelations of the
differenced regression residuals.
}
\item{sigma.innov}{
a robust estimate of the innovation scale.
}
\item{sigma.regresid}{
an estimate of the scale of the differenced regression residuals.
}
\item{sigma.first}{
the first estimate of the innovation scale based only on the scale of the
differenced model and the ARMA parameters.
}
\item{tuning.c}{
the bandwidth of the robust filter.
}
\item{y.robust}{
the response series cleaned of outliers by the robust filter.
}
\item{y.cleaned}{
the response series cleaned of additive outliers and level shifts after
the outliers detection procedure.
}
\item{predict.error}{
the fitted and predicted regression errors.
}
\item{predict.scales}{
the standard deviations of the fitted and predicted regression errors.
}
\item{n.predict}{
the number of predicted observations, which is equal to the
\code{n.predict} argument
passed to the \code{arima.rob} function that produced
the \code{"arima.rob"} object.
}
\item{tauef}{
the inverse of the estimated efficiency factor of the tau-estimate with
respect to the LS-estimate.
}
\item{inf}{
information about the outcome of the last optimization procedure: \code{inf=1}
indicates that the procedure converged, and \code{inf=0} that the procedure did
not converge.
}
\item{innov.outlier}{
logical flag, the same as the \code{innov.outlier}
argument passed to the \code{arima.rob} function
that produced the \code{"arima.rob"} object.
}
\item{outliers}{
an object of class \code{"outliers"},
which contains all the detected outliers
(and level shifts).
}
\item{outliers.iter}{
optionally a list of objects of class \code{"outliers"},
if the \code{iter} argument
passed to the \code{arima.rob} function that
produced the \code{"arima.rob"} object is
non-zero.
}
\item{n0}{
the number of missing innovations at the beginning.
}
\item{call}{
an image of the call that produced the object, but with the arguments all
named and with the actual formula included as
the \code{formula} argument.
}
\item{assign}{
the same as the \code{assign} component of
an \code{"lm"} object.
}
\item{contrasts}{
the same as the \code{contrasts} component of
an \code{"lm"} object.
}
\item{terms}{
the same as the \code{terms} component of
an \code{"lm"} object.
}
\item{rank}{
the same as the \code{rank} component of
an \code{"lm"} object.
}
}
\seealso{
\code{\link{arima.rob}},
\code{\link{outliers}},
\code{\link{outliers.object}}.
}
\keyword{ts}
\keyword{robust}
\keyword{regression}
\keyword{methods}
% docclass is function
% Converted by mySd2Rd version 47442.
