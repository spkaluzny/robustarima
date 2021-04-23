\name{arima.rob}
\alias{arima.rob}
\title{
  Robust Fit of a REGARIMA Model and Outliers Detection
}
\description{
Returns an object of class \code{"arima.rob"}
that represents a robust fit of a linear
regression model with ARIMA errors using a filtered tau-estimate. The error
model may have seasonal differences and one seasonal moving average
parameter. It also returns the detected outliers and level shifts.
}
\usage{
arima.rob(formula, data, contrasts=NULL, start=NULL, end=NULL,
      p=0, q=0, d=0, sd=0, freq=1, sfreq=NULL, sma=FALSE,
      max.p=NULL, auto.ar=FALSE, n.predict=20, tol=10^(-6),
      max.fcal=2000, iter=FALSE, innov.outlier=FALSE, critv=NULL, \dots)
}
\arguments{
\item{formula}{
a \code{formula} object, with the response on
the left of a \code{~} operator, and
the terms, separated by \code{+} operators,
on the right.
}
\item{data}{
a data frame or a \code{"timeSeries"} object with a data frame
in the data slot, which is used to interpret the variables named in
\code{formula}. If this
is missing, then the variables in \code{formula}
should be on the search list. Missing values are not allowed.
}
\item{contrasts}{
the same as the \code{contrasts} argument
for \code{lm} function.
}
\item{start}{
a character string which can be passed to \code{timeDate}
function to specify the starting date for the estimation.
This can only be used if the \code{data} argument is
a \code{"timeSeries"} data frame. The default is
\code{NULL}.
}
\item{end}{
a character string which can be passed to \code{timeDate}
function to specify the ending date for the estimation.
This can only be used if the \code{data} argument is
a \code{"timeSeries"} data frame. The default is
\code{NULL}.
}
\item{p}{
the autoregressive order of the errors model.
The default is \code{0}.
}
\item{q}{
the moving average order of the errors model.
The default is \code{0}.
}
\item{d}{
the number of regular differences in the ARIMA model.
It must be \code{0}, \code{1} or \code{2}.
The default is \code{0}.
}
\item{sd}{
the number of seasonal differences. It must be \code{0},
\code{1} or \code{2}. The default is \code{0}.
}
\item{freq}{
the frequency of \code{data}. The default is \code{1}.
}
\item{sfreq}{
the seasonality frequency of \code{data}. If \code{NULL},
it is set to be equal to \code{freq}. The default is
\code{NULL}.
}
\item{sma}{
logical flag: if \code{TRUE}, the errors model
includes a seasonal moving average parameter. The default is
\code{FALSE}.
}
\item{auto.ar}{
logical flag: If \code{TRUE} an AR(p) model
is selected automatically using a robust AIC criterion.
The default is \code{FALSE}.
}
\item{max.p}{
the maximum order of the autoregressive stationary model that approximates
the ARMA stationary model. If \code{NULL},
\code{max.p=max(p+q,5)}.
If \code{q=0},
then \code{max.p} is not necessary.
The default is \code{NULL}.
}
\item{n.predict}{
the maximum number of future periods for which we wish to compute the
predictions. The default is \code{20}.
}
\item{tol}{
the tolerance for convergence.
}
\item{max.fcal}{
the maximum number of function evaluations.
}
\item{iter}{
a logical flag or the number of iterations to
execute \code{arima.rob} with.
}
\item{innov.outlier}{
logical flag: if \code{TRUE}, the function
\code{arima.rob} looks for innovation outliers
in addition to additive outliers and level shifts;
otherwise, \code{arima.rob} only looks for additive outliers and level shifts.
The default is \code{FALSE}.
}
\item{critv}{
the critical value for detecting outliers.
If \code{NULL}, it assumes the
following default values:
\code{critv=3} if the length of the time series is less than 200;
\code{critv=3.5} if it is between 200 and 500, and
\code{critv=4} if it is greater than 500.
}
\item{\dots}{
extra arguments passed to or from other methods.
}
}
\value{
an object of class \code{"arima.rob"}
representing the fit and the outliers detected.
See \code{arima.rob.object} for details of
the components of the object.
}
\section{Warning}{
When either \code{d} or
\code{sd} is greater than zero, the interpretation
of the intercept in the \code{formula} is
different from its usual interpretation:
it represents the coefficient of the lowest order power of the time trend
which can be identified. For example, if
\code{d=2} and
\code{sd=0}, the
intercept represents the coefficient of the
term \code{t^2}.
}
\section{References}{
Bianco, A., Garcia Ben, M., Martinez, E., and Yohai, V. (1996).
Robust procedures for regression models with ARIMA errors.
\emph{COMPSTAT 96, Proceedings in Computational Statistics}.
Ed. Albert Prat, pages. 27-38. Physica-Verlag, Heidelberg.

Bianco, A., Garcia Ben, M., Martinez, E., and Yohai, V. (1997).
Outlier detection in regression models with ARIMA errors
using robust estimates.
mimeo.

Chang, I., Tiao, G. C., and Chen, C. (1988).
Estimation of time series parameters in the presence of outliers.
\emph{Technometrics},
30:193-204.

Martin, R. D., Samarov, A., and Vandaele, W. (1983).
Robust methods for ARIMA models. in
\emph{Applied Time Series Analysis of Economic Data},
E. Zellner, ed.

Yohai, V. Y., and Zamar, R. H. (1988).
High breakdown-point estimates of
regression by means of the minimization of an efficient scale.
\emph{Journal of the American Statistical Association},
83:406-413.
}
\seealso{
\code{\link{arima.rob.object}}.
}
\examples{
frip.rr <- arima.rob(log(frip.dat) ~ 1, p=2, d=1)
}
\keyword{ts}
\keyword{robust}
\keyword{models}
\keyword{regression}
% docclass is function
% Converted by mySd2Rd version 47442.
