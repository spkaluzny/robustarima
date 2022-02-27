\name{arima.rob.fit}
\alias{arima.rob.fit}
\title{
Robust Fitting of a REGARIMA Model
}
\description{
This is an internal function called iteratively by \code{arima.fit}.
Users are not likely to ever call it directly.
}
\usage{
arima.rob.fit(formula, data, contrasts = NULL, start = NULL, end = NULL,
     p = 0, q = 0, d = 0, sd = 0, freq = 1, sfreq = NULL, sma = FALSE,
     max.p = NULL, auto.ar = FALSE, n.predict = 20, y.cleaned = NULL,
     tol = 10^(-6), max.fcal = 2000, method = "fit", ...)
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
\item{y.cleaned}{
the cleaned response series.
}
\item{tol}{
the tolerance for convergence.
}
\item{max.fcal}{
the maximum number of function evaluations.
}
\item{method}{
a string that specifies the method.
}
\item{\dots}{
extra arguments passed to or from other methods.
}
}
\value{
an object of class \code{"arima.rob"}.
See \code{link{arima.rob.object}}.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\seealso{
\code{\link{arima.rob}}
}
\examples{
# See arima.rob for example
}
% \keyword{internal}
