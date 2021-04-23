\name{outliers.object}
\alias{outliers.object}
\title{
  Outliers Objects Extracted from <code>arima.rob</code> Objects
}
\description{
These are objects of class \code{"outliers"}
extracted an object of class \code{"arima.rob"},
which contain information about the detected outliers (and level shifts).
}
\section{Generation}{
This class of objects is returned from the
\code{outliers} function.
}
\section{Methods}{
\code{print},
\code{summary}.
}
\arguments{
The following components must be included in a
legitimate \code{"outliers"} object:

\item{nout}{
the number of outliers (and level shifts) detected.
}
\item{outlier.index}{
the index of each detected outlier (or level shift).
}
\item{outlier.type}{
the type of each detected outlier (or level shift): 1=innovation outlier,
2=additive outlier, 3=level shift.
}
\item{outlier.impact}{
the size of each detected outlier (or level shift).
}
\item{outlier.t.statistics}{
the t-statistics for each detected outlier (or level shift).
}
\item{outlier.positions}{
the \code{"timeDate"} objects associated with the
detected outliers if the original data is a \code{"timeSeries"}
object.
}
\item{sigma0}{
the estimate of the innovation scale before correcting the outliers.
}
\item{sigma}{
the estimate of the innovation scale after correcting the outliers.
}
\item{ierror}{
the error indicator: if \code{ierror > 0} the
search for outliers was stopped, the
algorithm detected too many outliers.
}
}
\seealso{
\code{\link{outliers}}
\code{\link{arima.rob}},
\code{\link{arima.rob.object}}.
}
\keyword{ts}
\keyword{methods}
% docclass is function
% Converted by mySd2Rd version 47442.
