\name{arima.rob.fit}
\alias{arima.rob.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
arima.rob.fit(formula, data, contrasts = NULL, start = NULL, end = NULL, p = 0, q = 0, d = 0, sd = 0, freq = 1, sfreq = NULL, sma = FALSE, max.p = NULL, auto.ar = FALSE, n.predict = 20, y.cleaned = NULL, tol = 10^(-6), max.fcal = 2000, method = "fit", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{contrasts}{
%%     ~~Describe \code{contrasts} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{end}{
%%     ~~Describe \code{end} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
}
  \item{freq}{
%%     ~~Describe \code{freq} here~~
}
  \item{sfreq}{
%%     ~~Describe \code{sfreq} here~~
}
  \item{sma}{
%%     ~~Describe \code{sma} here~~
}
  \item{max.p}{
%%     ~~Describe \code{max.p} here~~
}
  \item{auto.ar}{
%%     ~~Describe \code{auto.ar} here~~
}
  \item{n.predict}{
%%     ~~Describe \code{n.predict} here~~
}
  \item{y.cleaned}{
%%     ~~Describe \code{y.cleaned} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{max.fcal}{
%%     ~~Describe \code{max.fcal} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, contrasts = NULL, start = NULL, end = NULL, 
    p = 0, q = 0, d = 0, sd = 0, freq = 1, sfreq = NULL, sma = FALSE, 
    max.p = NULL, auto.ar = FALSE, n.predict = 20, y.cleaned = NULL, 
    tol = 10^(-6), max.fcal = 2000, method = "fit", ...) 
{
    call = match.call()
    mm = match.call(expand.dots = FALSE)
    mm$d = mm$sd = mm$freq = mm$sfreq = mm$sma = mm$max.p = mm$auto.ar = mm$n.predict = mm$y.cleaned = mm$tol = mm$max.fcal = mm$start = mm$end = mm$p = mm$q = mm$method = mm$... = NULL
    mm$na.action = as.name("na.pass")
    mm[[1]] = as.name("model.frame")
    isTS = FALSE
    if (missing(data)) {
        if (is(tmp <- eval(formula[[2]], envir = environment(formula)), 
            "timeSeries")) {
            isTS = TRUE
            pos = positions(tmp)
            formula[[2]] = bquote(splusTimeSeries::seriesData(.(f2)), 
                list(f2 = formula[[2]]))
        }
    }
    else {
        if (is(data, "timeSeries")) {
            isTS = TRUE
            pos = positions(data)
        }
    }
    if (isTS) {
        if (!is.null(start) || !is.null(end)) {
            tmp = tssub(pos, start, end, ...)
            pos = tmp$td
            mm$subset = tmp$idx
        }
    }
    Terms = terms(formula)
    mm$formula = formula
    mm = eval(mm, environment(formula))
    y = model.extract(mm, "response")
    if (!is.null(y) && any(is.na(y))) 
        stop("There are missing data in response.")
    x = model.matrix(Terms, mm, contrasts)
    if (method == "model.matrix") {
        return(x)
    }
    if (any(is.na(x))) 
        stop("There are missing data in model matrix.")
    tmp = dim(x)
    if (!is.null(tmp)) {
        n = tmp[1]
        m = tmp[2]
    }
    else {
        n = length(y)
        m = 0
        x = 0
    }
    interc = attr(Terms, "intercept")
    auto.ar = as.numeric(auto.ar)
    sma = as.numeric(sma)
    if (auto.ar == 1 && is.null(max.p)) 
        max.p = 5
    if (auto.ar == 0 && is.null(max.p)) 
        max.p = max(p + q, 5)
    if (auto.ar == 0 && !is.null(max.p)) 
        max.p = max(p + q, max.p)
    if (auto.ar == 0 && q == 0) 
        max.p = p
    if (is.null(sfreq)) 
        sfreq = freq
    if (p < 0) 
        stop("The order of the AR model must be >= 0.")
    if (q < 0) 
        stop("The order of the MA model must be >= 0.")
    ndim1 = max.p + q + m + sma + 1
    ndim2 = max(max.p + d + sfreq * sd, q + sma * sfreq + 1)
    nw1 = ndim2 * ndim2 + 7 * (ndim2 + 1) + 2 * m + (max.p + 
        1)^2 + (max.p + 1) + d + (n + 7) * ndim1 + 8 * n + n.predict + 
        max.p + q + sfreq * sd + 1 + (n + n.predict) * ndim2
    niw1 = ndim1 + max.p + q
    nw2 = 16 * n + 9 * m + ndim1 * (7 + n) + 30 * ndim2 + 19 * 
        ndim2 * ndim2 + 4 * m * m + 16 + 5 * max.p + d + sfreq * 
        sd + 3 * q + 2 * sfreq + 2 * (n + n.predict) * ndim2 + 
        n.predict + n * m
    niw2 = max(max.p + q + ndim1, ndim2, m) + max(ndim2, 2 * 
        m) + 2 * ndim2 + 3
    work1 = double(nw1)
    work2 = double(nw2)
    iwork1 = integer(niw1)
    iwork2 = integer(niw2)
    tmpn = double(n)
    tmpnnp = double(n + n.predict)
    storage.mode(x) = "double"
    if (is.null(y.cleaned)) 
        storage.mode(y) = "double"
    else {
        if (is(y.cleaned, "timeSeries")) {
            pos = positions(y.cleaned)
            y.cleaned = seriesData(y.cleaned)
        }
        storage.mode(y.cleaned) = "double"
    }
    if (d < 0 || d > 2) 
        stop("The number of regular differences must be 0, 1 or 2.")
    if (sd < 0 || sd > 2) 
        stop("The number of seasonal differences must be 0, 1 or 2.")
    dift = d + sd
    if (dift > 0 && interc == 1) {
        aux = 1:n
        x[, 1] = (aux^dift)/(prod(1:dift) * (sfreq^sd))
    }
    one = as.double(1)
    regcoef.cov = matrix(double(1), max(m, 1), max(m, 1))
    xy = matrix(double(1), n, m + 1)
    double.eps = as.double(5.18065378653631e-318)
    double.xmin = as.double(5.02523400001724e-315)
    idx = c("x", "y", "popt", "phi", "theta", "sma", "regcoef", 
        "sigma.innov", "regcoef.cov", "innov.acf", "regresid.acf", 
        "sigma.regresid", "tuning.c", "sigma.first", "y.robust", 
        "innov", "regresid", "predict.scales", "predict.error", 
        "n.predict", "tauef", "inf", "n0")
    tt = .Fortran(F_s_regafe, x = x, y = if (is.null(y.cleaned)) y else y.cleaned, 
        n = as.integer(n), m = as.integer(m), d = as.integer(d), 
        sfreq = as.integer(sfreq), sd = as.integer(sd), max.p = as.integer(max.p), 
        auto.ar = as.integer(auto.ar), p = as.integer(p), q = as.integer(q), 
        seasonal.ma = as.integer(sma), interc = as.integer(interc), 
        popt = integer(1), phi = double(ndim2), theta = double(2 * 
            ndim1), sma = double(1), regcoef = double(max(m, 
            1)), sigma.innov = one, regcoef.cov = regcoef.cov, 
        innov.acf = tmpn, regresid.acf = tmpn, sigma.regresid = one, 
        tuning.c = one, sigma.first = one, xy, y.robust = tmpn, 
        innov = tmpnnp, regresid = tmpnnp, predict.scales = tmpnnp, 
        predict.error = tmpnnp, n.predict = as.integer(n.predict), 
        tauef = one, inf = as.integer(1), as.integer(ndim1), 
        as.integer(ndim2), work1, as.integer(nw1), iwork1, as.integer(niw1), 
        work2, as.integer(nw2), iwork2, as.integer(niw2), as.double(tol), 
        as.integer(max.fcal), double.eps, double.xmin, n0 = integer(1))[idx]
    n0 = tt$n0
    if (is.null(y.cleaned)) {
        tt$y.cleaned = NULL
    }
    else {
        tt$y.cleaned = tt$y
        tt$y = y
    }
    if (tt$inf <= 4 && tt$inf != 0) 
        tt$inf = 1
    else {
        if (tt$inf != 5) 
            warning("The last optimization procedure did not converge")
        if (tt$inf == 5) 
            warning("The last optimization procedure did not converge: ", 
                "increase max.fcal.")
        tt$inf = 0
    }
    tt$model = list(d = d, sd = sd, sfreq = sfreq, freq = freq)
    if (auto.ar == 0) {
        if (p > 0) {
            tt$model$ar = tt$phi[1:p]
            names(tt$model$ar) = paste("AR(", 1:p, ")", sep = "")
        }
        tt$phi = NULL
        if (q > 0) {
            tt$model$ma = tt$theta[1:q]
            names(tt$model$ma) = paste("MA(", 1:q, ")", sep = "")
        }
        tt$theta = NULL
        tt$popt = NULL
    }
    else {
        tt$model$ar = tt$phi[1:tt$popt]
        names(tt$model$ar) = paste("AR(", 1:tt$popt, ")", sep = "")
        tt$phi = tt$theta = tt$popt = NULL
    }
    if (sma > 0) {
        tt$model$sma = tt$sma
        names(tt$model$sma) = "SMA"
    }
    tt$sma = NULL
    if (m == 0) {
        tt$regcoef.cov = NULL
        tt$regcoef = NULL
    }
    else {
        names(tt$regcoef) = colIds(tt$x)
        tt$regcoef.cov = as.matrix(tt$regcoef.cov[1:m, 1:m])
        dimnames(tt$regcoef.cov) <- list(names(tt$regcoef), names(tt$regcoef))
    }
    tt$innov.acf = tt$innov.acf[1:(n - n0)]
    tt$regresid.acf = tt$regresid.acf[1:(n - d - sfreq * sd)]
    tt$innov[1:n0] = NA
    tt$innov = tt$innov[1:n]
    tt$predict.error[1:n0] = NA
    tt$predict.scales[1:n0] = NA
    tt$terms = Terms
    tt$assign = attr(x, "assign")
    tt$contrasts = attr(x, "contrasts")
    tt$rank = m
    tt$call = call
    if (isTS) {
        if (!is.null(tt$y.cleaned)) {
            tt$y.cleaned = timeSeries(tt$y.cleaned, positions. = pos)
        }
        tt$y = timeSeries(tt$y, positions. = pos)
        tt$innov = timeSeries(tt$innov, positions. = pos)
    }
    oldClass(tt) = "arima.rob"
    tt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
