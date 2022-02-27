arima.invoper <-
function(phi, ep=0.0001, nterm=1000) 
{ 
  p = length(phi) 
  co = c(1, -phi) 
  roots = rep(1,p)/polyroot(co)
  mroots = Mod(roots) 
  M = matrix(roots, p, p, byrow=TRUE)
  for(i in 1:p) 
    M[i,] = M[i,]^i 
  coef = rep(0, p+1) 
  coef[1] = 1 
  for (i in 2:(p+1)) { 
    coef[i] = 0 
    for (j in 1:(i-1)) {
      coef[i] = coef[i] + coef[i-j]*phi[j]
    }
  } 
  const = solve(M, coef[2:(p+1)]) 
  M1 = max(Mod(const)) 
  M2 = max(mroots) 
  N  = min(nterm, log((ep*(1-M2))/(p*M1))/log(M2)) 
  coef1 = rep(0, N) 
  coef1[1] = 1 
  for (i in 2:N) 
    coef1[i] = sum(const*(roots^(i-1))) 
  coef1 = Re(coef1) 
  coef1
}
arima.rob <-
function(formula, data, contrasts=NULL, start=NULL, end=NULL, 
         p=0, q=0, d=0, sd=0, freq=1, sfreq=NULL, sma=FALSE, max.p=NULL, 
         auto.ar=FALSE, n.predict=20, tol=10^(-6), max.fcal=2000, 
         iter=FALSE, innov.outlier=FALSE, critv=NULL, ...) 
{
  call = match.call()
  m = match.call(expand.dots=TRUE)
  m$iter = m$innov.outlier = m$critv = NULL
  m[[1]] = quote(robustarima::arima.rob.fit)
  reg.rr = if (missing(data)) eval(m, envir=environment(formula)) else eval(m, envir=data, enclos=environment(formula))
  reg.rr$outliers = outliers.rr(reg.rr, innov.outlier, critv)
  y.cleaned = outliers.clean(reg.rr$outliers, reg.rr$y)
  if (iter) {
    outliers.iter = list()
    outliers.iter[[1]] = reg.rr$outliers
    if (is.logical(iter)) {
      iter = 1
      while (reg.rr$outliers$nout > 0) {
        iter = iter + 1
        m$y.cleaned = y.cleaned
        reg.rr = eval(m)
        reg.rr$outliers = outliers.rr(reg.rr, innov.outlier, critv)
        outliers.iter[[iter]] = reg.rr$outliers
        y.cleaned = outliers.clean(reg.rr$outliers, y.cleaned)
      }
      outliers.iter[[iter]] = NULL
    }
    else {
      for (i in 2:iter) {
        m$y.cleaned = y.cleaned
        reg.rr = eval(m)
        reg.rr$outliers = outliers.rr(reg.rr, innov.outlier, critv)
        outliers.iter[[i]] = reg.rr$outliers
        y.cleaned = outliers.clean(reg.rr$outliers, y.cleaned)
      }
    }
    reg.rr$outliers = outliers.combine(outliers.iter)
    reg.rr$outliers.iter = outliers.iter
  }
  reg.rr$y.cleaned = y.cleaned
  reg.rr$call = call
  reg.rr$innov.outlier = innov.outlier
  reg.rr
}
arima.rob.fit <-
function(formula, data, contrasts=NULL, start=NULL, end=NULL, 
         p=0, q=0, d=0, sd=0, freq=1, sfreq=NULL, sma=FALSE, max.p=NULL, 
         auto.ar=FALSE, n.predict=20, y.cleaned=NULL,
         tol=10^(-6), max.fcal=2000, method="fit", ...)
{ 
##
## 1. update formula for ARMA terms and exogenous variables.
##
  call = match.call()
  mm   = match.call(expand.dots=FALSE) 
  mm$d = mm$sd = mm$freq = mm$sfreq = mm$sma = mm$max.p = 
         mm$auto.ar = mm$n.predict = mm$y.cleaned = mm$tol = 
         mm$max.fcal = mm$start = mm$end = mm$p = mm$q = 
         mm$method = mm$... = NULL 
  mm$na.action = as.name("na.pass")
  mm[[1]] = as.name("model.frame")
  isTS = FALSE
  if (missing(data)) {
    if (is(tmp <- eval(formula[[2]], envir=environment(formula)), "timeSeries")) {
      isTS = TRUE
      pos = positions(tmp)
      formula[[2]] = bquote(splusTimeSeries::seriesData(.(f2)), list(f2=formula[[2]]))
    }
  }
  else {
    if (is(data, "timeSeries")) {
      isTS = TRUE
      pos = positions(data)
    }
  }
##
## 2. update subset for timeSeries objects.
##
  if (isTS) { 
    if (!is.null(start) || !is.null(end)) {
      tmp = tssub(pos, start, end, ...)
      pos = tmp$td
      mm$subset = tmp$idx
    }
  }
##
## 3. evaluate model frame and extract Y and X
##
  Terms = terms(formula)
  mm$formula = formula
  mm = eval(mm,  environment(formula))
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
  if(!is.null(tmp)) { 
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
  if(auto.ar == 1 && is.null(max.p)) 
    max.p = 5 
  if(auto.ar == 0 && is.null(max.p)) 
    max.p = max(p+q, 5) 
  if(auto.ar == 0 && !is.null(max.p)) 
    max.p = max(p+q, max.p) 
  if(auto.ar == 0 && q == 0) 
    max.p = p 
  if(is.null(sfreq)) 
    sfreq = freq 
  if (p < 0) 
    stop("The order of the AR model must be >= 0.") 
  if (q < 0) 
    stop("The order of the MA model must be >= 0.") 
  ndim1 = max.p+q+m+sma+1 
  ndim2 = max(max.p+d+sfreq*sd, q+sma*sfreq+1) 
  nw1  = ndim2*ndim2 + 7*(ndim2+1) + 2*m + (max.p+1)^2 +
         (max.p+1) + d + (n+7)*ndim1 + 8*n + n.predict + max.p + q + 
         sfreq*sd + 1 + (n+n.predict)*ndim2 
  niw1 = ndim1 + max.p + q 
  nw2  = 16*n + 9*m + ndim1*(7+n) + 30*ndim2 + 19*ndim2*ndim2 + 
         4*m*m + 16 + 5*max.p + d + sfreq*sd + 
         3*q + 2*sfreq + 2*(n+n.predict)*ndim2 + n.predict + n*m 
  niw2 = max(max.p+q+ndim1,ndim2,m) + max(ndim2,2*m) + 2*ndim2 + 3 
  work1  = double(nw1) 
  work2  = double(nw2) 
  iwork1 = integer(niw1) 
  iwork2 = integer(niw2)
  tmpn   = double(n)
  tmpnnp = double(n+n.predict)
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
    x[,1] = (aux^dift)/(prod(1:dift)*(sfreq^sd))
  }
  one = as.double(1)
  regcoef.cov = matrix(double(1), max(m,1), max(m,1)) 
  xy = matrix(double(1), n, m+1) 
  double.eps  = as.double(5.180654e-318) 
  double.xmin = as.double(5.025234e-315)
  idx = c("x", "y", "popt", "phi", "theta", "sma", "regcoef", 
          "sigma.innov", "regcoef.cov", "innov.acf", "regresid.acf",
          "sigma.regresid", "tuning.c", "sigma.first",
          "y.robust", "innov", "regresid", "predict.scales",
          "predict.error", "n.predict", "tauef", "inf", "n0")
  tt = .Fortran(F_s_regafe,
                x = x,
                y = if (is.null(y.cleaned)) y else y.cleaned,
                n = as.integer(n),
                m = as.integer(m),
                d = as.integer(d),
                sfreq = as.integer(sfreq),
                sd = as.integer(sd),
                max.p = as.integer(max.p),
                auto.ar = as.integer(auto.ar),
                p = as.integer(p),
                q = as.integer(q),
                seasonal.ma = as.integer(sma),
                interc = as.integer(interc),
                popt = integer(1),
                phi = double(ndim2),
                theta = double(2*ndim1),
                sma = double(1),
                regcoef = double(max(m,1)),
                sigma.innov = one, 
                regcoef.cov = regcoef.cov,
                innov.acf = tmpn,
                regresid.acf = tmpn,
                sigma.regresid = one,
                tuning.c = one,
                sigma.first = one,
                xy,
                y.robust = tmpn,
                innov = tmpnnp,
                regresid = tmpnnp,
                predict.scales = tmpnnp,
                predict.error = tmpnnp,
                n.predict = as.integer(n.predict),
                tauef = one,
                inf = as.integer(1),
                as.integer(ndim1),
                as.integer(ndim2),
                work1,
                as.integer(nw1),
                iwork1,
                as.integer(niw1),
                work2,
                as.integer(nw2),
                iwork2,
                as.integer(niw2),
                as.double(tol),
                as.integer(max.fcal),
                double.eps,
                double.xmin,
                n0 = integer(1))[idx]
  n0 = tt$n0
  if (is.null(y.cleaned)) {
    tt$y.cleaned = NULL
  }
  else {
    tt$y.cleaned = tt$y
    tt$y = y
  }
  if(tt$inf <= 4 && tt$inf != 0) 
    tt$inf = 1
  else {    
    if(tt$inf != 5)
      warning("The last optimization procedure did not converge")
    if(tt$inf == 5)
      warning("The last optimization procedure did not converge: ",
              "increase max.fcal.")
    tt$inf = 0
  }
##
##  n0 = if (auto.ar == 0) max(p+d+sd*seasonal.freq, q+1)
##       else           max(popt+d+sd*seasonal.freq, q+1) 
##
  tt$model = list(d=d, sd=sd, sfreq=sfreq, freq=freq)
  if(auto.ar == 0) {
    if (p > 0) {
      tt$model$ar = tt$phi[1:p]
      names(tt$model$ar) = paste("AR(", 1:p, ")", sep="")
    }
    tt$phi = NULL
    if (q > 0) { 
      tt$model$ma = tt$theta[1:q]
      names(tt$model$ma) = paste("MA(", 1:q, ")", sep="")
    }
    tt$theta = NULL
    tt$popt = NULL
  }
  else { 
    tt$model$ar = tt$phi[1:tt$popt] 
    names(tt$model$ar) = paste("AR(", 1:tt$popt, ")", sep="")
    tt$phi = tt$theta = tt$popt = NULL 
  } 
  if(sma > 0) {
    tt$model$sma = tt$sma
    names(tt$model$sma) = "SMA"
  }
  tt$sma = NULL 
  if(m == 0) { 
    tt$regcoef.cov = NULL 
    tt$regcoef = NULL
  } 
  else {
    names(tt$regcoef) = colIds(tt$x)
    tt$regcoef.cov = as.matrix(tt$regcoef.cov[1:m,1:m])
    dimnames(tt$regcoef.cov) <- list(names(tt$regcoef),
                                     names(tt$regcoef)) 
  }
  tt$innov.acf = tt$innov.acf[1:(n-n0)] 
  tt$regresid.acf = tt$regresid.acf[1:(n-d-sfreq*sd)] 
  tt$innov[1:n0] = NA 
  tt$innov = tt$innov[1:n] 
  tt$predict.error[1:n0]  = NA 
  tt$predict.scales[1:n0] = NA 
  tt$terms = Terms 
  tt$assign = attr(x, "assign") 
  tt$contrasts = attr(x, "contrasts") 
  tt$rank = m 
  tt$call = call
  if (isTS) {
    if (!is.null(tt$y.cleaned)) {
      tt$y.cleaned = timeSeries(tt$y.cleaned, positions.=pos)
    }
    tt$y = timeSeries(tt$y, positions.=pos)
    tt$innov = timeSeries(tt$innov, positions.=pos)
#   tt$regresid = timeSeries(tt$regresid, pos=pos)
  }
  oldClass(tt) = "arima.rob"
  tt 
}
coef.arima.rob <-
function(object, ...) 
{ 
  cf <- object$regcoef
  phi <- object$model$ar
  theta <- object$model$ma
  sma <- object$model$sma 
  c(cf, phi, theta, sma)
}
predict.arima.rob <-
function(object, n.predict=1, newdata=NULL, olddata=NULL, se.fit=FALSE, ...)
{ 
  n.pred0 <- object$n.predict
  if (n.predict > n.pred0) {
    warning("Only", n.pred0," predictions are available.")
    n.predict <- n.pred0
  }
  m <- length(object$regcoef) 
  if (m > 0) {
    coefs <- coef(object)[1:m]
  }
  int <- attr(object$terms, "intercept") 	
  is.x <- (m>0 && !(m==1 && int==1)) 
  if (is.x) {
    tmp <- as.character(unlist(object$call))
    if (match("tslag", tmp, nomatch=0) != 0) {
      if (is.null(olddata)) {
        stop("olddata must be specified if tslag() is used.")
      }
    }
    if (is.null(newdata)) {
      stop("newdata must be specified if regressors are used.")
    }
    if (is(newdata, "timeSeries")) {
      isNewTS <- TRUE
      new.pos <- positions(newdata)
      newdata <- seriesData(newdata)
    }
    else {
      isNewTS <- FALSE
    }
    if (is.null(dim(newdata))) {
      stop("newdata must be a data frame.")
    }
    if (is.null(colIds(newdata))) {
      stop("newdata must have column names matching the original data.")
    }
    if (!is.data.frame(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    if (nrow(newdata) < n.predict) {
      stop("Number of rows in newdata must be >= n.predict.") 
    }
    if (is.null(olddata)) {
      oc <- object$call
      oc[[1]] <- as.name("arima.rob.fit")
      oc$formula <- delete.response(object$terms)
      oc$method <- "model.matrix"
      oc$data <- newdata
      newdata <- eval(oc)
    }
    else {
      n.old <- nrow(olddata)
      n.new <- nrow(newdata)
      olddata <- rbind(olddata, matrix(1, n.new, ncol(olddata)))
      olddata[n.old+1:n.new, names(newdata)] <- newdata
      oc <- object$call
      oc[[1]] <- as.name("arima.rob.fit")
      oc$method <- "model.matrix"
      oc$na.rm <- TRUE
      oc$data <- olddata
      newdata <- eval(oc)
      n.old <- nrow(newdata)
      newdata <- newdata[(n.old-n.new+1):n.old, , drop=FALSE]
    }
    newdata <- newdata[1:n.predict, , drop=FALSE]
#   build terms
    asgn <- object$assign
    if(se.fit) { 
      pred <- .Build.terms(newdata, coefs, object$regcoef.cov, 
                          asgn, collapse = TRUE)
      names(pred) <- c("values", "std.err")
    } 
    else {
      pred <- .Build.terms(newdata, coefs, NULL, asgn, collapse=TRUE)  
    }
  }
  else {
    isNewTS <- FALSE
  }
  n <- length(object$innov)
  fit1 <- object$predict.error[(n+1):(n+n.predict)] 
  se1 <- object$predict.scales[(n+1):(n+n.predict)]
  if(se.fit) { 
    if(is.x) {
      pred$values  <- pred$values + fit1 
      pred$std.err <- sqrt((pred$std.err)^2 + se1^2)
    } 
    else {  
      pred <- list(values=fit1, std.err=se1) 
      if(m==1) { 
        pred$values  <- pred$values + coefs[1] 
        pred$std.err <- sqrt((pred$std.err)^2 + 
                        as.vector(object$regcoef.cov)) 
      } 
    }
  } 
  else {		
    if(is.x) {
      pred <- list(values=pred+fit1) 
    } 
    else {    			
      pred <- list(values=fit1) 
      if(m==1) pred$values <- pred$values + coefs[1] 
    } 
  } 
  oldClass(pred) <- "forecast"
  pred
}
print.arima.rob <-
function(x, digits=4, ...)  
{ 
  if(!is.null(cl <- x$call)) { 
    cat("\nCall:\n") 
    dput(cl) 
  }
  m <- length(x$regcoef) 
  if (m>0) { 
    cat("\nRegression Coefficients:\n") 
    print(format(round(x$regcoef, digits=digits)), quote=FALSE, ...) 
  } 
  tmp <- x$model$ar 
  if ((p <- length(tmp)) > 0) {  
    cat("\nAR Coefficients:\n") 
    print(format(round(tmp, digits=digits)), quote=FALSE, ...) 
  } 
  tmp <- x$model$ma 
  if ((q <- length(tmp)) > 0) {  
    cat("\nMA Coefficients:\n") 
    print(format(round(tmp, digits=digits)), quote=FALSE, ...) 
  } 
  tmp <- x$model$sma
  if ((indsma <- length(tmp)) > 0) {
    cat("\nSeasonal MA Coefficient:", 
        format(round(tmp, digits=digits)), "\n") 
  }
  n <- length(x$y)
  df <- n-m-x$model$d-p-q-x$model$sfreq*x$model$sd-indsma 
  cat("\nDegrees of freedom:", n, "total;", df,"residual\n") 
  tmp <- x$sigma.innov
  cat("Innovations standard deviation:", 
      format(round(tmp, digits=digits)), "\n")
  if (!is.null(x$outliers)) 
    print.outliers(x$outliers) 
  invisible(x) 
}
print.summary.arima.rob <-
function(x, digits=4, ...) 
{ 
  if(!is.null(cl <- x$call)) { 
    cat("\nCall:\n") 
    dput(cl) 
  } 
  dif <- x$ARIMA.model$d
  p <- length(x$ARIMA.model$ar)
  q <- length(x$ARIMA.model$ma)
  sma <- length(x$ARIMA.model$sma)
  sfreq <- x$ARIMA.model$sfreq
  sdif <- x$ARIMA.model$sd
#
  m <- length(x$reg.coef)
  if (m > 0) { 
    cat("\nRegression model:\n ") 
    print(x$call$formula) 
  } 
  if (dif>0 || p>0 || q>0 || sdif>0 || sma>0) { 
    cat("\nARIMA model:\n") 
    cat("Ordinary differences:",dif,"; AR order:",p, "; MA order:",q, "\n")  
    if (sdif>0 || sma>0) { 
      cat("Seasonal differences:",sdif,"; Seasonal period:",  
          sfreq, "; Seasonal MA:",sma,"\n") 
    } 
  } 
  if (m>0) { 
    cat("\nRegression Coefficients:\n") 
    print(format(round(x$reg.coef, digits=digits)), quote=FALSE, ...) 
  } 
  if (p>0) {  
    cat("\nAR Coefficients:\n") 
    print(format(round(x$AR.coef, digits=digits)), quote=FALSE, ...) 
  } 
  if (q>0) {  
    cat("\nMA Coefficients:\n") 
    print(format(round(x$MA.coef, digits=digits)), quote=FALSE, ...) 
  } 
  if (sma) { 
    cat("\nSeasonal MA Coefficient:\n") 
    print(format(round(x$sMA.coef, digits=digits)), quote=FALSE, ...) 
  } 
  cat("\nDegrees of freedom:", x$n, "total;", x$df, "residual\n") 
  cat("\nInnovations standard deviation:", 
      format(x$sigma, digits=digits), "\n") 
  if (!is.null(x$regcoef.corr)) { 
    cat("\nCorrelation Matrix of Regression Coefficients:\n")  
    print(format(round(x$regcoef.corr, digits=digits)), quote=FALSE, ...) 
  } 
  if (!is.null(x$ARIMA.corr)) {
    cat("\nCorrelation Matrix of ARIMA Coefficients:\n")  
    print(format(round(x$ARIMA.corr, digits=digits)), quote=FALSE, ...) 
  } 
  print.summary.outliers(x$outliers)
  invisible(x) 
}
summary.arima.rob <-
function(object, correlation=FALSE, ...)
{ 
  n <- length(object$y)
  m <- length(object$regcoef)
  dif <- object$model$d
  phi <- object$model$ar 
  p <- length(phi)
  theta <- object$model$ma
  q <- length(theta)
  sfreq <- object$model$sfreq
  sdif <- object$model$sd
  seasonal.theta <- object$model$sma
  sma <- length(seasonal.theta)
  df <- n-m-dif-p-q-sfreq*sdif-sma 
# 
  if (m > 0) { 
    coefi <- array(object$regcoef, c(m, 4)) 
    dimnames(coefi) <- list(names(object$regcoef),  
                       c("Value", "Std. Error", "t value", "Pr(>|t|)")) 
    stderr <- as.vector(sqrt(diag(object$regcoef.cov))) 
    coefi[, 2] <- stderr 
    coefi[, 3] <- coefi[, 1]/coefi[, 2] 
    coefi[, 4] <- if(df > 0) 2 * (1 - pt(abs(coefi[, 3]), df)) else NA 
    if (correlation) { 
      if (m==1) 
        sinv <- as.matrix(1/stderr) 
      else 
        sinv <- diag(1/(stderr)) 
      regcoef.corr <- sinv %*% object$regcoef.cov %*% sinv 
      dimnames(regcoef.corr) <- dimnames(object$regcoef.cov) 
    } 
  }
  object <- list(call=object$call, ARIMA.model=object$model, 
                 reg.coef=coefi, regcoef.cov=object$regcoef.cov,
                 sigma=object$sigma.innov,
                 outliers=summary(object$outliers),
                 ARIMA.cov=vcov.arima.rob(object, nterm=n))
  ARIMA.stderr <- as.vector(sqrt(diag(object$ARIMA.cov))) 
  if (correlation) {
    object$regcoef.corr <- regcoef.corr
    if (p+q+sma==1) 
      ARIMA.sinv <- as.matrix(1/ARIMA.stderr) 
    else 
      ARIMA.sinv <- diag(1/(ARIMA.stderr)) 
    object$ARIMA.corr <- ARIMA.sinv %*% object$ARIMA.cov %*% ARIMA.sinv 
    dimnames(object$ARIMA.corr) <- dimnames(object$ARIMA.cov) 
  }
  if (p>0) { 
    AR.coef <- array(phi, c(p,4)) 
    dimnames(AR.coef) <- list(names(phi),  
    c("Value", "Std. Error", "t value", "Pr(>|t|)")) 
    AR.coef[, 2] <- ARIMA.stderr[1:p] 
    AR.coef[, 3] <- AR.coef[, 1]/AR.coef[, 2] 
    AR.coef[, 4] <- if(df>0) 2 * (1 - pt(abs(AR.coef[, 3]), df)) else NA 
    object$AR.coef <- AR.coef
  }
  if (q>0) { 
    MA.coef <- array(theta, c(q,4)) 
    dimnames(MA.coef) <- list(names(theta),  
    c("Value", "Std. Error", "t value", "Pr(>|t|)")) 
    MA.coef[, 2] <- ARIMA.stderr[(p+1):(p+q)] 
    MA.coef[, 3] <- MA.coef[, 1]/MA.coef[, 2] 
    MA.coef[, 4] <- if(df>0) 2 * (1 - pt(abs(MA.coef[, 3]), df)) else NA 
    object$MA.coef <- MA.coef
  }
  if (sma) { 
    SMA.coef <- array(seasonal.theta, c(1,4)) 
    dimnames(SMA.coef) <- list(names(SMA.coef),  
    c("Value", "Std. Error", "t value", "Pr(>|t|)")) 
    SMA.coef[, 2] <- ARIMA.stderr[(p+q+1):(p+q+1)] 
    SMA.coef[, 3] <- SMA.coef[, 1]/SMA.coef[, 2] 
    SMA.coef[, 4] <- if(df>0) 2 * (1 - pt(abs(SMA.coef[, 3]), df)) else NA 
    object$sMA.coef <- SMA.coef
  }
  object$n <- n 
  object$df <- df 
  oldClass(object) <- "summary.arima.rob" 
  object         
}
vcov.arima.rob <-
function(object, nterm=1000, ...) 
{ 
  phi = object$model$ar
  theta = object$model$ma
  seasonal.theta = object$model$sma
  seasonal.freq = object$model$sfreq
  p  = length(phi) 
  q  = length(theta) 
  qs = length(seasonal.theta) 
  d1 = d2 = d3 = rep(0, p+q+qs) 
  mat = matrix(NA, p+q+qs, p+q+qs) 
  c.names = NULL 
  if (p > 0) {
    coef1 = -arima.invoper(phi, 1e-6, nterm) 
    n1 = length(coef1) 
    d1[1:p] = 1:p 
    d2[1:p] = 1 
    d3[1:p] = n1
    c.names = c(c.names, names(phi))
  }
  if(q > 0) { 
    coef2 = arima.invoper(theta, 1e-6, nterm) 
    n2 = length(coef2) 
    d1[(p+1):(p+q)] = 1:q 
    d2[(p+1):(p+q)] = 2 
    d3[(p+1):(p+q)] = n2 
    c.names = c(c.names, names(theta))
  }
  if(qs > 0) { 
    coef3p = -arima.invoper(seasonal.theta, 1e-6, nterm) 
    n3p = length(coef3p) 
    n3 = (n3p-1) * seasonal.freq + 1 
    d1[p+q+1] = seasonal.freq 
    d2[p+q+1] = 3 
    d3[p+q+1] = n3 
    coef3 = rep(0, n3) 
    coef3[1:n3] = 0 
    coef3[1] = coef3p[1] 
    coef3[seasonal.freq*(1:(n3p-1))+1] = coef3p[2:n3p] 
    c.names = c(c.names, names(seasonal.theta))
  }
  for (i in 1:(p+q+qs)) { 
    for (j in i:(p+q+qs)) { 
      ii = i 
      jj = j 
      if(d1[ii] > d1[jj]) { 
        ii = j 
        jj = i
      } 
      switch(d2[ii],
             cc1 <- coef1,
             cc1 <- coef2,
             cc1 <- coef3) 
      switch(d2[jj], 
             cc2 <- coef1,
             cc2 <- coef2,
             cc2 <- coef3) 
      nn = min(d3[ii], d3[jj] - d1[jj] + d1[ii]) 
      mat[ii,jj] = sum(cc1[1:nn]*cc2[(1+d1[jj]-d1[ii]):(nn+d1[jj]-d1[ii])]) 
      mat[jj,ii] = mat[ii,jj]
    }
  }
  mat = solve(mat) 
  mat = object$tauef * mat / length(object$y)
  dimnames(mat) = list(c.names, c.names)
  mat
}
