outliers <-
function(object, iter=NULL) 
{
  if (is.null(iter)) 
    x <- object$outliers
  else {
    x <- object$outliers.iter[[iter]]
  }
  x
}
outliers.clean <-
function(object, y)
{
  out <- outliers.series(object, length(y))
  if (!is.null(out$ls.series)) 
    y <- y - out$ls.series
  if (!is.null(out$ao.series))
    y <- y - out$ao.series
  y
}
outliers.combine <-
function(object)
{
  iter <- length(object)
  out2 <- object[[1]]
  out2$sigma <- object[[iter]]$sigma
  hasPos <- !is.null(out2$outlier.positions)
  if (iter > 1) {
    for (i in 2:iter) {
      if (object[[i]]$nout > 0) {
        out2$outlier.type <- c(out2$outlier.type,
                               object[[i]]$outlier.type)
        out2$outlier.index <- c(out2$outlier.index,
                                object[[i]]$outlier.index)
        out2$outlier.impact <- c(out2$outlier.impact,
                                 object[[i]]$outlier.impact)
        out2$outlier.t.statistics <- c(out2$outlier.t.statistics,
                                       object[[i]]$outlier.t.statistics)
        if (hasPos) {
          out2$outlier.positions <- c(out2$outlier.positions,
                                      object[[i]]$outlier.positions)
        }
      }
    }
  }
  sort.index <- sort.list(out2$outlier.index)
  out2$outlier.type <- out2$outlier.type[sort.index]
  out2$outlier.index <- out2$outlier.index[sort.index]
  out2$outlier.impact <- out2$outlier.impact[sort.index]
  out2$outlier.t.statistics <- out2$outlier.t.statistics[sort.index]
  if (hasPos) {
    out2$outlier.positions <- out2$outlier.positions[sort.index]
  }
  out2$nout <- length(sort.index)
  out2
}
outliers.iter.series <-
function(object, n) 
{
  ao.series <- ls.series <- io.series <- numeric(n)
  iter <- length(object)
  for (i in 1:iter) {
    oType <- object[[i]]$outlier.type 
    oImpact <- object[[i]]$outlier.impact
    oIndex <- object[[i]]$outlier.index
    ls.index <- oIndex[oType == 3]
    ao.index <- oIndex[oType == 2]
    io.index <- oIndex[oType == 1]
    n.ls <- length(ls.index)
    n.ao <- length(ao.index)
    n.io <- length(io.index)
    if (n.ao != 0) 
      ao.series[ao.index] <- ao.series[ao.index]+oImpact[oType == 2]
    if (n.io != 0) 
      io.series[io.index] <- io.series[io.index]+oImpact[oType == 1]
    if (n.ls != 0) {
      tmp <- oImpact[oType == 3]
      if (n.ls == 1) 
        ls.series[(ls.index+1):n] <- ls.series[(ls.index+1):n]+tmp
      else {
        for (i in 1:(n.ls-1)) 
          ls.series[(ls.index[i]+1):(ls.index[i+1]-1)] <- 
          ls.series[(ls.index[i]+1):(ls.index[i+1]-1)]+tmp[i]
        ls.series[(ls.index[n.ls]+1):n] <- 
        ls.series[(ls.index[n.ls]+1):n]+tmp[n.ls]
      }
    }
  }
  list(ao.series=ao.series, ls.series=ls.series, io.series=io.series)
}
outliers.rr <-
function(object, innov.outlier=FALSE, critv=NULL)
{
  n = length(object$y)
  m = length(object$regcoef)
  if (is.null(critv)) {
    if(n <= 200) {
      critv = 3
    }
    else if(n <= 500) {
      critv = 3.5
    }
    else {
      critv = 4
    }
  }
  y = object$y.cleaned
  if (is.null(y)) 
    y = object$y
  if (is(y, "timeSeries")) {
    isTS = TRUE
    pos = positions(y)
    y = seriesData(y)
  }
  else {
    isTS = FALSE
  }
  idif = object$model$d
  sfreq = object$model$sfreq
  nsd = object$model$sd
#
  phi = object$model$ar
  p = length(phi)
  theta = object$model$ma
  q = length(theta)
  seasonal.theta = object$model$sma
  indth = length(seasonal.theta)
  if (indth == 0) 
    seasonal.theta = double(1)
  beta = object$regcoef
  if (is.null(beta)) {
    beta = double(1)
  }
  idim = max(q+sfreq+1, p+idif+sfreq*nsd)
  if (idim > p) {
    phi = c(phi, rep(0, idim-p))
  }
  storage.mode(phi) = "double"
  if (idim > q) {
    theta = c(theta, rep(0, idim-q))
  }
  storage.mode(theta) = "double"
  idimw = n*(18+idim+max(m,1)+1) + 7*idim + 7 + 2*(q+sfreq) + 1 + 
          5*idim + 5*idim**2 + 3*(q+sfreq) + 2 + n + 
          max(4*idim+4*idim**2,n)
  idimiw = 2*n+idim+1
  idx = c("nout", "outlier.type", "outlier.index", "outlier.impact",
          "outlier.t.statistics", "sigma0", "sigma", "ierror")
  out = .Fortran(F_s_outlfe,
                 object$x,
                 y,
                 n,
                 m,
                 as.integer(idif),
                 as.integer(sfreq),
                 as.integer(nsd),
                 as.integer(p),
                 as.integer(q),
                 as.integer(indth),
                 beta,
                 phi,
                 theta,
                 as.double(seasonal.theta),
                 object$sigma.regresid,
                 as.integer(innov.outlier),
                 cck=object$tuning.c,
                 object$sigma.first,
                 as.double(critv),
                 nout=integer(1),
                 outlier.type=integer(n),
                 outlier.index=integer(n),
                 outlier.impact=double(n),
                 outlier.t.statistics=double(n),
                 sigma0=object$sigma.innov,
                 sigma=double(1),
                 as.integer(idim),
                 work=double(idimw),
                 as.integer(idimw),
                 iwork=integer(idimiw),
                 as.integer(idimiw),
                 ierror=integer(1),
                 object$n0)[idx]
  if (out$nout > 0) {
    orden <- order(out$outlier.index[1:out$nout])
    out$outlier.type <- out$outlier.type[orden]
    out$outlier.impact <- out$outlier.impact[orden]
    out$outlier.index <- sort(out$outlier.index[1:out$nout])
    out$outlier.t.statistics <- out$outlier.t.statistics[orden]
    if (isTS) {
      out$outlier.positions <- pos[out$outlier.index]
    }
  }
  else {
    out <- list(nout=0)
  }
  oldClass(out) <- "outliers"
  out
}
outliers.series <-
function(object, n)
{
  ls.series <- ao.series <- io.series <- numeric(n)
  oType <- object$outlier.type 
  oImpact <- object$outlier.impact
  oIndex <- object$outlier.index
  ls.index <- oIndex[oType == 3]
  ao.index <- oIndex[oType == 2]
  io.index <- oIndex[oType == 1]
  n.ls <- length(ls.index)
  n.ao <- length(ao.index)
  n.io <- length(io.index)
  if (n.ao != 0) {
    ao.series[ao.index] <- oImpact[oType == 2]
  }
  if (n.io != 0) {
    io.series[io.index] <- oImpact[oType == 1]
  }
  if (n.ls != 0) {
    ls.series[ls.index] <- oImpact[oType == 3]
    if (n.ls == 1) 
      ls.series[(ls.index+1):n] <- ls.series[ls.index]
    else {
      for (i in 1:(n.ls-1)) 
        ls.series[(ls.index[i]+1):(ls.index[i+1]-1)] <- ls.series[ls.index[i]]
        ls.series[(ls.index[n.ls]+1):n] <- ls.series[ls.index[n.ls]]
    }
  }
  list(ao.series=ao.series, ls.series=ls.series, io.series=io.series)
}
print.outliers <-
function(x, digits=4, ...) 
{
  if (x$nout>0) {
    cat("\n Number of outliers detected: ", length(unique(x$outlier.index)))
  }
  else {
    cat("\n Number of outliers detected: 0\n")
    return(invisible(x))
  }
  types <- c("IO", "AO", "LS")
  cat("\n\nOutlier index\n")
  print(x$outlier.index)
  cat("\nOutlier type\n")
  print(types[x$outlier.type])
  cat("\nOutlier impact\n")
  print(format(round(x$outlier.impact, digits=digits)), quote=FALSE, ...)
  cat("\nOutlier t-statistics\n")
  print(format(round(x$outlier.t.statistics, digits=digits)), quote=FALSE, ...)
  invisible(x)
}
print.summary.outliers <-
function(x, digits=4, ...) 
{
  if (x$nout > 0) { 
    if (is.null(x$nout.unique))
      cat("\n Number of outliers detected: ", x$nout, "\n")
    else
      cat("\n Number of outliers detected: ", x$nout.unique, "\n")
  }
  else {
    cat("\n Number of outliers detected: 0\n")
    return(invisible(x))
  }
  cat ("\nOutliers detected:\n\n")
  print(x$outliers.table)
  cat("\nInnovation scale estimate before correcting outliers:\n ", 
      format(x$sigma0, digits=digits), "\n")
  cat("\nInnovation scale estimate after correcting outliers:\n ", 
      format(x$sigma, digits=digits), "\n")
  invisible(x)
}
summary.outliers <-
function(object, ...)
{
  nout <- object$nout
  if (nout > 0) {
    index <- object$outlier.index
    nout.unique <- length(unique(index))
    aux <- c("IO", "AO", "LS")
    type2 <- aux[object$outlier.type]
    outliers.table <- array(index, c(nout, 4))
    oldClass(outliers.table) <- "char.matrix"
    if (!is.null(object$outlier.positions)) {
      dimnames(outliers.table) <- list(1:nout, 
      c("Time", "Type", "Impact", "t-value"))
      outliers.table[, 1] <- as.character(object$outlier.positions)
    }
    else {
      dimnames(outliers.table) <- list(1:nout, 
      c("Index", "Type", "Impact", "t-value")) 
    }
    outliers.table[, 2] <- type2
    outliers.table[, 3] <- format(object$outlier.impact, digits=4)
    outliers.table[, 4] <- format(object$outlier.t.statistics, digits=4)
    object$outliers.table <- outliers.table
    object$outlier.impact <- object$outlier.t.statistics <- 
    object$outlier.index <- object$outlier.type <- NULL
    object$outlier.positions <- NULL
    if (nout.unique != nout) 
      object$nout.unique <- nout.unique
  }
  else {
    object <- NULL
    object$nout <- nout
  }
  oldClass(object) <- "summary.outliers"
  object
}
