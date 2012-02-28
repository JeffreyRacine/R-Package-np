npregression = 
  function(bws, eval, mean, merr = NA, grad = NA, gerr = NA,
           resid = NA,
           ntrain, trainiseval = FALSE, gradients = FALSE, residuals = FALSE,
           xtra = rep(NA, 6),
           rows.omit = NA){

    if (missing(bws) | missing(eval) | missing(ntrain))
      stop("improper invocation of npregression constructor")

    if (length(rows.omit) == 0)
      rows.omit <- NA

    d = list(
      bw = bws$bw,
      bws = bws,
      pregtype = bws$pregtype,
      xnames = bws$xnames,
      ynames = bws$ynames,
      nobs = dim(eval)[1],
      ndim = bws$ndim,
      nord = bws$nord,
      nuno = bws$nuno,
      ncon = bws$ncon,
      pscaling = bws$pscaling,
      ptype = bws$ptype,
      pckertype = bws$pckertype,
      pukertype = bws$pukertype,
      pokertype = bws$pokertype,
      eval = eval,
      mean = mean,
      merr = merr,
      grad = grad,
      gerr = gerr,
      resid = resid,
      ntrain = ntrain,
      trainiseval = trainiseval,
      gradients = gradients,
      residuals = residuals,
      R2 = xtra[1],
      MSE = xtra[2],
      MAE = xtra[3],
      MAPE = xtra[4],
      CORR = xtra[5],
      SIGN = xtra[6],
      rows.omit = rows.omit,
      nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))

    class(d) = "npregression"

    d
  }

print.npregression <- function(x, digits=NULL, ...){
  cat("\nRegression Data: ", x$ntrain, " training points,",
      ifelse(x$trainiseval, "", paste(" and ", x$nobs," evaluation points,",
                                      sep="")),
      " in ",x$ndim," variable(s)\n",sep="")

  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$xnames)))
  
  cat(genRegEstStr(x))
  cat(genBwKerStrs(x$bws))
  
  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

fitted.npregression <- function(object, ...){
 object$mean 
}
residuals.npregression <- function(object, ...) {
 if(object$residuals) { return(object$resid) } else { return(npreg(bws = object$bws, residuals =TRUE)$resid) } 
}
se.npregression <- function(x) { x$merr }
gradients.npregression <- function(x, errors = FALSE, ...) {
  if(!errors)
    return(x$grad)
  else
    return(x$gerr)
}
predict.npregression <- function(object, se.fit = FALSE, ...) {
  tr <- eval(npreg(bws = object$bws, ...), envir = parent.frame())
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs, residual.scale = tr$MSE))
  else
    return(fitted(tr))
}
plot.npregression <- function(x, ...) { npplot(bws = x$bws, ...) }

summary.npregression <- function(object, ...) {
  cat("\nRegression Data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "", paste(" and ", object$nobs," evaluation points,",
                                      sep="")),
      " in ",object$ndim," variable(s)\n",sep="")

  cat(genOmitStr(object))

  print(matrix(object$bw,ncol=object$ndim,dimnames=list(paste(object$pscaling,":",sep=""),object$xnames)))

  cat(genRegEstStr(object))
  cat(genGofStr(object))

  cat(genBwKerStrs(object$bws))
  cat('\n\n')  
}
