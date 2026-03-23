singleindex = 
  function(bws,  betavcov = NULL, index, mean, merr = NA,
           grad = NA, gerr = NA,
           mean.grad = NA, mean.gerr = NA,
           resid = NA,
           ntrain, trainiseval = FALSE, residuals = FALSE,
           gradients = FALSE, xtra = NA,
           confusion.matrix = NA, CCR.overall = NA,
           CCR.byoutcome =  NA, fit.mcfadden = NA
           , timing = NA, total.time = NA,
           optim.time = NA, fit.time = NA
           ){

    if (missing(bws) || missing(index) || missing(mean) || missing(ntrain))
      stop("improper invocation of singleindex constructor")

    d = list(
      beta = bws$beta,
      betavcov = betavcov,
      bw = bws$bw,
      bws = bws,
      pregtype = bws$pregtype,
      pmethod = bws$pmethod,
      xnames = bws$xnames,
      ynames = bws$ynames,
      nobs = dim(index)[1],
      ndim = bws$ndim,
      nord = bws$nord,
      nuno = bws$nuno,
      ncon = bws$ncon,
      pckertype = bws$pckertype,
      index = index,
      mean = mean,
      merr = merr,
      grad = grad,
      gerr = gerr,
      mean.grad = mean.grad,
      mean.gerr = mean.gerr,
      resid = resid,
      ntrain = ntrain,
      trainiseval = trainiseval,
      residuals = residuals,
      gradients = gradients,
      R2 = xtra[1],
      MSE = xtra[2],
      MAE = xtra[3],
      MAPE = xtra[4],
      CORR = xtra[5],
      SIGN = xtra[6],
      confusion.matrix = confusion.matrix,
      CCR.byoutcome = CCR.byoutcome,
      CCR.overall = CCR.overall,
      fit.mcfadden = fit.mcfadden,
      timing = timing,
      total.time = total.time,
      optim.time = optim.time,
      fit.time = fit.time)

    class(d) = "singleindex"

    d
  }

print.singleindex <- function(x, digits=NULL, ...){
  cat("\nSingle Index Model",
      "\nRegression data: ", x$ntrain, " training points,",
      if (x$trainiseval) "" else paste(" and ", x$nobs,
                                      " evaluation points,", sep=""),
      " in ",x$ndim," variable(s)\n",sep="")

  print(matrix(x$beta,ncol=x$ndim,dimnames=list(paste("Beta",":",sep=""),x$xnames)))
  
  cat("Bandwidth:", x$bw)
  cat(genRegEstStr(x))
  cat(genBwKerStrs(x$bws))
  
  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

coef.singleindex <- function(object, ...) {
  tc <- object$beta
  names(tc) <- object$xnames
  return(tc)
}
vcov.singleindex <- function(object, ...) {
  tc <- object$betavcov
  if(!is.null(tc)) {
    return(tc)
  } else {
    .np_warning("variance-covariance matrix does not exist: verify gradients=TRUE")
  }
}
fitted.singleindex <- function(object, ...){
 object$mean 
}
residuals.singleindex <- function(object, ...) {
 if(object$residuals) { return(object$resid) } else { return(npindex(bws = object$bws, residuals =TRUE)$resid) } 
}
predict.singleindex <- function(object, se.fit = FALSE, ...) {
  dots <- list(...)
  has.formula.route <- !is.null(object$bws$formula)

  if (!has.formula.route && is.null(dots$exdat) && !is.null(dots$newdata)) {
    dots$exdat <- dots$newdata
    dots$newdata <- NULL
  }

  ## When no new evaluation inputs are supplied, reuse stored fit
  ## directly. This avoids unnecessary recomputation and sidesteps
  ## autodispatch temporary call symbols that may no longer exist.
  if (length(dots) == 0L) {
    if (se.fit) {
      return(list(
        fit = fitted(object),
        se.fit = se(object),
        df = object$nobs,
        residual.scale = object$MSE
      ))
    }
    return(fitted(object))
  }

  tr <- do.call(npindex, c(list(bws = object$bws, errors = se.fit), dots))
  if (se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr),
                df = tr$nobs, residual.scale = tr$MSE))
  fitted(tr)
}
se.singleindex <- function(x){ x$merr }
gradients.singleindex <- function(x, errors = FALSE, ...) {
  gout <- if (!errors) x$grad else x$gerr
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop(if (!errors)
      "gradients are not available: fit the model with gradients=TRUE"
    else
      "gradient standard errors are not available: fit the model with gradients=TRUE")
  gout
}

summary.singleindex <- function(object, ...){
  cat("\nSingle Index Model",
      "\nRegression Data: ", object$ntrain, " training points,",
      if (object$trainiseval) "" else paste(" and ", object$nobs,
                                      " evaluation points,", sep=""),
      " in ",object$ndim," variable(s)\n\n",sep="")

  cat(genOmitStr(object))
  print(matrix(object$beta,ncol=object$ndim,dimnames=list(paste("Beta",":",sep=""),object$xnames)))
  
  cat("Bandwidth:", object$bw)
  cat(genRegEstStr(object))
  cat("\n")
  cat(genGofStr(object))
  pCatGofStr(object)
  cat(genBwKerStrs(object$bws))
  cat(genTimingStr(object))
  
  cat("\n\n")
}
