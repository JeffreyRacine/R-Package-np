qregression <- 
    function(bws, xeval, tau, quantile, quanterr = NA, quantgrad = NA, ntrain, trainiseval = FALSE, gradients = FALSE,
             timing = NA, total.time = NA, optim.time = NA, fit.time = NA){

        if (missing(bws) | missing(xeval) | missing(tau) | missing(quantile) | missing(ntrain))
            stop("improper invocation of qregression constructor")

        d <- list(
            xbw = bws$xbw,
            ybw = bws$ybw,
            bws = bws,
            xnames = bws$xnames,
            ynames = bws$ynames,
            nobs = nrow(xeval),
            xndim = bws$xndim,
            yndim = bws$yndim,
            xnord = bws$xnord,
            xnuno = bws$xnuno,
            xncon = bws$xncon,
            ynord = bws$ynord,
            ynuno = bws$ynuno,
            yncon = bws$yncon,
            pscaling = bws$pscaling,
            ptype = bws$ptype,
            pcxkertype = bws$pcxkertype,
            puxkertype = bws$puxkertype,
            poxkertype = bws$poxkertype,
            pcykertype = bws$pcykertype,
            puykertype = bws$puykertype,
            poykertype = bws$poykertype,
            xeval = xeval,
            tau = tau,
            quantile = quantile,
            quanterr = quanterr,
            quantgrad = quantgrad,
            ntrain = ntrain,
            trainiseval = trainiseval,
            gradients = gradients,
            timing = timing, total.time = total.time,
            optim.time = optim.time, fit.time = fit.time)

        class(d) <- "qregression"

        return(d)
    }

print.qregression <- function(x, digits=NULL, ...){
  cat("\nQuantile regression data: ", x$ntrain, " training points,",
      ifelse(x$trainiseval, "", paste(" and ", x$nobs, " evaluation points,\n", sep="")),
      " in ", x$xndim + x$yndim, " variable(s)",
      "\n(", x$yndim, " dependent variable(s), and ", x$xndim, " explanatory variable(s))\n\n",
      sep="")
  print(matrix(x$ybw,ncol=x$yndim,dimnames=list(paste("Dep. Var. ",x$pscaling,":",sep=""),x$ynames)))

  print(matrix(x$xbw,ncol=x$xndim,dimnames=list(paste("Exp. Var. ",x$pscaling,":",sep=""),x$xnames)))

  cat(genRegEstStr(x))
  cat(genBwKerStrs(x$bws))
  
  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

fitted.qregression <- function(object, ...){
 object$quantile 
}
quantile.qregression <- function(x, ...){ x$quantile }
predict.qregression <- function(object, se.fit = FALSE, ...) {
  tr <- eval(npqreg(bws = object$bws, ...), envir = parent.frame())
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs, residual.scale = NA))
  else
    return(fitted(tr))
}

se.qregression <- function(x) { x$quanterr }
gradients.qregression <- function(x, errors = FALSE, ...) {
  if (errors)
    stop("gradient standard errors are not available for qregression")
  gout <- x$quantgrad
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop("gradients are not available: fit the model with gradients=TRUE")
  gout
}

summary.qregression <- function(object, ...) {
  cat("\nQuantile Regression Data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "", paste(" and ", object$nobs, " evaluation points,\n", sep="")),
      " in ", object$xndim + object$yndim, " variable(s)",
      "\n(", object$yndim, " dependent variable(s), and ", object$xndim, " explanatory variable(s))\n\n",
      sep="")

  cat(genOmitStr(object))

  print(matrix(object$ybw,ncol=object$yndim,dimnames=list(paste("Dep. Var. ",object$pscaling,":",sep=""),object$ynames)))

  print(matrix(object$xbw,ncol=object$xndim,dimnames=list(paste("Exp. Var. ",object$pscaling,":",sep=""),object$xnames)))

  cat(genRegEstStr(object))
  cat(genBwKerStrs(object$bws))
  cat(genTimingStr(object))
  
  cat("\n\n")
}
