npregression <- 
    function(bws, eval, mean, merr = NA, grad = NA, gerr = NA,
             resid = NA,
             ntrain, trainiseval = FALSE, gradients = FALSE, residuals = FALSE,
             gradient.order = NULL,
             xtra = rep(NA, 6),
             rows.omit = NA){

        if (missing(bws) | missing(eval) | missing(ntrain))
            stop("improper invocation of npregression constructor")

        if (length(rows.omit) == 0)
            rows.omit <- NA

        d <- list(
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
            gradient.order = gradient.order,
            residuals = residuals,
            R2 = xtra[1],
            MSE = xtra[2],
            MAE = xtra[3],
            MAPE = xtra[4],
            CORR = xtra[5],
            SIGN = xtra[6],
            rows.omit = rows.omit,
            nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))

        class(d) <- "npregression"

        return(d)
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
gradients.npregression <- function(x, errors = FALSE, gradient.order = NULL, ...) {
  gout <- if (!errors) x$grad else x$gerr
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop(if (!errors)
      "gradients are not available: fit the model with gradients=TRUE"
    else
      "gradient standard errors are not available: fit the model with gradients=TRUE")

  if (!identical(x$bws$regtype, "glp") || is.null(gradient.order))
    return(gout)

  if (!is.matrix(gout))
    return(gout)

  gorder <- npValidateGlpGradientOrder(regtype = x$bws$regtype,
                                       gradient.order = gradient.order,
                                       ncon = x$bws$ncon)
  gout.masked <- gout
  gout.masked[,] <- NA_real_
  cont.idx <- which(x$bws$icon)
  if (length(cont.idx)) {
    keep.cont <- (gorder == 1L) & (x$bws$glp.degree >= 1L)
    if (any(keep.cont)) {
      keep.idx <- cont.idx[keep.cont]
      gout.masked[, keep.idx] <- gout[, keep.idx, drop = FALSE]
    }
    if (any(gorder > x$bws$glp.degree))
      warning("some requested glp derivatives exceed polynomial degree; returning NA for those components")
    if (any(gorder > 1L))
      warning("higher-order glp derivatives are not yet available at C level; returning NA for requested orders > 1")
  }
  gout.masked
}
predict.npregression <- function(object, se.fit = FALSE, ...) {
  tr <- eval(npreg(bws = object$bws, ...), envir = parent.frame())
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs, residual.scale = tr$MSE))
  else
    return(fitted(tr))
}

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
