npregression <- 
    function(bws, eval, mean, merr = NA, grad = NA, gerr = NA,
             resid = NA,
             ntrain, trainiseval = FALSE, gradients = FALSE, residuals = FALSE,
             gradient.order = NULL,
             xtra = rep(NA, 6),
             rows.omit = NA,
             train.rows.omit = NULL,
             eval.rows.omit = NULL,
             timing = NA,
             total.time = NA,
             optim.time = NA,
             fit.time = NA){

        if (missing(bws) || missing(eval) || missing(ntrain))
            stop("improper invocation of npregression constructor")

        if (length(rows.omit) == 0)
            rows.omit <- NA
        if (is.null(train.rows.omit) || length(train.rows.omit) == 0)
            train.rows.omit <- NA
        if (is.null(eval.rows.omit) || length(eval.rows.omit) == 0)
            eval.rows.omit <- NA

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
            nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
            train.rows.omit = train.rows.omit,
            train.nobs.omit = if (identical(train.rows.omit, NA)) 0 else length(train.rows.omit),
            eval.rows.omit = eval.rows.omit,
            eval.nobs.omit = if (identical(eval.rows.omit, NA)) 0 else length(eval.rows.omit),
            timing = timing,
            total.time = total.time,
            optim.time = optim.time,
            fit.time = fit.time)

        class(d) <- "npregression"

        return(d)
    }

print.npregression <- function(x, digits=NULL, ...){
  cat("\nRegression Data: ", x$ntrain, " training points,",
      if (x$trainiseval) "" else paste(" and ", x$nobs," evaluation points,",
                                      sep=""),
      " in ",x$ndim," variable(s)\n",sep="")

  printSearchParameterSummary(x$bw, x$xnames, x$bws, vari = "x",
                              fallback.label = paste(x$pscaling, ":", sep = ""),
                              digits = digits)
  
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
  errors <- npValidateScalarLogical(errors, "errors")
  gout <- if (!errors) x$grad else x$gerr
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop(if (!errors)
      "gradients are not available: fit the model with gradients=TRUE"
    else
      "gradient standard errors are not available: fit the model with gradients=TRUE")

  if (identical(x$bws$regtype, "lc") && !is.null(gradient.order)) {
    npValidateLcGradientOrder(
      regtype = x$bws$regtype,
      gradient.order = gradient.order,
      ncon = x$bws$ncon,
      where = "gradients.npregression"
    )
  }
  if (isTRUE(x$bws$ncon == 0L)) {
    if (!is.null(gradient.order)) {
      npValidateCategoricalFirstDifferenceGradientOrder(
        gradient.order = gradient.order,
        where = "gradients.npregression"
      )
    }
    return(gout)
  }
  if (!identical(x$bws$regtype, "lp") || is.null(gradient.order))
    return(gout)

  if (!is.matrix(gout))
    return(gout)

  gorder <- npValidateGlpGradientOrder(regtype = x$bws$regtype,
                                       gradient.order = gradient.order,
                                       ncon = x$bws$ncon)
  stored.order <- x$gradient.order
  if (is.null(stored.order) && x$bws$ncon > 0L)
    stored.order <- rep.int(1L, x$bws$ncon)
  stored.order <- npValidateGlpGradientOrder(regtype = x$bws$regtype,
                                             gradient.order = stored.order,
                                             ncon = x$bws$ncon)
  available <- npGlpGradientAvailability(
    regtype.engine = x$bws$regtype,
    degree.engine = x$bws$degree,
    gradient.order = gorder,
    ncon = x$bws$ncon,
    where = "gradients.npregression"
  )
  if (!any(available) && (x$bws$nuno + x$bws$nord == 0L)) {
    npStopGlpGradientNoneAvailable(
      where = "gradients.npregression",
      action = "return",
      degree.engine = x$bws$degree,
      gradient.order = gorder,
      available = available,
      con.names = x$bws$xnames[x$bws$icon]
    )
  }
  if (any(!available)) {
    npWarnGlpGradientPartialAvailability(
      where = "gradients.npregression",
      degree.engine = x$bws$degree,
      gradient.order = gorder,
      available = available,
      con.names = x$bws$xnames[x$bws$icon]
    )
  }
  if (length(gorder)) {
    defined.request <- available
    if (any(defined.request & (gorder != stored.order))) {
      stop("requested gradient.order differs from the derivative order stored in this npregression object; refit or predict/evaluate with gradients=TRUE and the desired gradient.order",
           call. = FALSE)
    }
  }

  gout.masked <- gout
  cont.idx <- which(x$bws$icon)
  if (length(cont.idx)) {
    drop.idx <- cont.idx[!available]
    if (length(drop.idx))
      gout.masked[, drop.idx] <- NA_real_
  }
  gout.masked
}
predict.npregression <- function(object, se.fit = FALSE, ...) {
  se.fit <- npValidateScalarLogical(se.fit, "se.fit")
  dots <- list(...)
  has.formula.route <- !is.null(object$bws$formula)

  if (!is.null(dots$exdat) && !is.null(dots$newdata)) {
    dots$newdata <- NULL
  } else if (!has.formula.route && is.null(dots$exdat) && !is.null(dots$newdata)) {
    dots$exdat <- dots$newdata
    dots$newdata <- NULL
  }

  tr <- do.call(npreg, c(list(bws = object$bws), dots))
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs, residual.scale = tr$MSE))
  else
    return(fitted(tr))
}

summary.npregression <- function(object, ...) {
  cat("\nRegression Data: ", object$ntrain, " training points,",
      if (object$trainiseval) "" else paste(" and ", object$nobs," evaluation points,",
                                      sep=""),
      " in ",object$ndim," variable(s)\n",sep="")

  cat(genOmitStr(object))

  printSearchParameterSummary(object$bw, object$xnames, object$bws, vari = "x",
                              fallback.label = paste(object$pscaling, ":", sep = ""))

  cat(genRegEstStr(object))
  cat(genGofStr(object))

  cat(genBwKerStrs(object$bws))
  cat(genTimingStr(object))
  cat('\n\n')  
}
