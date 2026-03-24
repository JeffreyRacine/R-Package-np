condistribution <- 
    function(bws, xeval, yeval, condist, conderr = NA,
             congrad = NA, congerr = NA,           
             ntrain, trainiseval = FALSE, gradients = FALSE,
             proper.requested = FALSE,
             proper.applied = FALSE,
             proper.method = NULL,
             condist.raw = NULL,
             proper.info = NULL,
             rows.omit = NA,
             timing = NA, total.time = NA,
             optim.time = NA, fit.time = NA){

        if (missing(bws) || missing(xeval) || missing(yeval) || missing(condist) || missing(ntrain))
            stop("improper invocation of condistribution constructor")

        if (length(rows.omit) == 0)
            rows.omit <- NA

        d <- list(
            xbw = bws$xbw,
            ybw = bws$ybw,
            bws = bws,
            xnames = bws$xnames,
            ynames = bws$ynames,
            nobs = dim(xeval)[1],
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
            yeval = yeval,
            condist = condist,
            conderr = conderr,
            congrad = congrad,
            congerr = congerr,
            ntrain = ntrain,
            trainiseval = trainiseval,
            gradients = gradients,
            proper.requested = proper.requested,
            proper.applied = proper.applied,
            proper.method = proper.method,
            condist.raw = condist.raw,
            proper.info = proper.info,
            rows.omit = rows.omit,
            nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
            timing = timing, total.time = total.time,
            optim.time = optim.time, fit.time = fit.time)


        class(d) <- "condistribution"

        return(d)
    }

print.condistribution <- function(x, digits=NULL, ...){
  cat("\nConditional distribution data: ", x$ntrain, " training points,",
      if (x$trainiseval) "" else paste(" and ", x$nobs, " evaluation points,\n", sep=""),
      " in ", x$xndim + x$yndim, " variable(s)",
      "\n(", x$yndim, " dependent variable(s), and ", x$xndim, " explanatory variable(s))\n\n",
      sep="")
  print(matrix(x$ybw,ncol=x$yndim,dimnames=list(paste("Dep. Var. ",x$pscaling,":",sep=""),x$ynames)))

  print(matrix(x$xbw,ncol=x$xndim,dimnames=list(paste("Exp. Var. ",x$pscaling,":",sep=""),x$xnames)))

  cat(genDenEstStr(x))
  cat(genBwKerStrs(x$bws))
  if (!is.null(x$proper.requested) && !is.null(x$proper.applied)) {
    proper.state <- if (isTRUE(x$proper.applied)) {
      sprintf("requested and applied (%s)", x$proper.method)
    } else if (isTRUE(x$proper.requested)) {
      "requested but not applied"
    } else {
      "not requested"
    }
    cat("\nProper distribution repair:", proper.state)
  }


  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

fitted.condistribution <- function(object, ...){
 object$condist 
}
se.condistribution <- function(x){
  if (isTRUE(x$proper.applied)) {
    stop("standard errors are unavailable for repaired conditional distributions in tranche 1")
  }
  x$conderr
}
gradients.condistribution <- function(x, errors = FALSE, ...) {
  gout <- if (!errors) x$congrad else x$congerr
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop(if (!errors)
      "gradients are not available: fit the model with gradients=TRUE"
    else
      "gradient standard errors are not available: fit the model with gradients=TRUE")
  gout
}

predict.condistribution <- function(object, se.fit = FALSE, ...) {
  dots <- list(...)
  has.formula.route <- !is.null(object$bws$formula)
  proper_arg <- dots[["proper", exact = TRUE]]

  if (!has.formula.route &&
      is.null(dots$exdat) &&
      is.null(dots$eydat) &&
      !is.null(dots$newdata)) {
    nd <- toFrame(dots$newdata)
    req <- c(object$ynames, object$xnames)
    miss <- setdiff(req, names(nd))
    if (length(miss) > 0L) {
      stop(sprintf("'newdata' must include columns %s, or supply both 'exdat' and 'eydat'.",
                   paste(shQuote(req), collapse = ", ")))
    }
    dots$eydat <- nd[, object$ynames, drop = FALSE]
    dots$exdat <- nd[, object$xnames, drop = FALSE]
    dots$newdata <- NULL
  }

  if (is.null(proper_arg) && isTRUE(object$proper.requested)) {
    dots$proper <- TRUE
    proper_arg <- TRUE
  }
  if (isTRUE(proper_arg)) {
    proper.control <- dots[["proper.control", exact = TRUE]]
    if (is.null(proper.control))
      proper.control <- list()
    if (is.null(proper.control$fail.on.unsupported))
      proper.control$fail.on.unsupported <- TRUE
    dots$proper.control <- proper.control
  }

  tr <- do.call(npcdist, c(list(bws = object$bws), dots))
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs))
  else
    return(fitted(tr))
}


summary.condistribution <- function(object, ...){
  cat("\nConditional Distribution Data: ", object$ntrain, " training points,",
      if (object$trainiseval) "" else paste(" and ", object$nobs, " evaluation points,\n", sep=""),
      " in ", object$xndim + object$yndim, " variable(s)",
      "\n(", object$yndim, " dependent variable(s), and ", object$xndim, " explanatory variable(s))\n\n",
      sep="")

  cat(genOmitStr(object))
  print(matrix(object$ybw,ncol=object$yndim,dimnames=list(paste("Dep. Var. ",object$pscaling,":",sep=""),object$ynames)))

  print(matrix(object$xbw,ncol=object$xndim,dimnames=list(paste("Exp. Var. ",object$pscaling,":",sep=""),object$xnames)))

  cat(genDenEstStr(object))

  cat(genBwKerStrs(object$bws))
  if (!is.null(object$proper.requested) && !is.null(object$proper.applied)) {
    proper.state <- if (isTRUE(object$proper.applied)) {
      sprintf("requested and applied (%s)", object$proper.method)
    } else if (isTRUE(object$proper.requested)) {
      "requested but not applied"
    } else {
      "not requested"
    }
    cat("\nProper distribution repair:", proper.state)
  }
  cat(genTimingStr(object))
  cat('\n\n')  
}
