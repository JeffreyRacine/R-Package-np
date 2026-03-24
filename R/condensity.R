condensity <- 
    function(bws, xeval, yeval,
             condens, conderr = NA,
             congrad = NA, congerr = NA,
             ll = NA, ntrain, trainiseval = FALSE, gradients = FALSE,
             proper.requested = FALSE,
             proper.applied = FALSE,
             proper.method = NULL,
             condens.raw = NULL,
             proper.info = NULL,
             rows.omit = NA,
             timing = NA, total.time = NA,
             optim.time = NA, fit.time = NA){

        if (missing(bws) || missing(xeval) || missing(yeval) || missing(condens) || missing(ntrain))
            stop("improper invocation of condensity constructor")

        if (length(rows.omit) == 0)
            rows.omit <- NA

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
            yeval = yeval,
            condens = condens,
            conderr = conderr,
            congrad = congrad,
            congerr = congerr,
            log_likelihood = ll,
            ntrain = ntrain,
            trainiseval = trainiseval,
            gradients = gradients,
            proper.requested = proper.requested,
            proper.applied = proper.applied,
            proper.method = proper.method,
            condens.raw = condens.raw,
            proper.info = proper.info,
            rows.omit = rows.omit,
            nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
            timing = timing, total.time = total.time,
            optim.time = optim.time, fit.time = fit.time)
        
        class(d) <- "condensity"

        return(d)
    }

print.condensity <- function(x, digits=NULL, ...){
  cat("\nConditional Density Data: ", x$ntrain, " training points,",
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
    cat("\nProper density repair:", proper.state)
  }

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

fitted.condensity <- function(object, ...){
 object$condens 
}
se.condensity <- function(x){
  if (isTRUE(x$proper.applied)) {
    stop("standard errors are unavailable for repaired conditional densities in tranche 1")
  }
  x$conderr
}
gradients.condensity <- function(x, errors = FALSE, ...) {
  gout <- if (!errors) x$congrad else x$congerr
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop(if (!errors)
      "gradients are not available: fit the model with gradients=TRUE"
    else
      "gradient standard errors are not available: fit the model with gradients=TRUE")
  gout
}

predict.condensity <- function(object, se.fit = FALSE, ...) {
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

  tr <- do.call(npcdens, c(list(bws = object$bws), dots))
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs, log.likelihood = tr$log_likelihood))
  else
    return(fitted(tr))
}



summary.condensity <- function(object, ...){
  cat("\nConditional Density Data: ", object$ntrain, " training points,",
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
    cat("\nProper density repair:", proper.state)
  }
  cat(genTimingStr(object))
  cat('\n\n')  
}
