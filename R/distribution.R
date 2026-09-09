npdistribution <- 
    function(bws, eval, dist, derr = NA,
             ntrain, trainiseval = FALSE,
             rows.omit = NA,
             train.rows.omit = rows.omit,
             eval.rows.omit = if (trainiseval) rows.omit else NA,
             timing = NA, total.time = NA,
             optim.time = NA, fit.time = NA, se = TRUE){

        if (missing(bws) || missing(eval) || missing(dist) || missing(ntrain))
            stop("improper invocation of distribution constructor")

        if (length(rows.omit) == 0)
            rows.omit <- NA
        if (length(train.rows.omit) == 0)
            train.rows.omit <- NA
        if (length(eval.rows.omit) == 0)
            eval.rows.omit <- NA

        d <- list(
            bw = bws$bw,
            bws = bws,
            xnames = bws$xnames,
            nobs = nrow(eval),
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
            dist = dist,
            derr = derr,
            se = se,
            ntrain = ntrain,
            trainiseval = trainiseval,
            rows.omit = rows.omit,
            nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
            train.rows.omit = train.rows.omit,
            ntrain.omit = if (identical(train.rows.omit, NA)) 0 else length(train.rows.omit),
            eval.rows.omit = eval.rows.omit,
            neval.omit = if (identical(eval.rows.omit, NA)) 0 else length(eval.rows.omit),
            timing = timing, total.time = total.time,
            optim.time = optim.time, fit.time = fit.time)

        class(d) <- "npdistribution"

        return(d)
    }

print.npdistribution <- function(x, digits=NULL, ...){
  cat("\nDistribution Data: ", x$ntrain, " training points,",
      if (x$trainiseval) "" else paste(" and ", x$nobs,
                                      " evaluation points,", sep=""),
      " in ",x$ndim," variable(s)\n",sep="")

  printSearchParameterSummary(x$bw, x$xnames, x$bws, vari = "x",
                              fallback.label = paste(x$pscaling, ":", sep = ""),
                              digits = digits)

  
  cat(genDenEstStr(x))
  
  cat(genBwKerStrs(x$bws))
  
  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

fitted.npdistribution <- function(object, ...){
 object$dist 
}
se.npdistribution <- function(x) {
  .np_require_stored_se(
    x, x[["derr", exact = TRUE]], "npudist", expr = substitute(x),
    data.hint = if (is.null(x[["bws", exact = TRUE]][["call", exact = TRUE]]))
      "tdat = training_data" else NULL)
}

predict.npdistribution <- function(object, se.fit = FALSE, ...) {
  se.fit <- npValidateScalarLogical(se.fit, "se.fit")
  dots <- list(...)
  if ("se" %in% names(dots)) {
    requested.se <- npValidateScalarLogical(dots[["se"]], "se")
    if (!identical(requested.se, se.fit))
      stop("conflicting 'se' and 'se.fit' requests; use se.fit to request prediction standard errors",
           call. = FALSE)
    dots[["se"]] <- NULL
  }
  has.formula.route <- !is.null(object$bws$formula)

  if (!is.null(dots$edat) && !is.null(dots$newdata)) {
    dots$newdata <- NULL
  } else if (!has.formula.route && is.null(dots$edat) && !is.null(dots$newdata)) {
    dots$edat <- dots$newdata
    dots$newdata <- NULL
  }

  tr <- do.call(npudist, c(list(bws = object$bws, se = se.fit), dots))
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs))
  else
    return(fitted(tr))
}

summary.npdistribution <- function(object, ...) {
  cat("\nDistribution Data: ", object$ntrain, " training points,",
      if (object$trainiseval) "" else paste(" and ", object$nobs,
                                      " evaluation points,", sep=""),
      " in ",object$ndim," variable(s)\n",sep="")
  
  cat(genOmitStr(object))

  printSearchParameterSummary(object$bw, object$xnames, object$bws, vari = "x",
                              fallback.label = paste(object$pscaling, ":", sep = ""))

  cat(genDenEstStr(object))

  cat(genBwKerStrs(object$bws))
  cat(genTimingStr(object))
  cat('\n\n')

}
