npdensity <- 
    function(bws, eval, dens,
             derr = NA, ll = NA,
             ntrain, trainiseval = FALSE,
             rows.omit = NA){

        if (missing(bws) | missing(eval) | missing(dens) | missing(ntrain))
            stop("improper invocation of density constructor")

        if (length(rows.omit) == 0)
            rows.omit <- NA


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
            dens = dens,
            derr = derr,
            log_likelihood = ll,
            ntrain = ntrain,
            trainiseval = trainiseval,
            rows.omit = rows.omit,
            nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))


        class(d) <- "npdensity"

        return(d)
    }

print.npdensity <- function(x, digits=NULL, ...){
  cat("\nDensity Data: ", x$ntrain, " training points,",
      ifelse(x$trainiseval, "", paste(" and ", x$nobs,
                                      " evaluation points,", sep="")),
      " in ",x$ndim," variable(s)\n",sep="")

  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$xnames)))

  cat(genDenEstStr(x))
  
  cat(genBwKerStrs(x$bws))
  
  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

fitted.npdensity <- function(object, ...){
 object$dens 
}
se.npdensity <- function(x){ x$derr }
plot.npdensity <- function(x, ...) {
  dots <- list(...)

  ## For the simple 1D default plot, draw the stored eval/density pairs
  ## directly so overlays like points(x, fitted(obj)) align exactly.
  direct.args <- c("plot.behavior", "plot.errors.method", "plot.errors.type",
                   "plot.errors.boot.num", "plot.errors.boot.method",
                   "plot.errors.alpha", "perspective", "gradients",
                   "xdat", "data", "neval", "xtrim", "xq")
  use.direct <- isTRUE(x$ndim == 1) &&
    isTRUE(x$trainiseval) &&
    !any(direct.args %in% names(dots))

  if (use.direct) {
    ex <- x$eval[[1]]
    if (!is.factor(ex)) {
      ord <- order(ex)
      xlab.val <- if (!is.null(dots$xlab)) dots$xlab else gen.label(x$xnames[1], "X1")
      ylab.val <- if (!is.null(dots$ylab)) dots$ylab else "Density"
      type.val <- if (!is.null(dots$type)) dots$type else "l"

      dots$xlab <- NULL
      dots$ylab <- NULL
      dots$type <- NULL

      do.call(plot, c(list(x = as.numeric(ex[ord]),
                           y = as.numeric(x$dens[ord]),
                           type = type.val,
                           xlab = xlab.val,
                           ylab = ylab.val),
                      dots))
      return(invisible(x))
    }
  }

  npplot(bws = x$bws, ...)
}

predict.npdensity <- function(object, se.fit = FALSE, ...) {
  tr <- eval(npudens(bws = object$bws, ...), envir = parent.frame())
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs, log.likelihood = tr$ll))
  else
    return(fitted(tr))
}


summary.npdensity <- function(object, ...) {
  cat("\nDensity Data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "", paste(" and ", object$nobs,
                                      " evaluation points,", sep="")),
      " in ",object$ndim," variable(s)\n",sep="")

  cat(genOmitStr(object))

  print(matrix(object$bw,ncol=object$ndim,dimnames=list(paste(object$pscaling,":",sep=""),object$xnames)))

  cat(genDenEstStr(object))

  cat(genBwKerStrs(object$bws))
  cat('\n\n')  

}
