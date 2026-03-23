.np_plreg_bws <- function(object, where = "plregression") {
  if (!is.null(object$bws))
    return(object$bws)

  if (!is.null(object$bw) && inherits(object$bw, "plbandwidth"))
    return(object$bw)

  stop(sprintf("%s object does not contain a usable partially linear bandwidth object", where),
       call. = FALSE)
}

plregression = 
  function(bws, xcoef, xcoeferr = 0, xcoefvcov, evalx, evalz, mean, resid = NA,
           ntrain, trainiseval = FALSE, residuals = FALSE,
           xtra = double(6),
           timing = NA, total.time = NA,
           optim.time = NA, fit.time = NA){

    if (missing(bws) || missing(evalx) || missing(evalz) || missing(mean) || missing(ntrain) || missing(xcoef))
      stop("improper invocation of plregression constructor")

    d = list(
      bws = bws,
      bw = bws,
      xcoef = xcoef,
      xcoeferr = xcoeferr,
      xcoefvcov = xcoefvcov,      
      pregtype = bws$pregtype,
      data.znames = names(evalz),
      data.xnames = names(evalx),
      nobs = dim(evalz)[1],
      zndim = dim(evalz)[2],
      xndim = dim(evalx)[2],
      pscaling = bws$pscaling,
      ptype = bws$ptype,
      pckertype = bws$pckertype,
      pukertype = bws$pukertype,
      pokertype = bws$pokertype,
      evalx = evalx,
      evalz = evalz,
      mean = mean,
      resid = resid,
      ntrain = ntrain,
      trainiseval = trainiseval,
      residuals = residuals,
      R2 = xtra[1],
      MSE = xtra[2],
      MAE = xtra[3],
      MAPE = xtra[4],
      CORR = xtra[5],
      SIGN = xtra[6]
      ,
      timing = timing,
      total.time = total.time,
      optim.time = optim.time,
      fit.time = fit.time
      )

    names(d$xcoeferr) <- names(d$xcoef) <- d$data.xnames
    dimnames(d$xcoefvcov) <- list(d$data.xnames, d$data.xnames)
    
    class(d) = "plregression"
    
    d
  }

print.plregression <- function(x, digits=NULL, ...){
  obj.bws <- .np_plreg_bws(x, where = "print.plregression")

  cat("\nPartially Linear Model",
      "\nRegression data: ", x$ntrain, " training points,",
      if (x$trainiseval) "" else paste(" and ", x$nobs, " evaluation points,",
                                      sep = ""),
      " in ",(x$zndim+x$xndim)," variable(s)",
      "\nWith ", x$xndim, " linear parametric regressor(s), ",
      x$zndim, " nonparametric regressor(s)\n\n", sep="")

  bwmat = matrix(data = 0, nrow = x$xndim+1, ncol = obj.bws$bw$yzbw$ndim)
  
  for (i in seq_along(obj.bws$bw))
    bwmat[i,] = obj.bws$bw[[i]]$bw
  
  print(matrix(bwmat[1,], ncol=x$zndim,
               dimnames=list(paste(x$pscaling,":",sep=""),
                 c("y(z)", replicate(x$zndim-1,"")))))
  cat("\n")
  print(matrix(bwmat[2:(1+x$xndim),], ncol=x$zndim,
               dimnames=list(c(paste(x$pscaling,":",sep=""), replicate(x$xndim-1,"")),
                 c("x(z)", replicate(x$zndim-1,"")))))

  print(matrix(x$xcoef,ncol=x$xndim,
               dimnames=list("Coefficient(s):",x$data.xnames)))
  
  cat(genRegEstStr(x))

  cat(genBwKerStrs(obj.bws))
  cat('\n\n')  

  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

coef.plregression <- function(object, errors = FALSE, ...) {
  if(!errors)
    return(object$xcoef)
  else
    return(object$xcoeferr)
}

vcov.plregression <- function(object,...) {
  return(object$xcoefvcov)
}

fitted.plregression <- function(object, ...){
 object$mean 
}
residuals.plregression <- function(object, ...) {
 if(object$residuals) {
   return(object$resid)
 } else {
   return(npplreg(bws = .np_plreg_bws(object, where = "residuals.plregression"),
                  residuals = TRUE)$resid)
 }
}
predict.plregression <- function(object, se.fit = FALSE, ...) {
  if (isTRUE(se.fit))
    stop("se.fit = TRUE is not implemented for 'plregression' objects", call. = FALSE)

  obj.bws <- .np_plreg_bws(object, where = "predict.plregression")
  dots <- list(...)
  has.formula.route <- !is.null(obj.bws$formula)

  if (!has.formula.route && is.null(dots$exdat) && !is.null(dots$newdata)) {
    nd <- toFrame(dots$newdata)
    need <- c(obj.bws$xnames, obj.bws$znames)
    if (!all(need %in% names(nd)))
      stop("'newdata' must include columns: ", paste(need, collapse = ", "))
    dots$exdat <- nd[, obj.bws$xnames, drop = FALSE]
    dots$ezdat <- nd[, obj.bws$znames, drop = FALSE]
    dots$newdata <- NULL
  }

  tr <- do.call(npplreg, c(list(bws = obj.bws), dots))
  return(fitted(tr))
}
summary.plregression <- function(object, ...){
  obj.bws <- .np_plreg_bws(object, where = "summary.plregression")

  cat("\nPartially Linear Model",
      "\nRegression data: ", object$ntrain, " training points,",
      if (object$trainiseval) "" else paste(" and ", object$nobs, " evaluation points,",
                                      sep = ""),
      " in ",(object$zndim+object$xndim)," variable(s)",
      "\nWith ", object$xndim, " linear parametric regressor(s), ",
      object$zndim, " nonparametric regressor(s)\n", sep="")

  cat(genOmitStr(object))
  cat("\n")

  bwmat = matrix(data = 0, nrow = object$xndim+1, ncol = obj.bws$bw$yzbw$ndim)
  
  for (i in seq_along(obj.bws$bw))
    bwmat[i,] = obj.bws$bw[[i]]$bw
  
  print(matrix(bwmat[1,], ncol=object$zndim,
               dimnames=list(paste(object$pscaling,":",sep=""),
                 c("y(z)", replicate(object$zndim-1,"")))))
  cat("\n")
  print(matrix(bwmat[2:(1+object$xndim),], ncol=object$zndim,
               dimnames=list(c(paste(object$pscaling,":",sep=""), replicate(object$xndim-1,"")),
                 c("x(z)", replicate(object$zndim-1,"")))))
  cat("\n")

  print(matrix(object$xcoef,ncol=object$xndim,
               dimnames=list("Coefficient(s):",object$data.xnames)))
  
  cat(genRegEstStr(object))
  cat("\n")
  cat(genGofStr(object))

  cat(genBwKerStrs(obj.bws))
  cat(genTimingStr(object))
  cat('\n\n')  
}  
