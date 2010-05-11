plregression = 
  function(bws, xcoef, xcoeferr = 0, xcoefvcov, evalx, evalz, mean, resid = NA,
           ntrain, trainiseval = FALSE, residuals = FALSE,
           xtra = double(6)){

    if (missing(bws) | missing(evalx) | missing(evalz) | missing(mean) |
        missing(ntrain) | missing(xcoef))
      stop("improper invocation of plregression constructor")

    d = list(
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
      )

    names(d$xcoeferr) <- names(d$xcoef) <- d$data.xnames
    dimnames(d$xcoefvcov) <- list(d$data.xnames, d$data.xnames)
    
    class(d) = "plregression"
    
    d
  }

print.plregression <- function(x, digits=NULL, ...){
  cat("\nPartially Linear Model",
      "\nRegression data: ", x$ntrain, " training points,",
      ifelse(x$trainiseval, "", paste(" and ", x$nobs, " evaluation points,",
                                      sep = "")),
      " in ",(x$zndim+x$xndim)," variable(s)",
      "\nWith ", x$xndim, " linear parametric regressor(s), ",
      x$zndim, " nonparametric regressor(s)\n\n", sep="")

  bwmat = matrix(data = 0, nrow = x$xndim+1, ncol = x$bw$bw$yzbw$ndim)
  
  for (i in 1:length(x$bw$bw))
    bwmat[i,] = x$bw$bw[[i]]$bw
  
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

  cat(genBwKerStrs(x$bw))
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
 if(object$residuals) { return(object$resid) } else { return(npplreg(bws = object$bw, residuals =TRUE)$resid) } 
}
predict.plregression <- function(object, se.fit = FALSE, ...) {
  tr <- eval(npplreg(bws = object$bw, ...), env = parent.frame())
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs, residual.scale = tr$MSE))
  else
    return(fitted(tr))
}
plot.plregression <- function(x, ...) { npplot(bws = x$bw, ...) }
summary.plregression <- function(object, ...){
  cat("\nPartially Linear Model",
      "\nRegression data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "", paste(" and ", object$nobs, " evaluation points,",
                                      sep = "")),
      " in ",(object$zndim+object$xndim)," variable(s)",
      "\nWith ", object$xndim, " linear parametric regressor(s), ",
      object$zndim, " nonparametric regressor(s)\n", sep="")

  cat(genOmitStr(object))
  cat("\n")

  bwmat = matrix(data = 0, nrow = object$xndim+1, ncol = object$bw$bw$yzbw$ndim)
  
  for (i in 1:length(object$bw$bw))
    bwmat[i,] = object$bw$bw[[i]]$bw
  
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

  cat(genBwKerStrs(object$bw))
  cat('\n\n')  
}  
