condistribution =
  function(bws, xeval, yeval, condist, conderr = NA,
           congrad = NA, congerr = NA,           
           ntrain, trainiseval = FALSE,
           rows.omit = NA){

    if (missing(bws) | missing(xeval) | missing(yeval) | missing(condist) | missing(ntrain))
      stop("improper invocation of condistribution constructor")

    if (length(rows.omit) == 0)
      rows.omit <- NA

    d = list(
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
      rows.omit = rows.omit,
      nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))


    class(d) = "condistribution"

    d
  }

print.condistribution <- function(x, digits=NULL, ...){
  cat("\nConditional distribution data: ", x$ntrain, " training points,",
      ifelse(x$trainiseval, "", paste(" and ", x$nobs, " evaluation points,\n", sep="")),
      " in ", x$xndim + x$yndim, " variable(s)",
      "\n(", x$yndim, " dependent variable(s), and ", x$xndim, " explanatory variable(s))\n\n",
      sep="")
  print(matrix(x$ybw,ncol=x$yndim,dimnames=list(paste("Dep. Var. ",x$pscaling,":",sep=""),x$ynames)))

  print(matrix(x$xbw,ncol=x$xndim,dimnames=list(paste("Exp. Var. ",x$pscaling,":",sep=""),x$xnames)))

  cat(genDenEstStr(x))
  cat(genBwKerStrs(x$bws))


  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

fitted.condistribution <- function(object, ...){
 object$condist 
}
se.condistribution <- function(x){ x$conderr }
gradients.condistribution <- function(x, errors = FALSE, ...) {
  if(!errors)
    return(x$congrad)
  else
    return(x$congerr)
}

predict.condistribution <- function(object, se.fit = FALSE, ...) {
  tr <- eval(npcdist(bws = object$bws, ...), env = parent.frame())
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs))
  else
    return(fitted(tr))
}

plot.condistribution <- function(x, ...) { npplot(bws = x$bws, ...) }

summary.condistribution <- function(object, ...){
  cat("\nConditional Distribution Data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "", paste(" and ", object$nobs, " evaluation points,\n", sep="")),
      " in ", object$xndim + object$yndim, " variable(s)",
      "\n(", object$yndim, " dependent variable(s), and ", object$xndim, " explanatory variable(s))\n\n",
      sep="")

  cat(genOmitStr(object))
  print(matrix(object$ybw,ncol=object$yndim,dimnames=list(paste("Dep. Var. ",object$pscaling,":",sep=""),object$ynames)))

  print(matrix(object$xbw,ncol=object$xndim,dimnames=list(paste("Exp. Var. ",object$pscaling,":",sep=""),object$xnames)))

  cat(genDenEstStr(object))

  cat(genBwKerStrs(object$bws))
  cat('\n\n')  
}
