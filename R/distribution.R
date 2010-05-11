npdistribution =
  function(bws, eval, dist, derr = NA,
           ntrain, trainiseval = FALSE,
           rows.omit = NA){

    if (missing(bws) | missing(eval) | missing(dist) | missing(ntrain))
      stop("improper invocation of distribution constructor")

    if (length(rows.omit) == 0)
      rows.omit <- NA

    d = list(
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
      ntrain = ntrain,
      trainiseval = trainiseval,
      rows.omit = rows.omit,
      nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))


    class(d) = "npdistribution"

    d
  }

print.npdistribution <- function(x, digits=NULL, ...){
  cat("\nDistribution Data: ", x$ntrain, " training points,",
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

fitted.npdistribution <- function(object, ...){
 object$dist 
}
se.npdistribution <- function(x){ x$derr }
plot.npdistribution <- function(x, ...) { npplot(bws = x$bws, ...) }

predict.npdistribution <- function(object, se.fit = FALSE, ...) {
  tr <- eval(npudist(bws = object$bws, ...), env = parent.frame())
  if(se.fit)
    return(list(fit = fitted(tr), se.fit = se(tr), 
                df = tr$nobs))
  else
    return(fitted(tr))
}

summary.npdistribution <- function(object, ...) {
  cat("\nDistribution Data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "", paste(" and ", object$nobs,
                                      " evaluation points,", sep="")),
      " in ",object$ndim," variable(s)\n",sep="")
  
  cat(genOmitStr(object))

  print(matrix(object$bw,ncol=object$ndim,dimnames=list(paste(object$pscaling,":",sep=""),object$xnames)))

  cat(genDenEstStr(object))

  cat(genBwKerStrs(object$bws))
  cat('\n\n')

}
