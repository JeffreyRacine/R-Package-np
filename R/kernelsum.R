npkernelsum = 
  function(bws, eval, ksum,
           ntrain, trainiseval = FALSE){

    if (missing(bws) | missing(eval) | missing(ksum) | missing(ntrain))
      stop("improper invocation of npkernelsum constructor")

    d = list(
      bw = bws$bw,
      data.names = bws$data.names,
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
      ksum = ksum,
      ntrain = ntrain,
      trainiseval = trainiseval
      )

    class(d) = "npkernelsum"

    d
  }

print.npkernelsum <- function(x, digits=NULL, ...){
  cat("\nKernel Sum data: ", x$ntrain, " training points,",
      ifelse(x$trainiseval, "", " and "),
      ifelse(x$trainiseval, "", x$nobs),
      ifelse(x$trainiseval, ""," evaluation points,"),
      " in ",x$ndim," variable(s)\n",sep="")

  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$data.names)))
  
  cat("\nKernel Estimator:",x$pregtype,
      "\nBandwidth Type:",x$ptype,
      "\n\nKernel sum in component 'ksum'")

  if (x$ncon > 0)
    cat("\n\nContinuous Kernel Type:",x$pckertype,
        "\nNo. Continuous Vars.:",x$ncon)

  if (x$nuno > 0)
    cat("\n\nUnordered Categorical Kernel Type:",x$pukertype,
        "\nNo. Unordered Categorical Vars.:",x$nuno)

  if (x$nord > 0)
    cat("\n\nOrdered Categorical Kernel Type:",x$pokertype,
        "\nNo. Ordered Categorical Vars.:",x$nord)

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

