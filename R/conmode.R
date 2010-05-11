conmode =
  function(bws, xeval, yeval = NA,
           conmode,
           condens, conderr = NA,
           confusion.matrix = NA,
           CCR.overall = NA,
           CCR.byoutcome = NA,
           fit.mcfadden = NA,
           ntrain, trainiseval = FALSE){

    if (missing(bws) | missing(xeval) | missing(conmode) |
        missing(condens) | missing(ntrain))
      stop("improper invocation of conmode constructor")

    d = list(
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
      conmode = conmode,
      condens = condens,
      conderr = conderr,
      confusion.matrix = confusion.matrix,
      CCR.overall = CCR.overall,
      CCR.byoutcome = CCR.byoutcome,
      fit.mcfadden = fit.mcfadden,
      ntrain = ntrain,
      trainiseval = trainiseval)

    class(d) = "conmode"

    d
  }


print.conmode <- function(x, ...){
  cat("\nConditional Mode data: ", x$ntrain, " training points,",
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
    print(...)
  invisible(x)
}

fitted.conmode <- function(object, ...) {
 object$condens 
}
mode.conmode <- function(x) { x$conmode }

summary.conmode <- function(object, ...){
  cat("\nConditional Mode data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "", paste(" and ", object$nobs, " evaluation points,\n", sep="")),
      " in ", object$xndim + object$yndim, " variable(s)",
      "\n(", object$yndim, " dependent variable(s), and ", object$xndim, " explanatory variable(s))\n\n",
      sep="")

  cat(genOmitStr(object))
  print(matrix(object$ybw,ncol=object$yndim,dimnames=list(paste("Dep. Var. ",object$pscaling,":",sep=""),object$ynames)))

  print(matrix(object$xbw,ncol=object$xndim,dimnames=list(paste("Exp. Var. ",object$pscaling,":",sep=""),object$xnames)))

  cat(genDenEstStr(object))
  cat('\n')
  pCatGofStr(object)

  cat(genBwKerStrs(object$bws))
  cat('\n\n')
  
}
