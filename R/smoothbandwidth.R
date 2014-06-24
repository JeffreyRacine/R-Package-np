scbandwidth <-
  function(bw = stop("scbandwidth:argument 'bw' missing"),
           bwmethod = c("cv.ls", "manual"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ukertype = c("aitchisonaitken", "liracine"),
           okertype = c("liracine","wangvanryzin"),
           fval = NA,
           ifval = NA,
           nobs = NA,
           numimp = NA,
           fval.vector = NA,
           xdati = stop("scbandwidth:argument 'xdati' missing"),
           ydati = stop("scbandwidth:argument 'ydati' missing"),
           zdati = NULL,
           xnames = character(length(bw)),
           ynames = character(1),
           znames = NULL,
           sfactor = NA, bandwidth = NA,
           rows.omit = NA,
           nconfac = NA,
           ncatfac = NA,
           sdev = NA,
           bandwidth.compute = TRUE,
           optim.method = "NA",
           total.time = NA,
           ...){

  ndim = length(bw)
  bwmethod = match.arg(bwmethod)
  bwtype = match.arg(bwtype)
  ckertype = match.arg(ckertype)

  if(missing(ckerorder))
    ckerorder = 2
  else if (ckertype == "uniform")
    warning("ignoring kernel order specified with uniform kernel type")
  else {
    kord = eval(formals()$ckerorder) 
    if (!any(kord == ckerorder))
      stop("ckerorder must be one of ", paste(kord,collapse=" "))
  }

  if (ckertype == "truncated gaussian" && ckerorder != 2)
    warning("using truncated gaussian of order 2, higher orders not yet implemented")

  ukertype = match.arg(ukertype)
  okertype = match.arg(okertype)

  tdati <- zdati

  if(is.null(zdati))
    tdati <- xdati

  porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )

  if (!identical(sfactor,NA)){
    sumNum <- sapply(1:ndim, function(i) {
      if (tdati$icon[i])
        return(sfactor[i])

      if (tdati$iord[i])
        return(oMaxL(tdati$all.nlev[[i]], kertype = okertype))
      
      if (tdati$iuno[i])
        return(uMaxL(tdati$all.nlev[[i]], kertype = ukertype))
    })
  } else {
    sumNum <- NA
  }

  if (length(rows.omit) == 0)
    rows.omit <- NA

  
  mybw = list(
    bw=bw,
    method = bwmethod,
    pmethod = bwmToPrint(bwmethod),
    pomethod = switch(optim.method,
      "Nelder-Mead" = "Nelder-Mead",
      "BFGS" = "BFGS",
      "CG" = "CG", "NA"),
    fval = fval,
    ifval = ifval,
    scaling = bwscaling,
    pscaling = ifelse(bwscaling, "Scale Factor(s)", "Bandwidth(s)"),
    type = bwtype,
    ptype = bwtToPrint(bwtype),
    ckertype = ckertype,    
    ckerorder = ckerorder,
    pckertype = cktToPrint(ckertype, order = porder),
    ukertype = ukertype,
    pukertype = uktToPrint(ukertype),
    okertype = okertype,
    pokertype = oktToPrint(okertype),
    nobs = nobs,
    ndim = ndim,
    ncon = sum(tdati$icon),
    nuno = sum(tdati$iuno),
    nord = sum(tdati$iord),
    icon = tdati$icon,
    iuno = tdati$iuno,
    iord = tdati$iord,
    xnames = xnames,
    ynames = ynames,
    znames = znames,
    xdati = xdati,
    ydati = ydati,
    zdati = zdati,
    xmcv = mcvConstruct(xdati),
    sfactor = list(sfactor),
    bandwidth = list(bandwidth),
    sumNum = list(sumNum),
    dati = list(x = xdati, y = ydati, z = zdati),
    varnames = list(x = xnames, y = ynames, z = znames),
    vartitle = list(x = "Explanatory", y = "Dependent", z = "Explanatory"),
    vartitleabb = list(x = "Exp.", y = "Dep.", z = "Exp."),
    rows.omit = rows.omit,
    nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)),
    total.time = total.time)

  mybw$klist <-
      list(list(ckertype = ckertype,
                pckertype = mybw$pckertype,
                ukertype = ukertype,
                pukertype = mybw$pukertype,
                okertype = okertype,
                pokertype = mybw$pokertype))

  zorx <- ifelse(is.null(zdati), "x", "z")
  names(mybw$sfactor) <- zorx
  names(mybw$bandwidth) <- zorx 
  names(mybw$sumNum) <- zorx
  names(mybw$klist) <- zorx
  
  if(!bandwidth.compute)
    mybw$pmethod <- "Manual"

  class(mybw) = "scbandwidth"
  if(!any(is.na(mybw$bandwidth)))
    validateBandwidth(mybw)
  mybw
  
}

as.double.scbandwidth <- function(x, ...){
  x$bw
}

print.scbandwidth <- function(x, digits=NULL, ...){
  cat("\nSmooth Coefficient Model",
      "\nRegression Data (",x$nobs," observations, ",x$ndim," variable(s)):\n\n",sep="")
  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),
                                  if(is.null(x$znames)) x$xnames else x$znames)))

  cat(genBwSelStr(x))
  cat(genBwKerStrs(x))
  
  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

plot.scbandwidth <- function(...) { npplot(...)  }
predict.scbandwidth <- function(...) { eval(npscoef(...), envir = parent.frame()) }

summary.scbandwidth <- function(object, ...){
  cat("\nSmooth Coefficient Regression",
      "\nRegression Data (",object$nobs," observations, ",object$ndim," variable(s)):\n",sep="")

  cat("\nOptimisation Method: ", object$pomethod)

  cat(genOmitStr(object))
  cat(genBwSelStr(object))

  cat('\n')
  cat(genBwScaleStrs(object))

  cat(genBwKerStrs(object))
  cat(genTimingStr(object))
  cat("\n\n")

}
