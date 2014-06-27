dbandwidth <-
  function(bw = stop("dbandwidth: argument 'bw' missing"),
           bwmethod = c("cv.cdf","normal-reference"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ukertype = c("aitchisonaitken","liracine"),
           okertype = c("liracine","wangvanryzin"),
           fval = NA,
           ifval = NA,
           fval.history = NA,
           nobs = NA,
           xdati = stop("dbandwidth:argument 'xdati' missing"),
           xnames = character(length(bw)),
           sfactor = NA, bandwidth = NA,
           rows.omit = NA,
           nconfac = NA,
           ncatfac = NA,
           sdev = NA,
           bandwidth.compute = TRUE,
           timing = NA,
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

    if (bwmethod == "normal-reference" && (ckertype != "gaussian" || bwtype != "fixed")){    
      warning("normal-reference bandwidth selection assumes gaussian kernel with fixed bandwidth")
      bwtype = "fixed"
      ckertype = "gaussian"
    }

    if (ckertype == "truncated gaussian" && ckerorder != 2)
      warning("using truncated gaussian of order 2, higher orders not yet implemented")


    ukertype = match.arg(ukertype)
    okertype = match.arg(okertype)

    porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )

    if (!identical(sfactor,NA)){
      sumNum <- sapply(1:ndim, function(i) {
        if (xdati$icon[i])
          return(sfactor[i])

        if (xdati$iord[i])
          return(oMaxL(xdati$all.nlev[[i]], kertype = okertype))
        
        if (xdati$iuno[i])
          return(uMaxL(xdati$all.nlev[[i]], kertype = ukertype))
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
      fval = fval,
      ifval = ifval,
      fval.history = fval.history,
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
      pokertype = oktToPrint(okertype, normalized = TRUE),
      nobs = nobs,
      ndim = ndim,
      ncon = sum(xdati$icon),
      nuno = sum(xdati$iuno),
      nord = sum(xdati$iord),
      icon = xdati$icon,
      iuno = xdati$iuno,
      iord = xdati$iord,
      xnames = xnames,
      xdati = xdati,
      sfactor = list(x = sfactor),
      bandwidth = list(x = bandwidth),
      nconfac = nconfac,
      ncatfac = ncatfac,
      sdev = sdev,
      sumNum = list(x = sumNum),
      xmcv = mcvConstruct(xdati),
      dati = list(x = xdati),
      varnames = list(x = xnames),
      vartitle = list(x = ""),
      vartitleabb = list(x = ""),
      rows.omit = rows.omit,
      nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)),
      timing = timing,
      total.time = total.time)

    mybw$klist = list(
      x =
      list(ckertype = ckertype,
           pckertype = mybw$pckertype,
           ukertype = ukertype,
           pukertype = mybw$pukertype,
           okertype = okertype,
           pokertype = mybw$pokertype))

    if(!bandwidth.compute)
      mybw$pmethod <- "Manual"
    

    class(mybw) = "dbandwidth"
    if(!any(is.na(mybw$bandwidth)))
    validateBandwidth(mybw)
    mybw
    
  }

as.double.dbandwidth <- function(x, ...){
  x$bw
}

print.dbandwidth <- function(x, digits=NULL, ...){
  cat("\nData (",x$nobs," observations, ",x$ndim," variable(s)):\n\n",sep="")
  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$xnames)))

  cat(genBwSelStr(x))

  cat(genBwKerStrs(x))

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

summary.dbandwidth <- function(object, ...) {
  cat("\nData (",object$nobs," observations, ",object$ndim," variable(s)):\n",sep="")

  cat(genOmitStr(object))
  
  cat(genBwSelStr(object))

  cat('\n')
  cat(genBwScaleStrs(object))

  cat(genBwKerStrs(object))

  cat(genTimingStr(object))
  
  cat("\n\n")
}

plot.dbandwidth <- function(...) { npplot(...)  }

predict.dbandwidth <- function(...) { eval(npudist(...), envir =parent.frame()) }
