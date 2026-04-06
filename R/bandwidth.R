bandwidth <-
  function(bw = stop("bandwidth: argument 'bw' missing"),
           bwmethod = c("cv.ml","cv.ls","normal-reference"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ckerbound = c("none","range","fixed"),
           ckerlb = NULL,
           ckerub = NULL,
           ukertype = c("aitchisonaitken","liracine"),
           okertype = c("liracine","wangvanryzin","racineliyan"),
           fval = NA,
           ifval = NA,
           num.feval = NA,
           num.feval.fast = NA,
           fval.history = NA,
           eval.history = NA,
           invalid.history = NA,
           nobs = NA,
           xdati = stop("bandwidth:argument 'xdati' missing"),
           xnames = character(length(bw)),
           sfactor = NA, bandwidth = NA,
           rows.omit = NA, 
           nconfac = NA, ncatfac = NA, sdev = NA,
           bandwidth.compute = TRUE,
           timing = NA,
           total.time = NA,
           ...){
    
    ndim = length(bw)
    bwmethod = match.arg(bwmethod)
    bwtype = match.arg(bwtype)
    ckertype = match.arg(ckertype)
    ckerbound = match.arg(ckerbound)

    if(missing(ckerorder))
      ckerorder = 2
    else if (ckertype == "uniform")
      .np_warning("ignoring kernel order specified with uniform kernel type")
    else {
      kord = c(2,4,6,8) 
      if (!any(kord == ckerorder))
        stop("ckerorder must be one of ", paste(kord,collapse=" "))
    }

    if (ckertype == "truncated gaussian" && ckerorder != 2)
      .np_warning("using truncated gaussian of order 2, higher orders not yet implemented")

    if (bwmethod == "normal-reference" && (ckertype != "gaussian" || bwtype != "fixed")){    
      .np_warning("normal-reference bandwidth selection assumes gaussian kernel with fixed bandwidth")
      bwtype = "fixed"
      ckertype = "gaussian"
    }

    ukertype = match.arg(ukertype)
    okertype = match.arg(okertype)
    cbounds <- npKernelBoundsResolve(
      dati = xdati,
      varnames = xnames,
      kerbound = ckerbound,
      kerlb = ckerlb,
      kerub = ckerub,
      argprefix = "cker")
    bounded_nonfixed_supported <- bwtype %in% c("generalized_nn", "adaptive_nn")
    if (bwtype != "fixed" &&
        cbounds$bound != "none" &&
        !bounded_nonfixed_supported)
      stop("finite continuous kernel bounds require bwtype = \"fixed\"")

    porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )

    if (!identical(sfactor,NA)){
      sumNum <- sapply(seq_len(ndim), function(i) {
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
      num.feval = num.feval,
      num.feval.fast = num.feval.fast,
      fval.history = fval.history,
      eval.history = eval.history,
      invalid.history = invalid.history,
      scaling = bwscaling,
      pscaling = npBandwidthSummaryLabel(bwtype = bwtype, bwscaling = bwscaling),
      type = bwtype,
      ptype = bwtToPrint(bwtype),
      ckertype = ckertype,    
      ckerorder = ckerorder,
      ckerbound = cbounds$bound,
      ckerlb = cbounds$lb,
      ckerub = cbounds$ub,
      pckertype = cktToPrint(ckertype, order = porder, kerbound = cbounds$bound),
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
      nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
      timing = timing,
      total.time = total.time)

    mybw$klist = list(
      x =
      list(ckertype = ckertype,
           ckerbound = cbounds$bound,
           ckerlb = cbounds$lb,
           ckerub = cbounds$ub,
           pckertype = mybw$pckertype,
           ukertype = ukertype,
           pukertype = mybw$pukertype,
           okertype = okertype,
           pokertype = mybw$pokertype))

    if(!bandwidth.compute)
      mybw$pmethod <- "Manual"
    

    class(mybw) = "bandwidth"
    if(!any(is.na(mybw$bandwidth)))
    validateBandwidth(mybw)
    mybw
    
  }

as.double.bandwidth <- function(x, ...){
  x$bw
}

print.bandwidth <- function(x, digits=NULL, ...){
  cat("\nData (",x$nobs," observations, ",x$ndim," variable(s)):\n\n",sep="")
  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$xnames)))

  cat(genBwSelStr(x))

  cat(genBwKerStrs(x))

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

summary.bandwidth <- function(object, ...) {
  cat("\nData (",object$nobs," observations, ",object$ndim," variable(s)):\n",sep="")

  cat(genOmitStr(object))
  
  cat(genBwSelStr(object))

  cat('\n')
  cat(genBwScaleStrs(object))

  cat(genBwKerStrs(object))

  cat(genTimingStr(object))
  
  cat("\n\n")
}


predict.bandwidth <- function(...) { do.call(npudens, list(...)) }
