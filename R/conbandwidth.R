conbandwidth <-
  function(xbw,
           ybw,
           bwmethod = c("cv.ml","cv.ls","normal-reference", "cv.ls.np", "manual"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           cxkertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           cxkerorder = c(2,4,6,8),
           uxkertype = c("aitchisonaitken","liracine"),
           oxkertype = c("liracine","wangvanryzin"),
           cykertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           cykerorder = c(2,4,6,8),
           uykertype = c("aitchisonaitken","liracine"),
           oykertype = c("liracine","wangvanryzin"),
           fval = NA,
           ifval = NA,
           fval.history = NA,
           nobs = NA,
           xdati, ydati,
           xnames = character(length(xbw)),
           ynames = character(length(ybw)),
           sfactor = NA, bandwidth = NA,
           rows.omit = NA,
           nconfac = NA,
           ncatfac = NA,
           sdev = NA,
           bandwidth.compute = TRUE,
           timing = NA,
           total.time = NA,
           ...){

  if (missing(xbw) | missing(ybw))
    stop("improper invocation of conbandwidth constructor: 'bw' or i[cuo]* missing")
  
  xndim = length(xbw)
  yndim = length(ybw)
  
  bwmethod = match.arg(bwmethod)
  bwtype = match.arg(bwtype)

  cxkertype = match.arg(cxkertype)
  cykertype = match.arg(cykertype)

  if(missing(cxkerorder))
    cxkerorder = 2
  else if (cxkertype == "uniform")
    warning("ignoring kernel order specified with uniform kernel type")
  else {
    kord = eval(formals()$cxkerorder) 
    if (!any(kord == cxkerorder))
      stop("cxkerorder must be one of ", paste(kord,collapse=" "))
  }

  if (cxkertype == "truncated gaussian" && cxkerorder != 2)
    warning("using truncated gaussian of order 2, higher orders not yet implemented")

  if (bwmethod == "normal-reference" && (cxkertype != "gaussian" || bwtype != "fixed")){    
    warning("normal-reference bandwidth selection assumes gaussian kernel with fixed bandwidth")
    bwtype = "fixed"
    cxkertype = "gaussian"
  }

  if(missing(cykerorder))
    cykerorder = 2
  else if (cykertype == "uniform")
    warning("ignoring kernel order specified with uniform kernel type")
  else {
    kord = eval(formals()$cykerorder) 
    if (!any(kord == cykerorder))
      stop("cykerorder must be one of ", paste(kord,collapse=" "))
  }

  if (cykertype == "truncated gaussian" && cykerorder != 2)
    warning("using truncated gaussian of order 2, higher orders not yet implemented")

  if (bwmethod == "normal-reference" && (cykertype != "gaussian" || bwtype != "fixed")){    
    warning("normal-reference bandwidth selection assumes gaussian kernel with fixed bandwidth")
    bwtype = "fixed"
    cykertype = "gaussian"
  }

  if (cxkerorder != cykerorder & bwscaling)
    stop("scale factors with different order kernels for dependent and explanatory variables is unsupported")
  
  uxkertype = match.arg(uxkertype)
  uykertype = match.arg(uykertype)
  
  oxkertype = match.arg(oxkertype)
  oykertype = match.arg(oykertype)

  pxorder = switch( cxkerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )
  pyorder = switch( cykerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )

  dati <- list(x = xdati, y = ydati)
  
  if (!identical(sfactor,NA)){
    ## using the new model for generically accessing bandwidth objects

    okertype <- list(x = oxkertype, y = oykertype)
    ukertype <- list(x = uxkertype, y = uykertype)

    scaleOrMax <- function(i, j) {
      if (dati[[j]]$icon[i])
        return((sfactor[[j]])[i])

      if (dati[[j]]$iord[i])
        return(oMaxL(dati[[j]]$all.nlev[[i]], kertype = okertype[[j]]))
      
      if (dati[[j]]$iuno[i])
        return(uMaxL(dati[[j]]$all.nlev[[i]], kertype = ukertype[[j]]))
    }

    sumNum <- list(x = NA, y = NA)
    sumNum[] <- lapply(1:length(dati), function(i) {
      sapply(1:length(dati[[i]]$icon), scaleOrMax, j = i)
    })
  } else {
    sumNum <- NA
  }

  if (length(rows.omit) == 0)
    rows.omit <- NA

  mybw = list(
    xbw=xbw,
    ybw=ybw,
    method = bwmethod,
    pmethod = bwmToPrint(bwmethod),
    fval = fval,
    ifval = ifval,
    fval.history = fval.history,
    scaling = bwscaling,
    pscaling = ifelse(bwscaling, "Scale Factor(s)", "Bandwidth(s)"),
    type = bwtype,
    ptype = bwtToPrint(bwtype),
    cxkertype = cxkertype,
    cykertype = cykertype,
    cxkerorder = cxkerorder,
    cykerorder = cykerorder,
    pcxkertype = cktToPrint(cxkertype, order = pxorder),
    pcykertype = cktToPrint(cykertype, order = pyorder),
    uxkertype = uxkertype,
    uykertype = uykertype,
    puxkertype = uktToPrint(uxkertype),
    puykertype = uktToPrint(uykertype),
    oxkertype = oxkertype,
    oykertype = oykertype,
    poxkertype = oktToPrint(oxkertype),
    poykertype = oktToPrint(oykertype, normalized = TRUE),
    nobs = nobs,
    xndim = xndim,
    yndim = yndim,
    ndim = xndim + yndim,
    xncon = sum(xdati$icon),
    xnuno = sum(xdati$iuno),
    xnord = sum(xdati$iord),
    yncon = sum(ydati$icon),
    ynuno = sum(ydati$iuno),
    ynord = sum(ydati$iord),
    ncon = sum(c(xdati$icon, ydati$icon)),
    ixcon = xdati$icon,
    ixuno = xdati$iuno,
    ixord = xdati$iord,
    iycon = ydati$icon,
    iyuno = ydati$iuno,
    iyord = ydati$iord,
    xnames = xnames,
    ynames = ynames,
    xdati = xdati,
    ydati = ydati,
    xmcv = mcvConstruct(xdati),
    ymcv = mcvConstruct(ydati),
    sfactor = sfactor,
    bandwidth = bandwidth,
    nconfac = nconfac,
    ncatfac = ncatfac,
    sdev = sdev,
    sumNum = sumNum,
    dati = dati, 
    varnames = list(x = xnames, y = ynames),
    vartitle = list(x = "Explanatory", y = "Dependent"),
    vartitleabb = list(x = "Exp.", y = "Dep."),
    rows.omit = rows.omit,
    nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)),
    timing = timing,
    total.time = total.time)

  mybw$klist = list(
    x =
    list(ckertype = cxkertype,
         pckertype = mybw$pcxkertype,
         ukertype = uxkertype,
         pukertype = mybw$puxkertype,
         okertype = oxkertype,
         pokertype = mybw$poxkertype),
    y =
    list(ckertype = cykertype,
         pckertype = mybw$pcykertype,
         ukertype = uykertype,
         pukertype = mybw$puykertype,
         okertype = oykertype,
         pokertype = mybw$poykertype))

  if(!bandwidth.compute)
    mybw$pmethod <- "Manual"


  class(mybw) = "conbandwidth"
  if(!any(is.na(mybw$bandwidth)))
    validateBandwidth(mybw)
  mybw
}

print.conbandwidth <- function(x, digits=NULL, ...){
  cat("\nConditional density data (",x$nobs," observations, ",
      (x$xndim+x$yndim)," variable(s))",
      "\n(", x$yndim, " dependent variable(s), and ", x$xndim, " explanatory variable(s))\n\n",
      sep="")
  print(matrix(x$ybw,ncol=x$yndim,dimnames=list(paste("Dep. Var. ",x$pscaling,":",sep=""),x$ynames)))

  print(matrix(x$xbw,ncol=x$xndim,dimnames=list(paste("Exp. Var. ",x$pscaling,":",sep=""),x$xnames)))

  cat(genBwSelStr(x))
  cat(genBwKerStrsXY(x))
  
  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

plot.conbandwidth <- function(...) { npplot(...) }

summary.conbandwidth <- function(object, ...) {
  cat("\nConditional density data (",object$nobs," observations, ",
      (object$xndim+object$yndim)," variable(s))",
      "\n(", object$yndim, " dependent variable(s), and ", object$xndim, " explanatory variable(s))\n",
      sep="")

  cat(genOmitStr(object))
  cat(genBwSelStr(object))

  cat(paste("\n", genBwScaleStrs(object), sep=""))
  cat(genBwKerStrs(object))

  cat(genTimingStr(object))
  
  cat("\n\n")
}

predict.conbandwidth <- function(...) { eval(npcdens(...), envir = parent.frame()) }
