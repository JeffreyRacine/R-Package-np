rbandwidth <-
  function(bw = stop("rbandwidth:argument 'bw' missing"),
           regtype = c("lc","ll"),
           bwmethod = c("cv.ls","cv.aic"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian", "epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ukertype = c("aitchisonaitken", "liracine"),
           okertype = c("wangvanryzin", "liracine"),
           fval = NA,
           ifval = NA,
           nobs = NA,
           xdati = stop("rbandwidth:argument 'xdati' missing"),
           ydati = stop("rbandwidth:argument 'ydati' missing"),
           xnames = character(length(bw)),
           ynames = character(1),
           sfactor = NA, bandwidth = NA,
           rows.omit = NA,
           bandwidth.compute = TRUE,...){

  ndim = length(bw)
  regtype = match.arg(regtype)
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

  ukertype = match.arg(ukertype)
  okertype = match.arg(okertype)

  porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )
  ## calculate some info to be pretty-printed

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
    regtype = regtype,
    pregtype = switch(regtype,
      lc = "Local-Constant",
      ll = "Local-Linear"),
    method = bwmethod,
    pmethod = switch( bwmethod,
      cv.ls = "Least Squares Cross-Validation",
      cv.aic = "Expected Kullback-Leibler Cross-Validation"
##      ,"normal-reference" = "Normal Reference"
      ),
    fval = fval,
    ifval = ifval,
    scaling = bwscaling,
    pscaling = ifelse(bwscaling, "Scale Factor(s)", "Bandwidth(s)"),
    type = bwtype,
    ptype = switch( bwtype,
      fixed = "Fixed",
      generalized_nn = "Generalized Nearest Neighbour",
      adaptive_nn = "Adaptive Nearest Neighbour" ),
    ckertype = ckertype,    
    ckerorder = ckerorder,
    pckertype = switch(ckertype,
      gaussian = paste(porder,"Gaussian"),
      epanechnikov =  paste(porder,"Epanechnikov"),
      uniform = "Uniform"),
    ukertype = ukertype,
    pukertype = switch( ukertype,
      aitchisonaitken = "Aitchison and Aitken",
      liracine = "Li and Racine"),
    okertype = okertype,
    pokertype = switch( okertype,
      wangvanryzin = "Wang and Van Ryzin",
      liracine = "Li and Racine"),
    nobs = nobs,
    ndim = ndim,
    ncon = sum(xdati$icon),
    nuno = sum(xdati$iuno),
    nord = sum(xdati$iord),
    icon = xdati$icon,
    iuno = xdati$iuno,
    iord = xdati$iord,
    xnames = xnames,
    ynames = ynames,
    xdati = xdati,
    ydati = ydati,
    xmcv = mcvConstruct(xdati),
    sfactor = list(x = sfactor),
    bandwidth = list(x = bandwidth),
    sumNum = list(x = sumNum),
    dati = list(x = xdati, y = ydati),
    varnames = list(x = xnames, y = ynames),
    vartitle = list(x = "Explanatory", y = "Dependent"),
    vartitleabb = list(x = "Exp.", y = "Dep."),
    rows.omit = rows.omit,
    nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))

  
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

  class(mybw) = "rbandwidth"
  if(!any(is.na(mybw$bandwidth)))
    validateBandwidth(mybw)
  mybw
  
}

as.double.rbandwidth <- function(x, ...){
  x$bw
}

## feature: when using dataframe interface, summary and print methods don't 
## provide info on the dependent variable
print.rbandwidth <- function(x, digits=NULL, ...){
  cat("\nRegression Data (",x$nobs," observations, ",x$ndim," variable(s)):\n\n",sep="")
  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$xnames)))

  cat(genBwSelStr(x))
  cat(genBwKerStrs(x))

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}


plot.rbandwidth <- function(...) { npplot(...)  }
predict.rbandwidth <- function(...) { eval(npreg(...), envir = parent.frame()) }

summary.rbandwidth <- function(object, ...){
  cat("\nRegression Data (",object$nobs," observations, ",object$ndim," variable(s)):\n",sep="")

  cat(genOmitStr(object))
  cat(genBwSelStr(object))

  cat('\n')
  cat(genBwScaleStrs(object))
  cat(genBwKerStrs(object))

  cat("\n\n")

}
