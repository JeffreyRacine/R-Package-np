sibandwidth <-
  function(beta, h,
           method=c("ichimura","kleinspady"),
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian", "epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           fval = NA,
           ifval = NA,
           numimp = NA,
           fval.vector = NA,
           nobs = NA,
           xdati, ydati, idati = untangle(data.frame(double(1))),
           xnames = character(length(beta)),
           ynames = character(1),
           sfactor = NA, bandwidth = NA,
           rows.omit = NA, bandwidth.compute = TRUE,
           optim.method = "NA",
           only.optimize.beta = FALSE,
           ...){

  ndim = length(beta)
  regtype = "lc"
  method = match.arg(method)
  ckertype = match.arg(ckertype)
  bwtype <- match.arg(bwtype)
  
  if(missing(ckerorder))
    ckerorder = 2
  else if (ckertype == "uniform")
    warning("ignoring kernel order specified with uniform kernel type")
  else {
    kord = eval(formals()$ckerorder) 
    if (!any(kord == ckerorder))
      stop("ckerorder must be one of ", paste(kord,collapse=" "))
  }

  porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )

  sumNum <- sfactor
  ##idati <- NA
  
  if (length(rows.omit) == 0)
    rows.omit <- NA
  
  mybw = list(
    beta=beta,
    bw = h,
    regtype = regtype,
    pregtype = switch(regtype,
      lc = "Local-Constant",
      ll = "Local-Linear"),
    method = method,
    pmethod = switch( method,
      ichimura = "Ichimura",
      "kleinspady" = "Klein and Spady"
      ),
    pomethod = switch(optim.method,
      "Nelder-Mead" = "Nelder-Mead",
      "BFGS" = "BFGS",
      "CG" = "CG", "NA"),
    fval = fval,
    ifval = ifval,
    numimp = numimp,
    fval.vector = fval.vector,
    pscaling = "Bandwidth(s)",
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
    nobs = nobs,
    ndim = ndim,
    ncon = sum(xdati$icon),
    nuno = sum(xdati$iuno),
    nord = sum(xdati$iord),
    xdati = xdati,
    ydati = ydati,
    xnames = xnames,
    ynames = ynames,
    sfactor = list(index = sfactor),
    bandwidth = list(index = bandwidth),
    sumNum = list(index = sumNum),
    dati = list(x = xdati, y = ydati, index = idati),
    varnames = list(x = xnames, y = ynames, index = "index"),
    vartitle = list(x = "Explanatory", y = "Dependent", index = "Explanatory"),
    vartitleabb = list(x = "Exp.", y = "Dep.", index = "Exp."),
    rows.omit = rows.omit,
    nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))

  mybw$klist <- list(
    index =
    list(ckertype = ckertype,
         pckertype = mybw$pckertype))

  if(only.optimize.beta)
    mybw$pmethod <- ifelse(only.optimize.beta, paste("Pilot (bandwidth) +",mybw$pmethod, "(beta)"), mybw$pmethod)
  
  if(!bandwidth.compute)
    mybw$pmethod <- "Manual"

  class(mybw) = "sibandwidth"
  if(!any(is.na(mybw$bandwidth)))
    validateBandwidth(mybw)
  mybw
  
}

print.sibandwidth <- function(x, digits=NULL, ...){
  cat("\nSingle Index Model",
      "\nRegression data (",x$nobs,
      " observations, ",x$ndim," variable(s)):\n\n",sep="")

  print(matrix(x$beta,ncol=x$ndim,dimnames=list(paste("Beta",":",sep=""),x$xnames)))
  cat("Bandwidth: ",x$bw)
  cat(genBwSelStr(x))
  cat(genBwKerStrs(x))

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

coef.sibandwidth <- function(object, ...) {
 tc <- object$beta
 names(tc) <- object$xnames
 return(tc)
}
plot.sibandwidth <- function(...) { npplot(...) }
predict.sibandwidth <- function(...) { eval(npindex(...), env = parent.frame()) }

summary.sibandwidth <- function(object, ...){
  cat("\nSingle Index Model",
      "\nRegression data (",object$nobs,
      " observations, ",object$ndim," variable(s)):\n\n",sep="")

  print(matrix(object$beta,ncol=object$ndim,dimnames=list(paste("Beta",":",sep=""),object$xnames)))
  cat("Bandwidth: ",object$bw)
  cat("\nOptimisation Method: ", object$pomethod)

  cat(genOmitStr(object))
  cat(genBwSelStr(object))
  cat(genBwKerStrs(object))
  cat("\n\n")
}


