kbandwidth <-
  function(bw = stop("kbandwidth:argument 'bw' missing"), ...) {
    UseMethod("kbandwidth")
  }

kbandwidth.integer <-
  function(bw, ...) { kbandwidth.numeric(bw = bw, ...) }

kbandwidth.default <- function(bw, ...){
  kbandwidth.numeric(bw = unlist(bw$bandwidth),
                     bwscaling = FALSE,
                     bwtype = bw$type,
                     ckertype = bw$ckertype,
                     ckerorder = bw$ckerorder,
                     ukertype = bw$ukertype,
                     okertype = bw$okertype,
                     nobs = bw$nobs,
                     xdati = if(is.null(bw$zdati)) bw$xdati else bw$zdati,
                     ydati = bw$ydati,
                     xnames = if(is.null(bw$zdati)) bw$xnames else bw$znames,
                     ynames = bw$ynames,
                     ...)
}

kbandwidth.numeric <-
  function(bw,
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"),
           ckerorder = c(2,4,6,8),
           ukertype = c("aitchisonaitken", "liracine"),
           okertype = c("liracine","wangvanryzin", "nliracine"),
           nobs = NA,
           xdati = NULL,
           ydati = NULL,
           xnames = NULL,
           ynames = NULL,
           ...){

    ndim = length(bw)
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

    if(bwscaling != FALSE)
        stop("npksum only uses raw bandwidths, therefore bwscaling = TRUE is not allowed")
    
    ukertype = match.arg(ukertype)
    okertype = match.arg(okertype)

    porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )
    
    mybw = list(
      bw=bw,
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
      xmcv = mcvConstruct(xdati))

    class(mybw) = "kbandwidth"
    mybw
    
  }

as.double.kbandwidth <- function(x, ...){ x$bw }

print.kbandwidth <- function(x, digits=NULL, ...){
  cat("\nData (",x$nobs," observations, ",x$ndim," variable(s)):\n\n",sep="")
  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$xnames)))
  
  cat("\nKernel Sum Bandwidth",
      if (!identical(x$formula,NULL)) paste("\nFormula:",
                                          deparse(x$formula)),
      "\nBandwidth Type:",x$ptype)


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
