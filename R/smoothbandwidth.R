 .np_scbandwidth_manual_nn_validate <- function(bw, nobs, where = "scbandwidth") {
  vapply(bw, function(h) {
    .np_sibandwidth_manual_nn_validate(h = h, nobs = nobs, where = where)
  }, numeric(1))
}

scbandwidth <-
  function(bw = stop("scbandwidth:argument 'bw' missing"),
           regtype = c("lc","ll","lp"),
           basis = c("glp","additive","tensor"),
           degree = NULL,
           bernstein.basis = FALSE,
           bwmethod = c("cv.ls", "manual"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ckerbound = c("none","range","fixed"),
           ckerlb = NULL,
           ckerub = NULL,
           ukertype = c("aitchisonaitken", "liracine"),
           okertype = c("liracine","wangvanryzin","racineliyan"),
           fval = NA,
           ifval = NA,
           num.feval = NA,
           num.feval.fast = NA,
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
  npRejectLegacyLpArgs(names(list(...)), where = "scbandwidth")
  regtype = match.arg(regtype)
  basis <- npValidateLpBasis(regtype = regtype, basis = basis)
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

  ukertype = match.arg(ukertype)
  okertype = match.arg(okertype)

  tdati <- zdati

  if(is.null(zdati))
    tdati <- xdati

  cbounds <- npKernelBoundsResolve(
    dati = tdati,
    varnames = if(is.null(znames)) xnames else znames,
    kerbound = ckerbound,
    kerlb = ckerlb,
    kerub = ckerub,
    argprefix = "cker")
  bounded_nonfixed_supported <- bwtype %in% c("generalized_nn", "adaptive_nn")
  if (bwtype != "fixed" && cbounds$bound != "none" && !bounded_nonfixed_supported)
    stop("finite continuous kernel bounds require bwtype = \"fixed\"")
  if (bwtype != "fixed" && (!bandwidth.compute || any(bw != 0)))
    bw <- .np_scbandwidth_manual_nn_validate(bw = bw, nobs = nobs, where = "scbandwidth")
  ncon <- sum(tdati$icon)
  degree <- npValidateGlpDegree(regtype = regtype,
                                degree = degree,
                                ncon = ncon)
  bernstein.basis <- npValidateGlpBernstein(regtype = regtype,
                                            bernstein.basis = bernstein.basis)
  if (identical(regtype, "lp") && ncon > 0L && is.finite(nobs)) {
    lp.dim <- dim_basis(basis = basis,
                        kernel = TRUE,
                        degree = degree,
                        segments = rep.int(1L, ncon))
    if (is.finite(lp.dim) && lp.dim > (nobs - 1.0))
      stop(sprintf("LP basis dimension (%s) exceeds nobs - 1 (%s); reduce degree",
                   format(lp.dim, trim = TRUE, scientific = FALSE),
                   format(nobs - 1.0, trim = TRUE, scientific = FALSE)))
  }

  porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )

  if (!identical(sfactor,NA)){
    sumNum <- sapply(seq_len(ndim), function(i) {
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
    regtype = regtype,
    pregtype = switch(regtype,
      lc = "Local-Constant",
      ll = "Local-Linear",
      lp = "Local-Polynomial"),
    basis = basis,
    degree = degree,
    bernstein.basis = bernstein.basis,
    method = bwmethod,
    pmethod = bwmToPrint(bwmethod),
    pomethod = switch(optim.method,
      "Nelder-Mead" = "Nelder-Mead",
      "BFGS" = "BFGS",
      "CG" = "CG", "NA"),
    fval = fval,
    ifval = ifval,
    num.feval = num.feval,
    num.feval.fast = num.feval.fast,
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
    pokertype = oktToPrint(okertype),
    nobs = nobs,
    ndim = ndim,
    ncon = ncon,
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
    nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
    total.time = total.time)

  mybw$klist <-
      list(list(ckertype = ckertype,
                pckertype = mybw$pckertype,
                ukertype = ukertype,
                pukertype = mybw$pukertype,
                okertype = okertype,
                pokertype = mybw$pokertype))

  zorx <- if (is.null(zdati)) "x" else "z"
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

predict.scbandwidth <- function(...) { do.call(npscoef, list(...)) }

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
