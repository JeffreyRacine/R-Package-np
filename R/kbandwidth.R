# Numerical helpers consume physical smoothing parameters, never search scale
# factors. Conditional objects have two roles; other bandwidth families retain
# one vector (possibly wrapped in a named list for display).
.np_physical_bandwidth <- function(bws, role = NULL) {
  raw <- bws[[if (is.null(role)) "bw" else paste0(role, "bw"), exact = TRUE]]
  if (!isTRUE(bws[["scaling", exact = TRUE]]))
    return(raw)
  retained <- bws[["bandwidth", exact = TRUE]]
  physical <- if (is.null(role)) unlist(retained, use.names = FALSE) else
    retained[[role, exact = TRUE]]
  if (!is.numeric(physical) || length(physical) != length(raw) ||
      any(!is.finite(physical)))
    stop("numerical operator requires retained physical bandwidths", call. = FALSE)
  physical
}

kbandwidth <-
  function(bw = stop("kbandwidth:argument 'bw' missing"), ...) {
    UseMethod("kbandwidth")
  }

kbandwidth.integer <-
  function(bw, ...) { kbandwidth.numeric(bw = bw, ...) }

# Density/distribution families and raw kernel sums use different spellings
# for the same ordered mass. Translate semantics at the object boundary, not
# at the numerical kernel or by rescaling a finished estimate.
.np_density_okertype <- function(okertype) {
  if (identical(okertype, "liracine")) "nliracine" else okertype
}

.np_kbandwidth_okertype <- function(bw) {
  if (inherits(bw, c("bandwidth", "dbandwidth")))
    .np_density_okertype(bw[["okertype", exact = TRUE]])
  else
    bw[["okertype", exact = TRUE]]
}

kbandwidth.default <- function(bw, ...){
  kbandwidth.numeric(bw = unlist(bw$bandwidth),
                     bwscaling = FALSE,
                     bwtype = bw$type,
                     ckertype = bw$ckertype,
                     ckerorder = bw$ckerorder,
                     ckerbound = if (!is.null(bw$ckerbound)) bw$ckerbound else "none",
                     ckerlb = if (!is.null(bw$ckerlb)) bw$ckerlb else NULL,
                     ckerub = if (!is.null(bw$ckerub)) bw$ckerub else NULL,
                     ukertype = bw$ukertype,
                     okertype = .np_kbandwidth_okertype(bw),
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
           ckertype = c("gaussian","epanechnikov","uniform","beta"),
           ckerorder = c(2,4,6,8),
           ckerbound = c("none","range","fixed"),
           ckerlb = NULL,
           ckerub = NULL,
           ukertype = c("aitchisonaitken", "liracine"),
           okertype = c("liracine","wangvanryzin","racineliyan","nliracine"),
           nobs = NA,
           xdati = NULL,
           ydati = NULL,
           xnames = NULL,
           ynames = NULL,
           ...){

    ndim = length(bw)
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

    if(bwscaling != FALSE)
        stop("npksum only uses raw bandwidths, therefore bwscaling = TRUE is not allowed")
    
    ukertype = match.arg(ukertype)
    okertype = match.arg(okertype)
    cbounds <- npKernelBoundsResolve(
      dati = xdati,
      varnames = xnames,
      kerbound = ckerbound,
      kerlb = ckerlb,
      kerub = ckerub,
      argprefix = "cker",
      range.policy = if (identical(ckertype, "beta"))
        "beta_half_spacing" else "exact")
    npValidateBetaKernelSpecification(
      ckertype = ckertype,
      ckerorder = ckerorder,
      bwtype = bwtype,
      ckerbound = cbounds$bound,
      ckerlb = cbounds$lb,
      ckerub = cbounds$ub,
      dati = xdati,
      bw = bw,
      bandwidth.compute = FALSE,
      where = "beta kernel sums",
      allow.categorical = TRUE
    )
    bounded_nonfixed_supported <- bwtype %in% c("generalized_nn", "adaptive_nn")
    if (bwtype != "fixed" && cbounds$bound != "none" && !bounded_nonfixed_supported)
      stop("finite continuous kernel bounds require bwtype = \"fixed\"")

    porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )
    
    mybw = list(
      bw=bw,
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
