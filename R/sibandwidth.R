.np_round_half_to_even <- function(x) {
  intpart <- trunc(x)
  fracpart <- x - intpart

  if (fracpart < 0.5) {
    intpart
  } else if (fracpart > 0.5) {
    intpart + 1
  } else if ((intpart %% 2) != 0) {
    intpart + 1
  } else {
    intpart
  }
}

.np_sibandwidth_manual_nn_validate <- function(h, nobs, where = "sibandwidth") {
  if (!is.finite(h))
    stop(sprintf("%s: nearest-neighbor bandwidth must be finite", where), call. = FALSE)

  upper <- max(1L, as.integer(nobs) - 1L)
  tol <- sqrt(.Machine$double.eps)
  rounded <- .np_round_half_to_even(h)
  if ((h < 1) || (h > upper) || (abs(h - rounded) > tol)) {
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth must be an integer in [1, %d]",
        where,
        upper
      ),
      call. = FALSE
    )
  }

  invisible(as.double(rounded))
}

sibandwidth <-
  function(beta, h,
           method=c("ichimura","kleinspady"),
           regtype = c("lc", "ll", "lp"),
           basis = c("glp", "additive", "tensor"),
           degree = NULL,
           bernstein.basis = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ckerbound = c("none","range","fixed"),
           ckerlb = NULL,
           ckerub = NULL,
           fval = NA,
           ifval = NA,
           num.feval = NA,
           num.feval.fast = NA,
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
           total.time = NA,
           ...){

  ndim = length(beta)
  spec <- npCanonicalConditionalRegSpec(
    regtype = regtype,
    basis = basis,
    degree = degree,
    bernstein.basis = bernstein.basis,
    ncon = 1L,
    where = "sibandwidth"
  )
  regtype <- spec$regtype
  basis <- spec$basis
  degree <- spec$degree
  bernstein.basis <- spec$bernstein.basis
  method = match.arg(method)
  ckertype = match.arg(ckertype)
  ckerbound = match.arg(ckerbound)
  bwtype <- match.arg(bwtype)
  
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

  porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )
  cbounds <- npKernelBoundsResolve(
    dati = xdati,
    varnames = xnames,
    kerbound = ckerbound,
    kerlb = ckerlb,
    kerub = ckerub,
    argprefix = "cker")
  bounded_nonfixed_supported <- bwtype %in% c("generalized_nn", "adaptive_nn")
  if (bwtype != "fixed" && cbounds$bound != "none" && !bounded_nonfixed_supported)
    stop("finite continuous kernel bounds require bwtype = \"fixed\"")
  if (bwtype != "fixed" && (!bandwidth.compute || h != 0))
    .np_sibandwidth_manual_nn_validate(h = h, nobs = nobs, where = "sibandwidth")

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
      ll = "Local-Linear",
      lp = "Local-Polynomial"),
    basis = basis,
    degree = degree,
    bernstein.basis = bernstein.basis,
    regtype.engine = spec$regtype.engine,
    basis.engine = spec$basis.engine,
    degree.engine = spec$degree.engine,
    bernstein.basis.engine = spec$bernstein.basis.engine,
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
    num.feval = num.feval,
    num.feval.fast = num.feval.fast,
    numimp = numimp,
    fval.vector = fval.vector,
    pscaling = npBandwidthSummaryLabel(bwtype = bwtype),
    type = bwtype,
    ptype = switch( bwtype,
      fixed = "Fixed",
      generalized_nn = "Generalized Nearest Neighbour",
      adaptive_nn = "Adaptive Nearest Neighbour" ),
    ckertype = ckertype,    
    ckerorder = ckerorder,
    ckerbound = cbounds$bound,
    ckerlb = cbounds$lb,
    ckerub = cbounds$ub,
    pckertype = cktToPrint(ckertype, order = porder, kerbound = cbounds$bound),
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
    nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
    total.time = total.time)

  mybw$klist <- list(
    index =
    list(ckertype = ckertype,
         ckerbound = cbounds$bound,
         ckerlb = cbounds$lb,
         ckerub = cbounds$ub,
         pckertype = mybw$pckertype))

  if(only.optimize.beta)
    mybw$pmethod <- paste("Pilot (bandwidth) +", mybw$pmethod, "(beta)")
  
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
predict.sibandwidth <- function(...) { do.call(npindex, list(...)) }

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
  cat(genTimingStr(object))
  cat("\n\n")
}
