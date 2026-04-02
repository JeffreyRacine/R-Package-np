plbandwidth <-
  function(bws = stop("plbandwidth: bws missing"),
           regtype = c("lc","ll","lp"),
           basis = c("glp","additive","tensor"),
           degree = NULL,
           bernstein.basis = FALSE,
           bwmethod = c("cv.ls","cv.aic"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ckerbound = c("none","range","fixed"),
           ckerlb = NULL,
           ckerub = NULL,
           ukertype = c("aitchisonaitken", "liracine"),
           okertype = c("liracine","wangvanryzin","racineliyan"),
           xdati, ydati, zdati,
           xnames, ynames, znames,
           nobs = NA,
           fval = NA,
           ifval = NA,
           num.feval = NA,
           num.feval.fast = NA,
           rows.omit = NA,
           bandwidth.compute = TRUE,
           total.time = NA,...){

    npRejectLegacyLpArgs(names(list(...)), where = "plbandwidth")
    spec <- npCanonicalConditionalRegSpec(
      regtype = regtype,
      basis = basis,
      degree = degree,
      bernstein.basis = bernstein.basis,
      ncon = sum(zdati$icon),
      where = "plbandwidth"
    )
    regtype <- spec$regtype
    basis <- spec$basis
    degree <- spec$degree
    bernstein.basis <- spec$bernstein.basis
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
    cbounds <- npKernelBoundsResolve(
      dati = zdati,
      varnames = znames,
      kerbound = ckerbound,
      kerlb = ckerlb,
      kerub = ckerub,
      argprefix = "cker")
    bounded_nonfixed_supported <- bwtype %in% c("generalized_nn", "adaptive_nn")
    if (bwtype != "fixed" && cbounds$bound != "none" && !bounded_nonfixed_supported)
      stop("finite continuous kernel bounds require bwtype = \"fixed\"")
    ncon <- sum(zdati$icon)
    if (identical(spec$regtype.engine, "lp") && ncon > 0L && is.finite(nobs)) {
      lp.dim <- dim_basis(basis = spec$basis.engine,
                          kernel = TRUE,
                          degree = spec$degree.engine,
                          segments = rep.int(1L, ncon))
      if (is.finite(lp.dim) && lp.dim > (nobs - 1.0))
        stop(sprintf("LP basis dimension (%s) exceeds nobs - 1 (%s); reduce degree",
                     format(lp.dim, trim = TRUE, scientific = FALSE),
                     format(nobs - 1.0, trim = TRUE, scientific = FALSE)))
    }

    porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order",
      "Eighth-Order" )
    ## chug chug
    sfactor <- lapply(seq_along(bws), function(i) { unlist(bws[[i]]$sfactor) })
    bandwidth <- lapply(seq_along(bws), function(i) { unlist(bws[[i]]$bandwidth) })
    sumNum <- lapply(seq_along(bws), function(i) { unlist(bws[[i]]$sumNum) })

    names(sfactor) <- names(bandwidth) <- names(sumNum) <- rep("z", length(bws))

    if (length(rows.omit) == 0)
      rows.omit <- NA

    mybw = list(
      bw=bws,
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
      method = bwmethod,
      pmethod = bwmToPrint(bwmethod),
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
      zndim = bws$yzbw$ndim,
      xndim = length(bws)-1,
      ncon = ncon,
      nuno = sum(zdati$iuno),
      nord = sum(zdati$iord),
      icon = zdati$icon,
      iuno = zdati$iuno,
      iord = zdati$iord,
      xdati = xdati,
      ydati = ydati,
      zdati = zdati,
      xnames = xnames,
      ynames = ynames,
      znames = znames,
      fval = fval,
      ifval = ifval,
      num.feval = num.feval,
      num.feval.fast = num.feval.fast,
      sfactor = sfactor,
      bandwidth = bandwidth,
      sumNum = sumNum,
      varnames = list(x = xnames, y = ynames, z = znames),
      vartitle = list(x = "Explanatory", y = "Dependent", z = "Explanatory"),
      vartitleabb = list(x = "Exp.", y = "Dep.", z = "Exp."),
      dati = list(x = xdati, y = ydati, z = zdati),
      rows.omit = rows.omit,
      nobs.omit = if (identical(rows.omit, NA)) 0 else length(rows.omit),
      total.time = total.time)

    mybw$klist = list(z = list(ckertype = ckertype,
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

    class(mybw) = "plbandwidth"
    if(!any(is.na(mybw$bandwidth)))
      validateBandwidth(mybw)
    mybw
    
  }

print.plbandwidth <- function(x, digits=NULL, ...){
  cat("\nPartially Linear Model",
      "\nRegression data (",x$nobs," observations,\n",
      x$xndim," linear parametric variable(s) ",
      x$zndim," nonparametric variable(s):\n",sep="")

  bwmat = matrix(data = 0, nrow = x$xndim+1, ncol = x$bw$yzbw$ndim)
  
  for (i in seq_along(x$bw))
    bwmat[i,] = x$bw[[i]]$bw
  
  ## perhaps add a column for objective function value?
  print(matrix(bwmat[1,], ncol=x$zndim,
               dimnames=list(paste(x$pscaling,":",sep=""),
                 c("y(z)", replicate(x$zndim-1,"")))))
  cat("\n")
  print(matrix(bwmat[2:(1+x$xndim),], ncol=x$zndim,
               dimnames=list(c(paste(x$pscaling,":",sep=""), replicate(x$xndim-1,"")),
                 c("x(z)", replicate(x$zndim-1,"")))))

  cat(genBwSelStr(x))
  cat(genBwKerStrs(x))

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

predict.plbandwidth <- function(...) { do.call(npplreg, list(...)) }

summary.plbandwidth <- function(object, ...){
  cat("\nPartially Linear Model",
      "\nRegression data (",object$nobs," observations,\n",
      object$xndim," linear parametric variable(s) ",
      object$zndim," nonparametric variable(s):\n",sep="")

  cat(genOmitStr(object))
  cat(genBwSelStr(object))

  child.labels <- c(object$ynames, object$xnames)
  scale.blocks <- vapply(
    seq_along(object$bw),
    function(i) {
      paste0("\n\n", child.labels[[i]], " on z:",
             paste(genBwScaleStrs(object$bw[[i]]), collapse = ""))
    },
    character(1L)
  )

  cat(scale.blocks)
  
  cat(genBwKerStrs(object))
  cat(genTimingStr(object))
  cat("\n\n")  
}
