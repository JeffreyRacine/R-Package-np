npregbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npregbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat"),
                                     eval_env = parent.frame())
    UseMethod("npregbw", target)
  }

npregbw.formula <-
  function(formula, data, subset, na.action, call, ...){

    formula.terms <- terms(formula)

    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = formula.terms,
                        data = environment(formula),
                        eval_env = environment(formula))
    else .np_terms_ts_mask(terms_obj = formula.terms,
                           data = data,
                           eval_env = environment(formula))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    if(all(orig.ts)){
      args <- (as.list(attr(formula.terms, "variables"))[-1])
      formula <- formula.terms
      attr(formula, "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
      mf[["formula"]] <- formula
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(formula.terms, "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      formula <- formula.terms
      attr(formula, "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
      mf[["formula"]] <- formula
    }

    mf[[1]] <- as.name("model.frame")
    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())

    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    dots <- list(...)
    tbw <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat), dots))

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)

    tbw
  }

npregbw.NULL <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npregbw()")
    mc <- match.call(expand.dots = FALSE)
    dots <- list(...)
    dot.names <- names(dots)
    degree.select.value <- if ("degree.select" %in% dot.names) {
      match.arg(as.character(dots$degree.select[[1L]]), c("manual", "coordinate", "exhaustive"))
    } else {
      "manual"
    }
    automatic.degree.search <- isTRUE(dots$nomad) || !identical(degree.select.value, "manual")
    regtype.value <- if ("regtype" %in% dot.names) {
      match.arg(as.character(dots$regtype[[1L]]), c("lc", "ll", "lp"))
    } else if (isTRUE(dots$nomad)) {
      "lp"
    } else {
      "lc"
    }
    search.engine.value <- if ("search.engine" %in% dot.names) {
      match.arg(as.character(dots$search.engine[[1L]]), c("nomad+powell", "cell", "nomad"))
    } else if (isTRUE(dots$nomad)) {
      "nomad+powell"
    } else {
      "nomad+powell"
    }
    .np_nomad_validate_inner_multistart(
      call_names = names(mc),
      dot.args = dots,
      regtype = regtype.value,
      automatic.degree.search = automatic.degree.search,
      search.engine = search.engine.value
    )
    if (.npRmpi_autodispatch_active() && !isTRUE(automatic.degree.search))
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_as_generic_call("npregbw", mc),
        parent.frame()))

    xdat <- toFrame(xdat)

    bws = double(dim(xdat)[2])

    tbw <- npregbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)

    tbw
  }

npregbw.rbandwidth <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws,
           bandwidth.compute = TRUE,
           cfac.dir = 2.5*(3.0-sqrt(5)),
           scale.factor.init = 0.5,
           dfac.dir = 0.25*(3.0-sqrt(5)),
           dfac.init = 0.375,
           dfc.dir = 3,
           ftol = 1.490116e-07,
           scale.factor.init.upper = 2.0,
           hbd.dir = 1,
           hbd.init = 0.9,
           initc.dir = 1.0,
           initd.dir = 1.0,
           invalid.penalty = c("baseline","dbmax"),
           itmax = 10000,
           lbc.dir = 0.5,
           scale.factor.init.lower = 0.1,
           lbd.dir = 0.1,
           lbd.init = 0.1,
           nmulti,
           penalty.multiplier = 10,
           remin = TRUE,
           scale.init.categorical.sample = FALSE,
           scale.factor.search.lower = NULL,
           small = 1.490116e-05,
           tol = 1.490116e-04,
           transform.bounds = FALSE,
           ...){
    elapsed.start <- proc.time()[3]
    xdat <- toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(dim(xdat)[2])
    }
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    remin <- npValidateScalarLogical(remin, "remin")
    scale.init.categorical.sample <-
      npValidateScalarLogical(scale.init.categorical.sample, "scale.init.categorical.sample")
    transform.bounds <- npValidateScalarLogical(transform.bounds, "transform.bounds")
    itmax <- npValidatePositiveInteger(itmax, "itmax")
    ftol <- npValidatePositiveFiniteNumeric(ftol, "ftol")
    tol <- npValidatePositiveFiniteNumeric(tol, "tol")
    small <- npValidatePositiveFiniteNumeric(small, "small")
    penalty.multiplier <- npValidatePositiveFiniteNumeric(penalty.multiplier, "penalty.multiplier")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)
    .npRmpi_require_active_slave_pool(where = "npregbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_as_generic_call("npregbw", match.call()),
        parent.frame()))

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (length(bws$bw) != dim(xdat)[2])
      stop("length of bandwidth vector does not match number of columns of 'xdat'")

    npValidateRegressionNnLowerBound(
      bws,
      where = "npregbw",
      allow.zero.placeholder = TRUE
    )

    if ((any(bws$icon) &&
         !all(vapply(xdat[, bws$icon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iord) &&
         !all(vapply(xdat[, bws$iord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$iuno) &&
         !all(vapply(xdat[, bws$iuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'xdat' in type")

    if (dim(xdat)[1] != length(ydat))
      stop("number of regression data and response data do not match")

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- attr(na.omit(data.frame(xdat,ydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]

    nrow = dim(xdat)[1]
    ncol = dim(xdat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat = toMatrix(xdat)

    runo = xdat[, bws$iuno, drop = FALSE]
    rcon = xdat[, bws$icon, drop = FALSE]
    rord = xdat[, bws$iord, drop = FALSE]

    tbw <- bws
    tbw$basis <- npValidateLpBasis(regtype = tbw$regtype,
                                   basis = tbw$basis)
    tbw$degree <- npValidateGlpDegree(regtype = tbw$regtype,
                                          degree = tbw$degree,
                                          ncon = tbw$ncon)
    tbw$bernstein.basis <- npValidateGlpBernstein(regtype = tbw$regtype,
                                                bernstein.basis = tbw$bernstein.basis)

    mysd <- EssDee(rcon)
    nconfac <- nrow^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- nrow^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

    reg.c <- npRegtypeToC(regtype = tbw$regtype,
                          degree = tbw$degree,
                          ncon = tbw$ncon,
                          context = "npregbw")
    npCheckRegressionDesignCondition(reg.code = reg.c$code,
                                     xcon = rcon,
                                     basis = tbw$basis,
                                     degree = tbw$degree,
                                     bernstein.basis = tbw$bernstein.basis,
                                     where = "npregbw")
    degree.c <- if (tbw$ncon > 0) {
      as.integer(if (is.null(reg.c$degree)) rep.int(0L, tbw$ncon) else reg.c$degree)
    } else {
      integer(1)
    }

    if (bandwidth.compute){
      cont.start <- npContinuousSearchStartControls(
        scale.factor.init.lower,
        scale.factor.init.upper,
        scale.factor.init,
        scale.factor.search.lower,
        where = "npregbw"
      )
      myopti = list(num_obs_train = dim(xdat)[1],
        iMultistart = IMULTI_TRUE,
        iNum_Multistart = nmulti,
        int_use_starting_values = (if (all(bws$bw==0)) USE_START_NO else USE_START_YES),
        int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
        BANDWIDTH_reg_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=(if (remin) RE_MIN_TRUE else RE_MIN_FALSE),
        int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
        bwmethod = switch(bws$method,
          cv.aic = BWM_CVAIC,
          cv.ls = BWM_CVLS),
        kerneval = switch(bws$ckertype,
          gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
          uniform = CKER_UNI,
          "truncated gaussian" = CKER_TGAUSS),
        ukerneval = switch(bws$ukertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        okerneval = switch(bws$okertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_LR,
        "racineliyan" = OKER_RLY),
        nuno = bws$nuno,
        nord = bws$nord,
        ncon = bws$ncon,
        regtype = reg.c$code,
        int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
        scale.init.categorical.sample = scale.init.categorical.sample,
        dfc.dir = dfc.dir,
        transform.bounds = transform.bounds)

      myoptd = list(ftol=ftol, tol=tol, small=small,
        lbc.dir = lbc.dir, cfac.dir = cfac.dir, initc.dir = initc.dir,
        lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir,
        lbc.init = cont.start$scale.factor.init.lower,
        hbc.init = cont.start$scale.factor.init.upper,
        cfac.init = cont.start$scale.factor.init,
        lbd.init = lbd.init, hbd.init = hbd.init, dfac.init = dfac.init,
        nconfac = nconfac, ncatfac = ncatfac,
        scale.factor.lower.bound = scale.factor.search.lower)

      cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

      system.time(myout <-
        .Call("C_np_regression_bw",
              as.double(runo), as.double(rord), as.double(rcon), as.double(ydat),
              as.double(mysd),
              as.integer(myopti), as.double(myoptd),
              as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
              as.integer(nmulti),
              as.integer(penalty_mode),
              as.double(penalty.multiplier),
              as.integer(degree.c),
              as.integer(isTRUE(tbw$bernstein.basis)),
              as.integer(npLpBasisCode(tbw$basis)),
              as.double(cker.bounds.c$lb),
              as.double(cker.bounds.c$ub),
              PACKAGE = "npRmpi"))[1]


      rorder = numeric(ncol)
      ord_idx <- seq_len(ncol)
      rorder[c(ord_idx[bws$icon], ord_idx[bws$iuno], ord_idx[bws$iord])] <- ord_idx

      tbw$bw <- myout$bw[rorder]
      tbw$fval <- myout$fval[1]
      tbw$ifval <- myout$fval[2]
      tbw$initial.fval <- if (length(myout$fval) >= 3L) myout$fval[3L] else NA_real_
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
      tbw$fval.history <- myout$fval.history
      tbw$eval.history <- myout$eval.history
      tbw$invalid.history <- myout$invalid.history
      tbw$timing <- myout$timing
    }

    tbw$total.time <- proc.time()[3] - elapsed.start
    initial.fval <- tbw$initial.fval

    tbw$sfactor <- tbw$bandwidth <- tbw$bw

    if (tbw$nuno > 0){
      if(tbw$scaling){
        tbw$bandwidth[tbw$xdati$iuno] <- tbw$bandwidth[tbw$xdati$iuno]*ncatfac
      } else {
        tbw$sfactor[tbw$xdati$iuno] <- tbw$sfactor[tbw$xdati$iuno]/ncatfac
      }
    }

    if (tbw$nord > 0){
      if(tbw$scaling){
        tbw$bandwidth[tbw$xdati$iord] <- tbw$bandwidth[tbw$xdati$iord]*ncatfac
      } else {
        tbw$sfactor[tbw$xdati$iord] <- tbw$sfactor[tbw$xdati$iord]/ncatfac
      }
    }

    if (tbw$ncon > 0){
      dfactor <- mysd*nconfac

      if (tbw$scaling) {
        tbw$bandwidth[tbw$xdati$icon] <- tbw$bandwidth[tbw$xdati$icon]*dfactor
      } else {
        tbw$sfactor[tbw$xdati$icon] <- tbw$sfactor[tbw$xdati$icon]/dfactor
      }
    }



    tbw <- rbandwidth(bw = tbw$bw,
                      regtype = tbw$regtype,
                      basis = tbw$basis,
                      degree = tbw$degree,
                      bernstein.basis = tbw$bernstein.basis,
                      bwmethod = tbw$method,
                      bwscaling = tbw$scaling,
                      bwtype = tbw$type,
                      ckertype = tbw$ckertype,
                      ckerorder = tbw$ckerorder,
                      ckerbound = tbw$ckerbound,
                      ckerlb = tbw$ckerlb,
                      ckerub = tbw$ckerub,
                      ukertype = tbw$ukertype,
                      okertype = tbw$okertype,
                      fval = tbw$fval,
                      ifval = tbw$ifval,
                      num.feval = tbw$num.feval,
                      num.feval.fast = tbw$num.feval.fast,
                      fval.history = tbw$fval.history,
                      eval.history = tbw$eval.history,
                      invalid.history = tbw$invalid.history,
                      nobs = tbw$nobs,
                      xdati = tbw$xdati,
                      ydati = tbw$ydati,
                      xnames = tbw$xnames,
                      ynames = tbw$ynames,
                      sfactor = tbw$sfactor,
                      bandwidth = tbw$bandwidth,
                      rows.omit = rows.omit,
                      nconfac = nconfac,
                      ncatfac = ncatfac,
                      sdev = mysd,
                      bandwidth.compute = bandwidth.compute,
                      timing = tbw$timing,
                      total.time = tbw$total.time)
    tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)
    tbw$initial.fval <- if (!is.null(initial.fval)) initial.fval else NA_real_
    tbw
  }

.npregbw_build_rbandwidth <- function(xdat,
                                      ydat,
                                      bws,
                                      bandwidth.compute,
                                      reg.args,
                                      yname) {
  bw.args <- c(
    list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      xnames = names(xdat),
      ynames = yname,
      bandwidth.compute = bandwidth.compute
    ),
    reg.args
  )
  out <- do.call(rbandwidth, bw.args)
  if (!is.null(reg.args$scale.factor.search.lower)) {
    out$scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      reg.args$scale.factor.search.lower
    )
  }
  out
}

.npregbw_call_fixed_degree_core <- function(xdat,
                                            ydat,
                                            bws,
                                            nmulti = 1L,
                                            itmax = 0L,
                                            remin = FALSE,
                                            scale.init.categorical.sample = FALSE,
                                            ftol = 0,
                                            tol = 0,
                                            small = 0,
                                            lbc.dir = 0,
                                            dfc.dir = 0,
                                            cfac.dir = 0,
                                            initc.dir = 0,
                                            lbd.dir = 0,
                                            hbd.dir = 0,
                                            dfac.dir = 0,
                                            initd.dir = 0,
                                            scale.factor.init.lower = 0,
                                            scale.factor.init.upper = 0,
                                            scale.factor.init = 0,
                                            lbd.init = 0,
                                            hbd.init = 0,
                                            dfac.init = 0,
                                            invalid.penalty = c("baseline", "dbmax"),
                                            penalty.multiplier = 10,
                                            transform.bounds = FALSE,
                                            scale.factor.search.lower = NULL,
                                            eval.only = FALSE) {
  invalid.penalty <- match.arg(invalid.penalty)
  scale.factor.search.lower <- npResolveScaleFactorLowerBound(
    if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
  )
  if (!isTRUE(eval.only)) {
    cont.start <- npContinuousSearchStartControls(
      scale.factor.init.lower,
      scale.factor.init.upper,
      scale.factor.init,
      scale.factor.search.lower,
      where = "npregbw"
    )
    scale.factor.init.lower <- cont.start$scale.factor.init.lower
    scale.factor.init.upper <- cont.start$scale.factor.init.upper
    scale.factor.init <- cont.start$scale.factor.init
  }

  xdat <- toFrame(xdat)
  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector")

  if (length(bws$bw) != dim(xdat)[2])
    stop("length of bandwidth vector does not match number of columns of 'xdat'")

  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE

  xdat <- xdat[keep.rows,, drop = FALSE]
  ydat <- ydat[keep.rows]

  if (is.factor(ydat))
    ydat <- dlev(ydat)[as.integer(ydat)]
  else
    ydat <- as.double(ydat)

  xmat <- toMatrix(xdat)
  runo <- xmat[, bws$iuno, drop = FALSE]
  rcon <- xmat[, bws$icon, drop = FALSE]
  rord <- xmat[, bws$iord, drop = FALSE]

  mysd <- EssDee(rcon)
  nrow <- dim(xmat)[1L]
  nconfac <- nrow^(-1.0 / (2.0 * bws$ckerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$ckerorder + bws$ncon))

  penalty_mode <- if (invalid.penalty == "baseline") 1L else 0L
  reg.c <- npRegtypeToC(regtype = bws$regtype,
                        degree = bws$degree,
                        ncon = bws$ncon,
                        context = "npregbw")
  degree.c <- if (bws$ncon > 0) as.integer(bws$degree) else integer(1L)
  nmulti <- as.integer(nmulti[1L])

  myopti <- list(
    num_obs_train = nrow,
    iMultistart = IMULTI_TRUE,
    iNum_Multistart = nmulti,
    int_use_starting_values = if (all(bws$bw == 0)) USE_START_NO else USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_reg_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = itmax,
    int_RESTART_FROM_MIN = if (isTRUE(remin)) RE_MIN_TRUE else RE_MIN_FALSE,
    int_MINIMIZE_IO = if (isTRUE(eval.only) || !isTRUE(getOption("np.messages"))) IO_MIN_TRUE else IO_MIN_FALSE,
    bwmethod = switch(bws$method,
      cv.aic = BWM_CVAIC,
      cv.ls = BWM_CVLS),
    kerneval = switch(bws$ckertype,
      gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS),
    ukerneval = switch(bws$ukertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    okerneval = switch(bws$okertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR,
      racineliyan = OKER_RLY),
    nuno = bws$nuno,
    nord = bws$nord,
    ncon = bws$ncon,
    regtype = reg.c$code,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    scale.init.categorical.sample = scale.init.categorical.sample,
    dfc.dir = dfc.dir,
    transform.bounds = transform.bounds
  )

  myoptd <- list(
    ftol = ftol,
    tol = tol,
    small = small,
    lbc.dir = lbc.dir,
    cfac.dir = cfac.dir,
    initc.dir = initc.dir,
    lbd.dir = lbd.dir,
    hbd.dir = hbd.dir,
    dfac.dir = dfac.dir,
    initd.dir = initd.dir,
    lbc.init = scale.factor.init.lower,
    hbc.init = scale.factor.init.upper,
    cfac.init = scale.factor.init,
    lbd.init = lbd.init,
    hbd.init = hbd.init,
    dfac.init = dfac.init,
    nconfac = nconfac,
    ncatfac = ncatfac,
    scale.factor.lower.bound = scale.factor.search.lower
  )

  cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

  if (isTRUE(eval.only)) {
    out <- .npRmpi_with_local_regression(.Call(
      C_np_regression_bw_eval,
      as.double(runo),
      as.double(rord),
      as.double(rcon),
      as.double(ydat),
      as.double(mysd),
      as.integer(myopti),
      as.double(myoptd),
      as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
      as.integer(1L),
      as.integer(penalty_mode),
      as.double(penalty.multiplier),
      as.integer(degree.c),
      as.integer(isTRUE(bws$bernstein.basis)),
      as.integer(npLpBasisCode(bws$basis)),
      as.double(cker.bounds.c$lb),
      as.double(cker.bounds.c$ub),
      PACKAGE = "npRmpi"
    ))
  } else {
    out <- .npRmpi_with_local_regression(.Call(
      C_np_regression_bw,
      as.double(runo),
      as.double(rord),
      as.double(rcon),
      as.double(ydat),
      as.double(mysd),
      as.integer(myopti),
      as.double(myoptd),
      as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
      as.integer(max(1L, nmulti)),
      as.integer(penalty_mode),
      as.double(penalty.multiplier),
      as.integer(degree.c),
      as.integer(isTRUE(bws$bernstein.basis)),
      as.integer(npLpBasisCode(bws$basis)),
      as.double(cker.bounds.c$lb),
      as.double(cker.bounds.c$ub),
      PACKAGE = "npRmpi"
    ))
  }

  rorder <- numeric(length(bws$bw))
  ord.idx <- seq_along(bws$bw)
  rorder[c(ord.idx[bws$icon], ord.idx[bws$iuno], ord.idx[bws$iord])] <- ord.idx

  list(
    objective = as.numeric(out$fval[1L]),
    ifval = as.numeric(out$fval[2L]),
    bw = as.numeric(out$bw[rorder]),
    num.feval = sum(out$eval.history[is.finite(out$eval.history)]),
    num.feval.fast = as.numeric(out$fast.history[1L]),
    fval.history = out$fval.history,
    eval.history = out$eval.history,
    invalid.history = out$invalid.history,
    timing = out$timing
  )
}

.npregbw_finalize_fixed_degree_payload <- function(xdat,
                                                   ydat,
                                                   bws,
                                                   core,
                                                   total.time) {
  tbw <- bws
  tbw$bw <- core$bw
  tbw$fval <- core$objective
  tbw$ifval <- core$ifval
  tbw$num.feval <- core$num.feval
  tbw$num.feval.fast <- core$num.feval.fast
  tbw$fval.history <- core$fval.history
  tbw$eval.history <- core$eval.history
  tbw$invalid.history <- core$invalid.history
  tbw$timing <- core$timing

  payload <- npregbw.rbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = tbw,
    bandwidth.compute = FALSE
  )
  if (!is.null(payload$method) && length(payload$method))
    payload$pmethod <- bwmToPrint(as.character(payload$method[1L]))
  payload$timing <- core$timing
  payload$total.time <- total.time
  payload
}

.npregbw_run_fixed_degree_source_of_truth <- function(xdat,
                                                      ydat,
                                                      bws,
                                                      opt.args) {
  opt.value <- function(name, default) {
    if (is.null(opt.args[[name]])) default else opt.args[[name]]
  }

  elapsed.start <- proc.time()[3L]
  core <- .npregbw_call_fixed_degree_core(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    nmulti = opt.value("nmulti", npDefaultNmulti(dim(toFrame(xdat))[2L])),
    itmax = opt.value("itmax", 10000L),
    remin = opt.value("remin", TRUE),
    scale.init.categorical.sample = opt.value("scale.init.categorical.sample", FALSE),
    ftol = opt.value("ftol", 1.490116e-07),
    tol = opt.value("tol", 1.490116e-04),
    small = opt.value("small", 1.490116e-05),
    lbc.dir = opt.value("lbc.dir", 0.5),
    dfc.dir = opt.value("dfc.dir", 3),
    cfac.dir = opt.value("cfac.dir", 2.5 * (3.0 - sqrt(5))),
    initc.dir = opt.value("initc.dir", 1.0),
    lbd.dir = opt.value("lbd.dir", 0.1),
    hbd.dir = opt.value("hbd.dir", 1),
    dfac.dir = opt.value("dfac.dir", 0.25 * (3.0 - sqrt(5))),
    initd.dir = opt.value("initd.dir", 1.0),
    scale.factor.init.lower = opt.value("scale.factor.init.lower", 0.1),
    scale.factor.init.upper = opt.value("scale.factor.init.upper", 2.0),
    scale.factor.init = opt.value("scale.factor.init", 0.5),
    lbd.init = opt.value("lbd.init", 0.1),
    hbd.init = opt.value("hbd.init", 0.9),
    dfac.init = opt.value("dfac.init", 0.375),
    invalid.penalty = opt.value("invalid.penalty", "baseline"),
    penalty.multiplier = opt.value("penalty.multiplier", 10),
    transform.bounds = opt.value("transform.bounds", FALSE),
    scale.factor.search.lower = opt.value("scale.factor.search.lower", NULL),
    eval.only = FALSE
  )

  .npregbw_finalize_fixed_degree_payload(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    core = core,
    total.time = proc.time()[3L] - elapsed.start
  )
}

.npregbw_run_fixed_degree_source_of_truth_bcast_payload <- function(xdat, ydat, bws, opt.args) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.messages <- getOption("np.messages")
  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)

  options(npRmpi.autodispatch.disable = TRUE)
  if (!isTRUE(rank == 0L))
    options(np.messages = FALSE)

  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(np.messages = old.messages), add = TRUE)

  .npregbw_run_fixed_degree_source_of_truth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    opt.args = opt.args
  )
}

.npregbw_run_fixed_degree_source_of_truth_collective <- function(xdat,
                                                                 ydat,
                                                                 bws,
                                                                 opt.args,
                                                                 comm = 1L) {
  if (.npRmpi_has_active_slave_pool(comm = comm) &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    mc <- substitute(
      get(".npregbw_run_fixed_degree_source_of_truth_bcast_payload", envir = asNamespace("npRmpi"), inherits = FALSE)(
        XDAT,
        YDAT,
        BWS,
        OPTARGS
      ),
      list(
        XDAT = xdat,
        YDAT = ydat,
        BWS = bws,
        OPTARGS = opt.args
      )
    )

    return(.npRmpi_bcast_cmd_expr(mc, comm = comm, caller.execute = TRUE))
  }

  .npregbw_run_fixed_degree_source_of_truth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    opt.args = opt.args
  )
}

.npregbw_run_fixed_degree <- function(xdat, ydat, bws, reg.args, opt.args, yname) {
  tbw <- .npregbw_build_rbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args,
    yname = yname
  )

  do.call(npregbw.rbandwidth, c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args))
}

.npregbw_run_fixed_degree_bcast_payload <- function(xdat, ydat, bws, reg.args, opt.args, yname) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.messages <- getOption("np.messages")
  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)

  options(npRmpi.autodispatch.disable = TRUE)
  if (!isTRUE(rank == 0L))
    options(np.messages = FALSE)

  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(np.messages = old.messages), add = TRUE)

  .npregbw_run_fixed_degree(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    reg.args = reg.args,
    opt.args = opt.args,
    yname = yname
  )
}

.npregbw_run_fixed_degree_collective <- function(xdat,
                                                 ydat,
                                                 bws,
                                                 reg.args,
                                                 opt.args,
                                                 yname,
                                                 comm = 1L) {
  if (.npRmpi_has_active_slave_pool(comm = comm) &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    mc <- substitute(
      get(".npregbw_run_fixed_degree_bcast_payload", envir = asNamespace("npRmpi"), inherits = FALSE)(
        XDAT,
        YDAT,
        BWS,
        REGARGS,
        OPTARGS,
        YNAME
      ),
      list(
        XDAT = xdat,
        YDAT = ydat,
        BWS = bws,
        REGARGS = reg.args,
        OPTARGS = opt.args,
        YNAME = yname
      )
    )

    return(.npRmpi_bcast_cmd_expr(mc, comm = comm, caller.execute = TRUE))
  }

  .npregbw_run_fixed_degree(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    reg.args = reg.args,
    opt.args = opt.args,
    yname = yname
  )
}

.npregbw_eval_only <- function(xdat,
                               ydat,
                               bws,
                               invalid.penalty = c("baseline", "dbmax"),
                               penalty.multiplier = 10) {
  out <- .npregbw_call_fixed_degree_core(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    invalid.penalty = invalid.penalty,
    penalty.multiplier = penalty.multiplier,
    eval.only = TRUE
  )

  list(
    objective = out$objective,
    num.feval = out$num.feval,
    num.feval.fast = out$num.feval.fast
  )
}

npRmpiNomadShadowPrepareRegression <- function(runo,
                                               rord,
                                               rcon,
                                               yvec,
                                               mysd,
                                               myopti,
                                               myoptd,
                                               rbw,
                                               penalty.mode,
                                               penalty.multiplier,
                                               degree,
                                               bernstein,
                                               basis,
                                               ckerlb,
                                               ckerub) {
  if (length(myoptd) <= 18L) {
    rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
    if (isTRUE(rank == 0L))
      stop("resident npreg NOMAD shadow options are missing the scale-factor search lower bound", call. = FALSE)
    return(FALSE)
  }

  ok <- .Call(
    "C_np_regression_nomad_shadow_prepare",
    runo,
    rord,
    rcon,
    yvec,
    mysd,
    myopti,
    myoptd,
    rbw,
    penalty.mode,
    penalty.multiplier,
    degree,
    bernstein,
    basis,
    ckerlb,
    ckerub,
    PACKAGE = "npRmpi"
  )

  if (isTRUE(ok))
    return(TRUE)

  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
  if (isTRUE(rank == 0L))
    stop("failed to prepare resident npreg NOMAD shadow state", call. = FALSE)

  FALSE
}

npRmpiNomadShadowEvalRegression <- function(bw, degree) {
  .Call(
    "C_np_regression_nomad_shadow_eval",
    bw,
    degree,
    PACKAGE = "npRmpi"
  )
}

npRmpiNomadShadowClearRegression <- function() {
  .Call("C_np_regression_nomad_shadow_clear", PACKAGE = "npRmpi")
}

npRmpiNomadEvalOnlyRegression <- function(runo,
                                          rord,
                                          rcon,
                                          yvec,
                                          mysd,
                                          myopti,
                                          myoptd,
                                          bw,
                                          penalty.mode,
                                          penalty.multiplier,
                                          degree,
                                          bernstein,
                                          basis,
                                          ckerlb,
                                          ckerub) {
  .Call(
    "C_np_regression_bw_eval",
    runo,
    rord,
    rcon,
    yvec,
    mysd,
    myopti,
    myoptd,
    bw,
    as.integer(1L),
    penalty.mode,
    penalty.multiplier,
    degree,
    bernstein,
    basis,
    ckerlb,
    ckerub,
    PACKAGE = "npRmpi"
  )
}

.npregbw_nomad_shadow_prepare_args <- function(xdat,
                                               ydat,
                                               bws,
                                               start.bw = NULL,
                                               invalid.penalty = c("baseline", "dbmax"),
                                               penalty.multiplier = 10) {
  invalid.penalty <- match.arg(invalid.penalty)

  xdat <- toFrame(xdat)
  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector")

  if (length(bws$bw) != dim(xdat)[2])
    stop("length of bandwidth vector does not match number of columns of 'xdat'")

  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE

  xdat <- xdat[keep.rows,, drop = FALSE]
  ydat <- ydat[keep.rows]

  if (is.factor(ydat))
    ydat <- dlev(ydat)[as.integer(ydat)]
  else
    ydat <- as.double(ydat)

  xmat <- toMatrix(xdat)
  runo <- xmat[, bws$iuno, drop = FALSE]
  rcon <- xmat[, bws$icon, drop = FALSE]
  rord <- xmat[, bws$iord, drop = FALSE]

  mysd <- EssDee(rcon)
  nrow <- dim(xmat)[1L]
  nconfac <- nrow^(-1.0 / (2.0 * bws$ckerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$ckerorder + bws$ncon))
  scale.factor.search.lower <- npResolveScaleFactorLowerBound(
    npGetScaleFactorSearchLower(bws),
    argname = "bws$scale.factor.search.lower"
  )

  penalty_mode <- if (match.arg(invalid.penalty) == "baseline") 1L else 0L
  reg.c <- npRegtypeToC(regtype = bws$regtype,
                        degree = bws$degree,
                        ncon = bws$ncon,
                        context = "npregbw")
  degree.c <- if (bws$ncon > 0) as.integer(bws$degree) else integer(0L)

  myopti <- list(
    num_obs_train = nrow,
    iMultistart = IMULTI_FALSE,
    iNum_Multistart = 0L,
    int_use_starting_values = USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_reg_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = 0L,
    int_RESTART_FROM_MIN = RE_MIN_FALSE,
    int_MINIMIZE_IO = IO_MIN_TRUE,
    bwmethod = switch(bws$method,
      cv.aic = BWM_CVAIC,
      cv.ls = BWM_CVLS),
    kerneval = switch(bws$ckertype,
      gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS),
    ukerneval = switch(bws$ukertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    okerneval = switch(bws$okertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR,
      racineliyan = OKER_RLY),
    nuno = bws$nuno,
    nord = bws$nord,
    ncon = bws$ncon,
    regtype = reg.c$code,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    scale.init.categorical.sample = FALSE,
    dfc.dir = 0,
    transform.bounds = FALSE
  )

  myoptd <- list(
    ftol = 0,
    tol = 0,
    small = 0,
    lbc.dir = 0,
    cfac.dir = 0,
    initc.dir = 0,
    lbd.dir = 0,
    hbd.dir = 0,
    dfac.dir = 0,
    initd.dir = 0,
    lbc.init = 0,
    hbc.init = 0,
    cfac.init = 0,
    lbd.init = 0,
    hbd.init = 0,
    dfac.init = 0,
    nconfac = nconfac,
    ncatfac = ncatfac,
    scale.factor.lower.bound = scale.factor.search.lower
  )

  cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])
  if (is.null(start.bw)) {
    start.bw <- c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])
  } else {
    start.bw <- c(start.bw[bws$icon], start.bw[bws$iuno], start.bw[bws$iord])
  }

  list(
    runo = as.double(runo),
    rord = as.double(rord),
    rcon = as.double(rcon),
    ydat = as.double(ydat),
    mysd = as.double(mysd),
    myopti = as.integer(myopti),
    myoptd = as.double(myoptd),
    rbw = as.double(start.bw),
    penalty_mode = as.integer(penalty_mode),
    penalty_multiplier = as.double(penalty.multiplier),
    degree = as.integer(degree.c),
    bernstein = as.integer(isTRUE(bws$bernstein.basis)),
    basis = as.integer(npLpBasisCode(bws$basis)),
    ckerlb = as.double(cker.bounds.c$lb),
    ckerub = as.double(cker.bounds.c$ub)
  )
}

.npregbw_nomad_shadow_begin <- function(xdat,
                                        ydat,
                                        bws,
                                        start.bw = NULL,
                                        invalid.penalty = c("baseline", "dbmax"),
                                        penalty.multiplier = 10,
                                        comm = 1L,
                                        broadcast = TRUE) {
  prep <- .npregbw_nomad_shadow_prepare_args(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    start.bw = start.bw,
    invalid.penalty = invalid.penalty,
    penalty.multiplier = penalty.multiplier
  )

  mc <- substitute(
    get("npRmpiNomadShadowPrepareRegression", envir = asNamespace("npRmpi"), inherits = FALSE)(
      RUNO,
      RORD,
      RCON,
      YVEC,
      MYSD,
      MYOPTI,
      MYOPTD,
      RBW,
      PMODE,
      PMULT,
      DEGREE,
      BERN,
      BASIS,
      CKERLB,
      CKERUB
    ),
    list(
      RUNO = prep$runo,
      RORD = prep$rord,
      RCON = prep$rcon,
      YVEC = prep$ydat,
      MYSD = prep$mysd,
      MYOPTI = prep$myopti,
      MYOPTD = prep$myoptd,
      RBW = prep$rbw,
      PMODE = prep$penalty_mode,
      PMULT = prep$penalty_multiplier,
      DEGREE = prep$degree,
      BERN = prep$bernstein,
      BASIS = prep$basis,
      CKERLB = prep$ckerlb,
      CKERUB = prep$ckerub
    )
  )

  if (isTRUE(broadcast)) {
    .npRmpi_bcast_cmd_expr(mc, comm = comm, caller.execute = TRUE)
  } else {
    eval(mc, envir = parent.frame())
  }

  list(
    comm = comm,
    icon = bws$icon,
    iuno = bws$iuno,
    iord = bws$iord,
    active = TRUE
  )
}

.npregbw_nomad_shadow_eval <- function(shadow, bw, degree, broadcast = TRUE) {
  flat.bw <- c(bw[shadow$icon], bw[shadow$iuno], bw[shadow$iord])
  mc <- substitute(
    get("npRmpiNomadShadowEvalRegression", envir = asNamespace("npRmpi"), inherits = FALSE)(
      BW,
      DEGREE
    ),
    list(
      BW = as.double(flat.bw),
      DEGREE = as.integer(degree)
    )
  )

  if (isTRUE(broadcast)) {
    as.numeric(.npRmpi_bcast_cmd_expr(mc, comm = shadow$comm, caller.execute = TRUE))
  } else {
    as.numeric(eval(mc, envir = parent.frame()))
  }
}

.npregbw_nomad_shadow_end <- function(shadow, broadcast = TRUE) {
  if (is.null(shadow) || !isTRUE(shadow$active))
    return(invisible(NULL))

  mc <- quote(get("npRmpiNomadShadowClearRegression", envir = asNamespace("npRmpi"), inherits = FALSE)())
  if (isTRUE(broadcast)) {
    .npRmpi_bcast_cmd_expr(mc, comm = shadow$comm, caller.execute = TRUE)
  } else {
    eval(mc, envir = parent.frame())
  }
  invisible(NULL)
}

.npregbw_nomad_controls <- function(search.engine) {
  .np_degree_search_engine_controls(search.engine)
}

.npregbw_nomad_bw_setup <- function(xdat, template, bandwidth.scale.categorical = 1e4) {
  xdat <- toFrame(xdat)
  xmat <- toMatrix(xdat)
  rcon <- xmat[, template$icon, drop = FALSE]
  mysd <- EssDee(rcon)
  nrow <- dim(xmat)[1L]
  nconfac <- nrow^(-1.0 / (2.0 * template$ckerorder + template$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * template$ckerorder + template$ncon))

  cont_scale <- mysd * nconfac
  cont_idx <- which(template$icon)
  cat_idx <- c(which(template$iuno), which(template$iord))

  cat_upper <- numeric(length(cat_idx))
  for (k in seq_along(cat_idx)) {
    i <- cat_idx[k]
    if (template$iuno[i] && identical(template$ukertype, "aitchisonaitken")) {
      nlev <- length(unique(xdat[[i]]))
      cat_upper[k] <- (nlev - 1) / nlev
    } else {
      cat_upper[k] <- 1
    }
  }

  list(
    cont_scale = cont_scale,
    cont_idx = cont_idx,
    cat_idx = cat_idx,
    ncatfac = ncatfac,
    bandwidth.scale.categorical = bandwidth.scale.categorical,
    cat_upper = cat_upper
  )
}

.npregbw_nomad_point_to_bw <- function(point, template, setup) {
  point <- as.numeric(point)
  bws <- numeric(length(template$bw))
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)

  if (ncon > 0L) {
    gamma <- point[seq_len(ncon)]
    ext_bw <- gamma * setup$cont_scale
    bws[setup$cont_idx] <- if (isTRUE(template$scaling)) gamma else ext_bw
  }

  if (ncat > 0L) {
    lambda_scaled <- point[ncon + seq_len(ncat)]
    ext_bw <- lambda_scaled / setup$bandwidth.scale.categorical
    bws[setup$cat_idx] <- if (isTRUE(template$scaling)) ext_bw / setup$ncatfac else ext_bw
  }

  bws
}

.npregbw_nomad_bw_to_point <- function(bws, template, setup) {
  point <- numeric(length(setup$cont_idx) + length(setup$cat_idx))

  if (length(setup$cont_idx) > 0L) {
    raw <- bws[setup$cont_idx]
    point[seq_along(setup$cont_idx)] <- if (isTRUE(template$scaling)) {
      raw
    } else {
      raw / setup$cont_scale
    }
  }

  if (length(setup$cat_idx) > 0L) {
    raw <- bws[setup$cat_idx]
    ext_bw <- if (isTRUE(template$scaling)) raw * setup$ncatfac else raw
    point[length(setup$cont_idx) + seq_along(setup$cat_idx)] <- ext_bw * setup$bandwidth.scale.categorical
  }

  point
}

.npregbw_powell_progress_fields <- function(state,
                                            done = NULL,
                                            detail = NULL,
                                            now = .np_progress_now()) {
  fields <- character()
  elapsed <- max(0, now - state$started)

  fields <- c(fields, sprintf("elapsed %ss", .np_progress_fmt_num(elapsed)))

  if (!is.null(state$nomad_current_degree)) {
    fields <- c(
      fields,
      sprintf("degree %s", .np_degree_format_degree(state$nomad_current_degree))
    )
  }

  if (!is.null(done)) {
    done <- suppressWarnings(as.integer(done)[1L])
    if (!is.na(done) && done >= 1L) {
      fields <- c(fields, sprintf("iter %s", format(done)))
    }
  }

  fields
}

.npregbw_with_powell_refinement_progress <- function(degree, expr) {
  old.state <- .np_progress_runtime$bandwidth_state
  worker_silent <- FALSE
  active.state <- NULL

  if (isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) {
    attach.size <- tryCatch(as.integer(mpi.comm.size(1L)), error = function(e) NA_integer_)
    attach.rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) NA_integer_)

    worker_silent <- !is.na(attach.size) && attach.size > 1L &&
      !is.na(attach.rank) && attach.rank != 0L
  }

  on.exit({
    current.state <- .np_progress_runtime$bandwidth_state
    if (!is.null(active.state) && !is.null(current.state)) {
      .np_progress_end(current.state)
    }
    .np_progress_runtime$bandwidth_state <- old.state
  }, add = TRUE)

  if (!worker_silent) {
    active.state <- .np_progress_begin(
      label = .np_nomad_powell_progress_label(),
      domain = "general",
      surface = "bandwidth"
    )
    active.state$unknown_total_fields <- .npregbw_powell_progress_fields
    active.state$nomad_current_degree <- as.integer(degree)
    active.state <- .np_progress_show_now(active.state)
    .np_progress_runtime$bandwidth_state <- active.state
  } else {
    .np_progress_runtime$bandwidth_state <- NULL
  }

  value <- force(expr)

  if (!is.null(active.state) && is.list(value) && !is.null(value$num.feval)) {
    done <- suppressWarnings(as.integer(value$num.feval[1L]))
    if (!is.na(done) && done >= 1L && !is.null(.np_progress_runtime$bandwidth_state)) {
      .np_progress_runtime$bandwidth_state <- .np_progress_step_at(
        state = .np_progress_runtime$bandwidth_state,
        now = .np_progress_now(),
        done = done,
        force = TRUE
      )
    }
  }

  value
}

npRmpiNomadShadowSearchRegression <- function(xdat,
                                              ydat,
                                              template,
                                              setup,
                                              prep,
                                              reg.args,
                                              degree.search,
                                              x0,
                                              bbin,
                                              lb,
                                              ub,
                                              nomad.nmulti = 1L,
                                              nomad.inner.nmulti = 0L,
                                              penalty.multiplier = 10,
                                              random.seed = 42L) {
  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
  old.messages <- getOption("np.messages")
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)

  if (!isTRUE(rank == 0L))
    options(np.messages = FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  on.exit(options(np.messages = old.messages), add = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)

  set.seed(as.integer(random.seed))
  if (isTRUE(rank == 0L))
    .np_nomad_baseline_note(degree.search$start.degree)

  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)
  prepared <- npRmpiNomadShadowPrepareRegression(
    runo = prep$runo,
    rord = prep$rord,
    rcon = prep$rcon,
    yvec = prep$ydat,
    mysd = prep$mysd,
    myopti = prep$myopti,
    myoptd = prep$myoptd,
    rbw = prep$rbw,
    penalty.mode = prep$penalty_mode,
    penalty.multiplier = prep$penalty_multiplier,
    degree = prep$degree,
    bernstein = prep$bernstein,
    basis = prep$basis,
    ckerlb = prep$ckerlb,
    ckerub = prep$ckerub
  )
  if (!isTRUE(prepared))
    return(list(
      best_payload = NULL,
      powell.time = NA_real_,
      num.feval.total = 0,
      num.feval.fast.total = 0,
      method = degree.search$engine
    ))
  mpi.barrier(1L)
  on.exit({
    mpi.barrier(1L)
    npRmpiNomadShadowClearRegression()
  }, add = TRUE)
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[ncon + ncat + seq_len(ncon)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)
    flat.bw <- c(bw_vec[template$icon], bw_vec[template$iuno], bw_vec[template$iord])
    out <- npRmpiNomadShadowEvalRegression(
      bw = as.double(flat.bw),
      degree = as.integer(degree)
    )
    nomad.num.feval.total <<- nomad.num.feval.total + 1L
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out[2L])

    list(
      objective = as.numeric(out[1L]),
      degree = degree,
      num.feval = 1L,
      num.feval.fast = as.numeric(out[2L])
    )
  }

  search.engine.used <- if (identical(degree.search$engine, "nomad+powell")) {
    "nomad"
  } else {
    degree.search$engine
  }

  search.result <- .np_nomad_search(
    engine = search.engine.used,
    baseline_record = NULL,
    start_degree = degree.search$start.degree,
    x0 = x0,
    bbin = bbin,
    lb = lb,
    ub = ub,
    eval_fun = eval_fun,
    build_payload = function(point, best_record, solution, interrupted) {
      list(
        payload = NULL,
        objective = as.numeric(best_record$objective),
        powell.time = NA_real_
      )
    },
    direction = "min",
    objective_name = "fval",
    nmulti = nomad.nmulti,
    nomad.inner.nmulti = nomad.inner.nmulti,
    random.seed = random.seed,
    handoff_before_build = identical(degree.search$engine, "nomad+powell"),
    nomad.opts = list(
      DIRECTION_TYPE = "ORTHO 2N",
      QUAD_MODEL_SEARCH = "no",
      NM_SEARCH = "no",
      SPECULATIVE_SEARCH = "no",
      EVAL_OPPORTUNISTIC = "no"
    ),
    degree_spec = list(
      initial = degree.search$start.degree,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      user_supplied = degree.search$start.user
    )
  )
  if (!identical(search.engine.used, degree.search$engine)) {
    search.result$method <- degree.search$engine
  }

  search.result$best_payload <- NULL
  search.result$powell.time <- NA_real_
  search.result$num.feval.total <- as.numeric(nomad.num.feval.total)
  search.result$num.feval.fast.total <- as.numeric(nomad.num.feval.fast.total)
  search.result
}

.npregbw_nomad_search <- function(xdat,
                                  ydat,
                                  bws,
                                  reg.args,
                                  opt.args,
                                  yname,
                                  degree.search,
                                  nomad.inner.nmulti = 0L,
                                  random.seed = 42L) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template <- .npregbw_build_rbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args,
    yname = yname
  )

  if (!identical(template$type, "fixed"))
    stop("automatic degree search with search.engine='nomad' currently requires bwtype='fixed'")

  setup <- .npregbw_nomad_bw_setup(xdat = xdat, template = template)
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(dim(xdat)[2]) else npValidateNmulti(opt.args$nmulti[1L])

  cont_lower <- npGetScaleFactorSearchLower(template,
                                            argname = "template$scale.factor.search.lower")
  bw_lower <- c(rep.int(cont_lower, ncon), rep.int(0, ncat))
  bw_upper <- c(rep.int(1e6, ncon), setup$cat_upper * setup$bandwidth.scale.categorical)

  point.start <- if (all(template$bw == 0)) {
    NULL
  } else {
    .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup)
  }
  x0 <- c(
    .np_nomad_complete_start_point(
      point = point.start,
      lower = bw_lower,
      upper = bw_upper,
      ncont = ncon
    ),
    as.integer(degree.search$start.degree)
  )
  lb <- c(bw_lower, degree.search$lower)
  ub <- c(bw_upper, degree.search$upper)
  bbin <- c(rep.int(0L, ncon + ncat), rep.int(1L, ncon))
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    degree <- as.integer(best_record$degree)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)
    powell.elapsed <- NA_real_
    final.reg.args <- reg.args
    final.reg.args$regtype <- "lp"
    final.reg.args$degree <- degree
    final.reg.args$bernstein.basis <- degree.search$bernstein.basis
    final.tbw <- .npregbw_build_rbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = final.reg.args,
      yname = yname
    )

    direct.payload <- .npregbw_finalize_fixed_degree_payload(
      xdat = xdat,
      ydat = ydat,
      bws = final.tbw,
      core = list(
        bw = as.numeric(final.tbw$bw),
        objective = as.numeric(best_record$objective),
        ifval = as.numeric(best_record$objective),
        num.feval = as.numeric(nomad.num.feval.total),
        num.feval.fast = as.numeric(nomad.num.feval.fast.total),
        fval.history = as.numeric(best_record$objective),
        eval.history = if (isTRUE(nomad.num.feval.total > 0)) rep(1, max(1L, as.integer(nomad.num.feval.total))) else 1,
        invalid.history = 0,
        timing = NA_real_
      ),
      total.time = NA_real_
    )
    if (is.null(direct.payload$timing.profile) && is.list(best_record$timing.profile))
      direct.payload$timing.profile <- best_record$timing.profile
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.opt.args <- opt.args
      hot.opt.args$nmulti <- .np_nomad_powell_hotstart_nmulti("disable_multistart")
      powell.start <- proc.time()[3L]
      hot.payload <- .npregbw_with_powell_refinement_progress(degree, local({
        .npregbw_run_fixed_degree_source_of_truth_collective(
          xdat = xdat,
          ydat = ydat,
          bws = final.tbw,
          opt.args = hot.opt.args
        )
      }))
      powell.elapsed <- proc.time()[3L] - powell.start
      direct.payload$num.feval <- as.numeric(direct.payload$num.feval[1L]) + as.numeric(hot.payload$num.feval[1L])
      direct.payload$num.feval.fast <- as.numeric(direct.payload$num.feval.fast[1L]) + as.numeric(hot.payload$num.feval.fast[1L])
      hot.payload$num.feval <- direct.payload$num.feval
      hot.payload$num.feval.fast <- direct.payload$num.feval.fast
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min")) {
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  if (.npRmpi_has_active_slave_pool(comm = 1L) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    start.bw <- .npregbw_nomad_point_to_bw(x0[seq_len(ncon + ncat)], template = template, setup = setup)
    prep <- .npregbw_nomad_shadow_prepare_args(
      xdat = xdat,
      ydat = ydat,
      bws = template,
      start.bw = start.bw,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier
    )

    mc <- substitute(
      get("npRmpiNomadShadowSearchRegression", envir = asNamespace("npRmpi"), inherits = FALSE)(
        XDAT,
        YDAT,
        TEMPLATE,
        SETUP,
        PREP,
        REGARGS,
        DEGREESEARCH,
        X0,
        BBIN,
        LB,
        UB,
        NOMADNMULTI,
        INNERNMULTI,
        PENMULT,
        RSEED
      ),
      list(
        XDAT = xdat,
        YDAT = ydat,
        TEMPLATE = template,
        SETUP = setup,
        PREP = prep,
        REGARGS = reg.args,
        DEGREESEARCH = degree.search,
        X0 = x0,
        BBIN = bbin,
        LB = lb,
        UB = ub,
        NOMADNMULTI = nomad.nmulti,
        INNERNMULTI = nomad.inner.nmulti,
        PENMULT = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier,
        RSEED = random.seed
      )
    )

    if (isTRUE(.npRmpi_autodispatch_called_from_bcast())) {
      search.result <- eval(mc, envir = environment())
    } else {
      # Prime the first spawned-worker broadcast before the rank-wide NOMAD search.
      .npRmpi_bcast_cmd_expr(quote(invisible(NULL)), comm = 1L, caller.execute = TRUE)
      search.result <- .npRmpi_bcast_cmd_expr(mc, comm = 1L, caller.execute = TRUE)
    }
    if (!is.null(search.result$num.feval.total))
      nomad.num.feval.total <- as.numeric(search.result$num.feval.total[1L])
    if (!is.null(search.result$num.feval.fast.total))
      nomad.num.feval.fast.total <- as.numeric(search.result$num.feval.fast.total[1L])
    best.solution <- NULL
    if (!is.null(search.result$best.restart) &&
        is.finite(search.result$best.restart) &&
        length(search.result$restart.results) >= as.integer(search.result$best.restart)) {
      best.solution <- search.result$restart.results[[as.integer(search.result$best.restart)]]
    }

    payload_result <- build_payload(
      point = search.result$best_point,
      best_record = search.result$best,
      solution = best.solution,
      interrupted = !isTRUE(search.result$completed)
    )

    if (is.list(payload_result) && !is.null(payload_result$payload)) {
      search.result$best_payload <- payload_result$payload
      if (!is.null(payload_result$powell.time))
        search.result$powell.time <- as.numeric(payload_result$powell.time[1L])
      if (!is.null(payload_result$objective) &&
          .np_degree_better(payload_result$objective, search.result$best$objective, direction = "min")) {
        search.result$best$objective <- as.numeric(payload_result$objective[1L])
      }
    } else {
      search.result$best_payload <- payload_result
    }

    return(.npRmpi_reconcile_nomad_search_timing(search.result))
  }

  .np_nomad_baseline_note(degree.search$start.degree)

  start.bw <- .npregbw_nomad_point_to_bw(x0[seq_len(ncon + ncat)], template = template, setup = setup)
  prep <- .npregbw_nomad_shadow_prepare_args(
    xdat = xdat,
    ydat = ydat,
    bws = template,
    start.bw = start.bw,
    invalid.penalty = "baseline",
    penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier
  )

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[ncon + ncat + seq_len(ncon)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)
    flat.bw <- c(bw_vec[template$icon], bw_vec[template$iuno], bw_vec[template$iord])
    out <- npRmpiNomadEvalOnlyRegression(
      runo = prep$runo,
      rord = prep$rord,
      rcon = prep$rcon,
      yvec = prep$ydat,
      mysd = prep$mysd,
      myopti = prep$myopti,
      myoptd = prep$myoptd,
      bw = as.double(flat.bw),
      penalty.mode = prep$penalty_mode,
      penalty.multiplier = prep$penalty_multiplier,
      degree = as.integer(degree),
      bernstein = prep$bernstein,
      basis = prep$basis,
      ckerlb = prep$ckerlb,
      ckerub = prep$ckerub
    )
    nomad.num.feval.total <<- nomad.num.feval.total + 1L
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out$fast.history[1L])

    list(
      objective = as.numeric(out$fval[1L]),
      degree = degree,
      num.feval = 1L,
      num.feval.fast = as.numeric(out$fast.history[1L])
    )
  }

  search.engine.used <- if (identical(degree.search$engine, "nomad+powell")) {
    "nomad"
  } else {
    degree.search$engine
  }

  search.result <- .np_nomad_search(
    engine = search.engine.used,
    baseline_record = baseline.record,
    start_degree = degree.search$start.degree,
    x0 = x0,
    bbin = bbin,
    lb = lb,
    ub = ub,
    eval_fun = eval_fun,
    build_payload = build_payload,
    direction = "min",
    objective_name = "fval",
    nmulti = nomad.nmulti,
    nomad.inner.nmulti = nomad.inner.nmulti,
    random.seed = random.seed,
    handoff_before_build = identical(degree.search$engine, "nomad+powell"),
    nomad.opts = list(
      DIRECTION_TYPE = "ORTHO 2N",
      QUAD_MODEL_SEARCH = "no",
      NM_SEARCH = "no",
      SPECULATIVE_SEARCH = "no",
      EVAL_OPPORTUNISTIC = "no"
    ),
    degree_spec = list(
      initial = degree.search$start.degree,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      user_supplied = degree.search$start.user
    )
  )

  if (!identical(search.engine.used, degree.search$engine)) {
    search.result$method <- degree.search$engine
  }

  search.result
}

.npregbw_degree_search_controls <- function(regtype,
                                            regtype.named,
                                            ncon,
                                            nobs,
                                            basis,
                                            degree.select,
                                            search.engine,
                                            degree.min,
                                            degree.max,
                                            degree.start,
                                            degree.restarts,
                                            degree.max.cycles,
                                            degree.verify,
                                            bernstein.basis,
                                            bernstein.named) {
  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  if (identical(degree.select, "manual"))
    return(NULL)
  search.engine <- .npregbw_nomad_controls(search.engine)

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (ncon < 1L)
    stop("automatic degree search requires at least one continuous regressor")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = ncon,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  baseline.degree <- rep.int(0L, ncon)
  default.start.degree <- if (identical(search.engine, "cell")) {
    baseline.degree
  } else {
    rep.int(1L, ncon)
  }
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, default.start.degree))
  } else {
    start.raw <- npValidateGlpDegree(regtype = "lp", degree = degree.start, ncon = ncon, argname = "degree.start")
    out.of.range <- vapply(seq_len(ncon), function(j) !(start.raw[j] %in% bounds$candidates[[j]]), logical(1))
    if (any(out.of.range))
      stop("degree.start must lie within the searched degree candidates for every continuous predictor")
    start.raw
  }

  list(
    method = if (identical(search.engine, "cell")) degree.select else search.engine,
    engine = search.engine,
    lower = bounds$lower,
    upper = bounds$upper,
    candidates = bounds$candidates,
    baseline.degree = baseline.degree,
    start.degree = start.degree,
    start.user = !is.null(degree.start),
    basis = if (missing(basis) || is.null(basis)) "glp" else as.character(basis[1L]),
    nobs = as.integer(nobs[1L]),
    restarts = npValidateNonNegativeInteger(degree.restarts, "degree.restarts"),
    max.cycles = npValidatePositiveInteger(degree.max.cycles, "degree.max.cycles"),
    verify = npValidateScalarLogical(degree.verify, "degree.verify"),
    bernstein.basis = bern.auto
  )
}

.npregbw_attach_degree_search <- function(bws, search_result) {
  metadata <- list(
    mode = search_result$method,
    verify = isTRUE(search_result$verify),
    completed = isTRUE(search_result$completed),
    certified = isTRUE(search_result$certified),
    interrupted = isTRUE(search_result$interrupted),
    baseline.degree = search_result$baseline$degree,
    baseline.fval = search_result$baseline$objective,
    best.degree = search_result$best$degree,
    best.fval = search_result$best$objective,
    nomad.time = search_result$nomad.time,
    powell.time = search_result$powell.time,
    optim.time = search_result$optim.time,
    n.unique = search_result$n.unique,
    n.visits = search_result$n.visits,
    n.cached = search_result$n.cached,
    grid.size = search_result$grid.size,
    best.restart = search_result$best.restart,
    restart.starts = search_result$restart.starts,
    restart.degree.starts = search_result$restart.degree.starts,
    restart.bandwidth.starts = search_result$restart.bandwidth.starts,
    restart.start.info = search_result$restart.start.info,
    restart.results = search_result$restart.results,
    trace = search_result$trace
  )

  if (!is.null(search_result$nomad.time))
    bws$nomad.time <- as.numeric(search_result$nomad.time[1L])
  if (!is.null(search_result$powell.time))
    bws$powell.time <- as.numeric(search_result$powell.time[1L])
  if (!is.null(search_result$optim.time) && is.finite(search_result$optim.time))
    bws$total.time <- as.numeric(search_result$optim.time[1L])
  if (is.null(bws$timing.profile) && is.list(search_result$best$timing.profile))
    bws$timing.profile <- search_result$best$timing.profile
  bws <- .np_attach_nomad_restart_summary(bws, search_result)
  bws$degree.search <- metadata
  bws
}

npregbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws,
           bandwidth.compute = TRUE,
           basis,
           bernstein.basis,
           bwmethod,
           bwscaling,
           bwtype,
           cfac.dir,
           scale.factor.init,
           ckerbound,
           ckerlb,
           ckerorder,
           ckertype,
           ckerub,
           degree,
           degree.select = c("manual", "coordinate", "exhaustive"),
           search.engine = c("nomad+powell", "cell", "nomad"),
           nomad = FALSE,
           nomad.nmulti = 0L,
           degree.min = NULL,
           degree.max = NULL,
           degree.start = NULL,
           degree.restarts = 0L,
           degree.max.cycles = 20L,
           degree.verify = FALSE,
           dfac.dir,
           dfac.init,
           dfc.dir,
           ftol,
           scale.factor.init.upper,
           hbd.dir,
           hbd.init,
           initc.dir,
           initd.dir,
           invalid.penalty = c("baseline","dbmax"),
           itmax,
           lbc.dir,
           scale.factor.init.lower,
           lbd.dir,
           lbd.init,
           nmulti,
           okertype,
           penalty.multiplier = 10,
           regtype,
           remin,
           scale.init.categorical.sample,
           scale.factor.search.lower = NULL,
           small,
           tol,
           transform.bounds = FALSE,
           ukertype,
           ...){
    .npRmpi_require_active_slave_pool(where = "npregbw()")
    lp.dot.args <- list(...)
    npRejectLegacyLpArgs(names(lp.dot.args), where = "npregbw")
    random.seed.value <- .np_degree_extract_random_seed(lp.dot.args)
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)

    xdat <- toFrame(xdat)
    yname <- deparse(substitute(ydat))

    mc.names <- names(match.call(expand.dots = FALSE))
    nomad.shortcut <- .np_prepare_nomad_shortcut(
      nomad = nomad,
      call_names = mc.names,
      preset = list(
        regtype = "lp",
        search.engine = "nomad+powell",
        degree.select = "coordinate",
        bernstein.basis = TRUE,
        degree.min = 0L,
        degree.max = 10L,
        degree.verify = FALSE,
        bwtype = "fixed"
      ),
      values = list(
        regtype = if ("regtype" %in% mc.names) regtype else NULL,
        search.engine = if ("search.engine" %in% mc.names) search.engine else NULL,
        degree.select = if ("degree.select" %in% mc.names) degree.select else NULL,
        bernstein.basis = if ("bernstein.basis" %in% mc.names) bernstein.basis else NULL,
        degree.min = if ("degree.min" %in% mc.names) degree.min else NULL,
        degree.max = if ("degree.max" %in% mc.names) degree.max else NULL,
        degree.verify = if ("degree.verify" %in% mc.names) degree.verify else NULL,
        bwtype = if ("bwtype" %in% mc.names) bwtype else NULL,
        degree = if ("degree" %in% mc.names) degree else NULL
      ),
      where = "npregbw"
    )

    if (isTRUE(nomad.shortcut$enabled)) {
      if ("degree" %in% mc.names)
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$bwtype, c("fixed", "generalized_nn", "adaptive_nn")))[1L], "fixed"))
        stop("nomad=TRUE currently requires bwtype='fixed'")
      if ("degree.select" %in% mc.names &&
          identical(as.character(match.arg(nomad.shortcut$values$degree.select, c("manual", "coordinate", "exhaustive")))[1L], "manual"))
        stop("nomad=TRUE requires automatic degree search; use degree.select='coordinate' or 'exhaustive'")
      if ("search.engine" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$search.engine, c("nomad+powell", "cell", "nomad")))[1L] %in%
              c("nomad", "nomad+powell")))
        stop("nomad=TRUE requires search.engine='nomad' or 'nomad+powell'")
      if ("bernstein.basis" %in% mc.names &&
          !isTRUE(npValidateGlpBernstein(regtype = "lp",
                                        bernstein.basis = nomad.shortcut$values$bernstein.basis)))
        stop("nomad=TRUE currently requires bernstein.basis=TRUE")
      if ("degree.verify" %in% mc.names &&
          isTRUE(npValidateScalarLogical(nomad.shortcut$values$degree.verify, "degree.verify")))
        stop("nomad=TRUE currently requires degree.verify=FALSE")
    }

    search.mc.names <- mc.names
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    automatic.degree.search <- !identical(match.arg(degree.select.value, c("manual", "coordinate", "exhaustive")), "manual")
    search.engine.value <- if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell"
    regtype.value.early <- if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else "lc"
    nomad.inner <- .np_nomad_validate_inner_multistart(
      call_names = search.mc.names,
      dot.args = lp.dot.args,
      nomad.nmulti = nomad.nmulti,
      regtype = regtype.value.early,
      automatic.degree.search = automatic.degree.search,
      search.engine = search.engine.value
    )
    nomad.inner.named <- nomad.inner$named
    nomad.inner.nmulti <- nomad.inner$nmulti

    if (.npRmpi_autodispatch_active() && !isTRUE(automatic.degree.search))
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_as_generic_call("npregbw", match.call()),
        parent.frame()))

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    margs <- c("regtype", "basis", "degree", "bernstein.basis", "bwmethod", "bwscaling", "bwtype",
               "ckertype", "ckerorder", "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    bw.args <- list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      xnames = names(xdat),
      ynames = yname,
      bandwidth.compute = bandwidth.compute
    )
    if (any.m) {
      nms <- mc.names[m]
      bw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }

    margs <- c("nmulti", "remin", "itmax", "ftol", "tol",
               "small",
               "lbc.dir", "dfc.dir", "cfac.dir","initc.dir",
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.init.categorical.sample",
               "transform.bounds",
               "scale.factor.search.lower",
               "invalid.penalty",
               "penalty.multiplier")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (any.m) {
      nms <- mc.names[m]
      opt.args <- mget(nms, envir = environment(), inherits = FALSE)
    } else {
      opt.args <- list()
    }

    reg.args <- bw.args[setdiff(names(bw.args), c("bw", "nobs", "xdati", "ydati", "xnames", "ynames", "bandwidth.compute"))]
    opt.args <- c(list(bandwidth.compute = bandwidth.compute), opt.args)
    reg.args$scale.factor.search.lower <- scale.factor.search.lower
    opt.args$scale.factor.search.lower <- scale.factor.search.lower

    ncon <- sum(untangle(xdat)$icon)
    regtype.value <- if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else "lc"
    bernstein.value <- if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else TRUE
    degree.min.value <- nomad.shortcut$values$degree.min
    degree.max.value <- nomad.shortcut$values$degree.max
    degree.start.value <- if ("degree.start" %in% search.mc.names) degree.start else NULL
    degree.restarts.value <- if ("degree.restarts" %in% search.mc.names) degree.restarts else 0L
    degree.max.cycles.value <- if ("degree.max.cycles" %in% search.mc.names) degree.max.cycles else 20L
    degree.verify.value <- if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE
    degree.search <- .npregbw_degree_search_controls(
      regtype = regtype.value,
      regtype.named = isTRUE(nomad.shortcut$enabled) || ("regtype" %in% search.mc.names),
      ncon = ncon,
      nobs = NROW(xdat),
      basis = if ("basis" %in% search.mc.names) basis else "glp",
      degree.select = degree.select.value,
      search.engine = search.engine.value,
      degree.min = degree.min.value,
      degree.max = degree.max.value,
      degree.start = degree.start.value,
      degree.restarts = degree.restarts.value,
      degree.max.cycles = degree.max.cycles.value,
      degree.verify = degree.verify.value,
      bernstein.basis = bernstein.value,
      bernstein.named = isTRUE(nomad.shortcut$enabled) || ("bernstein.basis" %in% search.mc.names)
    )
    if (!is.null(degree.search) && is.null(reg.args$degree)) {
      reg.args$degree <- npSetupGlpDegree(
        regtype = regtype.value,
        degree = NULL,
        ncon = ncon,
        degree.select = degree.select.value
      )
    }

    if (!is.null(degree.search)) {
      if (identical(degree.search$engine, "cell")) {
        eval_fun <- function(degree.vec) {
          cell.reg.args <- reg.args
          cell.reg.args$regtype <- "lp"
          cell.reg.args$degree <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
          cell.bws <- .npregbw_run_fixed_degree(
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            reg.args = cell.reg.args,
            opt.args = opt.args,
            yname = yname
          )
          list(
            objective = as.numeric(cell.bws$fval[1L]),
            payload = cell.bws,
            num.feval = if (!is.null(cell.bws$num.feval)) as.numeric(cell.bws$num.feval[1L]) else NA_real_
          )
        }

        search.result <- .np_degree_search(
          method = degree.search$method,
          candidates = degree.search$candidates,
          baseline_degree = degree.search$baseline.degree,
          start_degree = degree.search$start.degree,
          restarts = degree.search$restarts,
          max_cycles = degree.search$max.cycles,
          verify = degree.search$verify,
          eval_fun = eval_fun,
          direction = "min",
          trace_level = "full",
          objective_name = "fval"
        )
      } else {
        search.result <- .npregbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          reg.args = reg.args,
          opt.args = opt.args,
          yname = yname,
          degree.search = degree.search,
          nomad.inner.nmulti = nomad.inner.nmulti,
          random.seed = random.seed.value
        )
      }
      tbw <- .npregbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
    } else {
      tbw <- .npregbw_build_rbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bws,
        bandwidth.compute = bandwidth.compute,
        reg.args = reg.args,
        yname = yname
      )

      bwsel.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting regression bandwidth",
        do.call(npregbw.rbandwidth, bwsel.args)
      )
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)

  }
