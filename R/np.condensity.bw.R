npcdensbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npcdensbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat"),
                                     eval_env = parent.frame())
    UseMethod("npcdensbw", target)
  }

npcdensbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = terms(formula),
                        data = environment(formula),
                        eval_env = environment(formula))
    else .np_terms_ts_mask(terms_obj = terms(formula),
                           data = data,
                           eval_env = environment(formula))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    formula.call <- .np_bw_formula_from_call(call_obj = call, eval_env = parent.frame())
    if (!is.null(formula.call))
      mf[[2]] <- formula.call


    mf[[1]] <- as.name("model.frame")
    formula.obj <- .np_bw_resolve_formula(formula_obj = formula,
                                        formula_call = formula.call,
                                        eval_env = parent.frame())

    variableNames <- explodeFormula(formula.obj)

    ## make formula evaluable, then eval
    varsPlus <- lapply(variableNames, paste, collapse=" + ")
    mf[["formula"]] <- as.formula(paste(" ~ ", varsPlus[[1]]," + ",
                                        varsPlus[[2]]),
                                  env = environment(formula))
    mf[["formula"]] <- terms(mf[["formula"]])
    if(all(orig.ts)){
      args <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      attr(mf[["formula"]], "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
    }

    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())

    ydat <- mf[, variableNames[[1]], drop = FALSE]
    xdat <- mf[, variableNames[[2]], drop = FALSE]

    dots <- list(...)
    tbw <- do.call(npcdensbw, c(list(xdat = xdat, ydat = ydat), dots))

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$variableNames <- variableNames

    tbw
  }

.npcdensbw_assert_bounded_cvls_supported <- function(bws,
                                                     where = "npcdensbw()") {
  method <- if (!is.null(bws$method) && length(bws$method)) {
    as.character(bws$method[1L])
  } else {
    "cv.ml"
  }

  if (!identical(method, "cv.ls"))
    return(invisible(TRUE))

  cykerlb <- if (is.null(bws$cykerlb)) numeric(0L) else bws$cykerlb[bws$iycon]
  cykerub <- if (is.null(bws$cykerub)) numeric(0L) else bws$cykerub[bws$iycon]
  bounded.y <- length(cykerlb) > 0L && any(is.finite(cykerlb) | is.finite(cykerub))

  if (!bounded.y)
    return(invisible(TRUE))

  if (bws$yncon < 1L || bws$yncon > 2L) {
    stop(
      sprintf(
        "%s bounded response cv.ls currently supports up to two continuous response variables with optional ordered/unordered discrete response components",
        where
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.npcdensbw_validate_scale_factor_lower_bound <- function(value,
                                                         argname = "scale.factor.search.lower") {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value < 0) {
    stop(sprintf("%s must be a single nonnegative finite numeric value", argname),
         call. = FALSE)
  }

  as.double(value)
}

.npcdensbw_resolve_scale_factor_lower_bound <- function(value,
                                                        fallback = 0.1,
                                                        argname = "scale.factor.search.lower") {
  if (is.null(value))
    return(as.double(fallback))

  .npcdensbw_validate_scale_factor_lower_bound(value, argname = argname)
}

.npcdensbw_validate_cvls_quadrature_extend_factor <- function(value,
                                                              argname = "cvls.quadrature.extend.factor") {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value <= 0) {
    stop(sprintf("%s must be a single positive finite numeric value", argname),
         call. = FALSE)
  }

  as.double(value)
}

.npcdensbw_resolve_cvls_quadrature_extend_factor <- function(value,
                                                             fallback = 1,
                                                             argname = "cvls.quadrature.extend.factor") {
  if (is.null(value))
    return(as.double(fallback))

  .npcdensbw_validate_cvls_quadrature_extend_factor(value, argname = argname)
}

.npcdensbw_validate_cvls_quadrature_points <- function(value,
                                                       argname = "cvls.quadrature.points") {
  if (length(value) != 2L || !is.numeric(value) || anyNA(value) ||
      any(!is.finite(value)) || any(value < 2) ||
      any(abs(value - round(value)) > sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a two-element finite integer vector with entries at least 2", argname),
         call. = FALSE)
  }

  as.integer(round(value))
}

.npcdensbw_resolve_cvls_quadrature_points <- function(value,
                                                      fallback = c(100L, 50L),
                                                      argname = "cvls.quadrature.points") {
  if (is.null(value))
    return(as.integer(fallback))

  .npcdensbw_validate_cvls_quadrature_points(value, argname = argname)
}

.npcdensbw_validate_cvls_quadrature_ratios <- function(value,
                                                       argname = "cvls.quadrature.ratios") {
  if (length(value) != 3L || !is.numeric(value) || anyNA(value) ||
      any(!is.finite(value)) || any(value < 0) ||
      !isTRUE(all.equal(sum(value), 1, tolerance = 1e-8))) {
    stop(sprintf("%s must be a three-element non-negative numeric vector summing to one",
                 argname),
         call. = FALSE)
  }

  as.double(value)
}

.npcdensbw_resolve_cvls_quadrature_ratios <- function(value,
                                                      fallback = c(0.20, 0.55, 0.25),
                                                      argname = "cvls.quadrature.ratios") {
  if (is.null(value))
    return(as.double(fallback))

  .npcdensbw_validate_cvls_quadrature_ratios(value, argname = argname)
}

.npcdensbw_resolve_cvls_quadrature_grid <- function(value,
                                                    fallback = "hybrid",
                                                    argname = "cvls.quadrature.grid") {
  if (is.null(value))
    value <- fallback

  if (length(value) != 1L || is.na(value))
    stop(sprintf("%s must be one of 'hybrid', 'uniform', or 'sample'", argname),
         call. = FALSE)

  value <- as.character(value)
  if (!value %in% c("hybrid", "uniform", "sample"))
    stop(sprintf("%s must be one of 'hybrid', 'uniform', or 'sample'", argname),
         call. = FALSE)

  value
}

.npcdensbw_cvls_quadrature_grid_fallback <- function(yncon) {
  if (as.integer(yncon) >= 2L) "uniform" else "hybrid"
}

.npcdensbw_validate_cvls_quadrature_grid_dimension <- function(value,
                                                               yncon,
                                                               argname = "cvls.quadrature.grid") {
  if (as.integer(yncon) >= 2L && !identical(value, "uniform")) {
    stop(sprintf("%s = '%s' is currently supported only for scalar continuous responses",
                 argname, value),
         call. = FALSE)
  }

  value
}

.npcdensbw_cvls_quadrature_grid_code <- function(value) {
  switch(.npcdensbw_resolve_cvls_quadrature_grid(value),
         uniform = 0L,
         hybrid = 1L,
         sample = 2L)
}


.npcdensbw_effective_cvls_quadrature_points <- function(points, yncon) {
  points <- .npcdensbw_resolve_cvls_quadrature_points(points)
  if (as.integer(yncon) >= 2L) points[[2L]] else points[[1L]]
}

.npcdensbw_apply_continuous_search_floor <- function(tbw,
                                                     xdat,
                                                     ydat,
                                                     scale.factor.search.lower) {
  tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)

  if (!isTRUE(tbw$bandwidth.compute) || !identical(tbw$type, "fixed") || tbw$ncon < 1L)
    return(tbw)

  if (isTRUE(tbw$scaling)) {
    floor.values <- rep(as.double(scale.factor.search.lower), tbw$ncon)
  } else {
    xcon <- xdat[, tbw$ixcon, drop = FALSE]
    ycon <- ydat[, tbw$iycon, drop = FALSE]
    mysd <- EssDee(data.frame(xcon, ycon))
    nconfac <- nrow(xdat)^(-1.0 / (2.0 * tbw$cxkerorder + tbw$ncon))
    floor.values <- as.double(scale.factor.search.lower) * mysd * nconfac
  }

  if (tbw$xncon > 0L) {
    tbw$xbw[tbw$ixcon] <- pmax(tbw$xbw[tbw$ixcon], floor.values[seq_len(tbw$xncon)])
  }
  if (tbw$yncon > 0L) {
    idx <- tbw$xncon + seq_len(tbw$yncon)
    tbw$ybw[tbw$iycon] <- pmax(tbw$ybw[tbw$iycon], floor.values[idx])
  }

  tbw
}

npcdensbw.conbandwidth <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
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
           memfac = 500.0,
           nmulti,
           penalty.multiplier = 10,
           remin = TRUE,
           scale.init.categorical.sample = FALSE,
           scale.factor.search.lower = NULL,
           cvls.quadrature.grid = NULL,
           cvls.quadrature.extend.factor = NULL,
           cvls.quadrature.points = NULL,
           cvls.quadrature.ratios = NULL,
           small = 1.490116e-05,
           tol = 1.490116e-04,
           transform.bounds = FALSE,
           ...){
    elapsed.start <- proc.time()[3]
    ydat = toFrame(ydat)
    xdat = toFrame(xdat)

    mc.expanded <- match.call(expand.dots = TRUE)
    if ("cvls.i1.rescue" %in% names(mc.expanded))
      stop("cvls.i1.rescue has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive.tol" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive.tol has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive.grid.hy.ratio" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive.grid.hy.ratio has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive.floor.tol" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive.floor.tol has been removed; use cvls.quadrature.grid",
           call. = FALSE)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(dim(ydat)[2]+dim(xdat)[2])
    }
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    remin <- npValidateScalarLogical(remin, "remin")
    scale.init.categorical.sample <-
      npValidateScalarLogical(scale.init.categorical.sample, "scale.init.categorical.sample")
    scale.factor.search.lower <- .npcdensbw_resolve_scale_factor_lower_bound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower,
      fallback = 0.1,
      argname = "scale.factor.search.lower"
    )
    cvls.quadrature.grid <- .npcdensbw_resolve_cvls_quadrature_grid(
      if (is.null(cvls.quadrature.grid)) bws$cvls.quadrature.grid else cvls.quadrature.grid,
      fallback = .npcdensbw_cvls_quadrature_grid_fallback(bws$yncon),
      argname = "cvls.quadrature.grid"
    )
    cvls.quadrature.grid <- .npcdensbw_validate_cvls_quadrature_grid_dimension(
      cvls.quadrature.grid,
      bws$yncon,
      argname = "cvls.quadrature.grid"
    )
    cvls.quadrature.extend.factor <- .npcdensbw_resolve_cvls_quadrature_extend_factor(
      if (is.null(cvls.quadrature.extend.factor)) bws$cvls.quadrature.extend.factor else cvls.quadrature.extend.factor,
      fallback = 1,
      argname = "cvls.quadrature.extend.factor"
    )
    cvls.quadrature.points <- .npcdensbw_resolve_cvls_quadrature_points(
      if (is.null(cvls.quadrature.points)) bws$cvls.quadrature.points else cvls.quadrature.points,
      fallback = c(100L, 50L),
      argname = "cvls.quadrature.points"
    )
    cvls.quadrature.ratios <- .npcdensbw_resolve_cvls_quadrature_ratios(
      if (is.null(cvls.quadrature.ratios)) bws$cvls.quadrature.ratios else cvls.quadrature.ratios,
      fallback = c(0.20, 0.55, 0.25),
      argname = "cvls.quadrature.ratios"
    )
    transform.bounds <- npValidateScalarLogical(transform.bounds, "transform.bounds")
    bws <- npSetScaleFactorSearchLower(bws, scale.factor.search.lower)
    bws$cvls.quadrature.grid <- cvls.quadrature.grid
    bws$cvls.quadrature.extend.factor <- cvls.quadrature.extend.factor
    bws$cvls.quadrature.points <- cvls.quadrature.points
    bws$cvls.quadrature.ratios <- cvls.quadrature.ratios
    itmax <- npValidatePositiveInteger(itmax, "itmax")
    ftol <- npValidatePositiveFiniteNumeric(ftol, "ftol")
    tol <- npValidatePositiveFiniteNumeric(tol, "tol")
    small <- npValidatePositiveFiniteNumeric(small, "small")
    memfac <- npValidatePositiveFiniteNumeric(memfac, "memfac")
    penalty.multiplier <- npValidatePositiveFiniteNumeric(penalty.multiplier, "penalty.multiplier")
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)
    if (length(bws$ybw) != dim(ydat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of", "'ydat'"))

    if (length(bws$xbw) != dim(xdat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of", "'xdat'"))

    if (dim(ydat)[1] != dim(xdat)[1])
      stop(paste("number of rows of", "'ydat'", "does not match", "'xdat'"))

    if ((any(bws$iycon) &&
         !all(vapply(as.data.frame(ydat[, bws$iycon]), inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iyord) &&
         !all(vapply(as.data.frame(ydat[, bws$iyord]), inherits, logical(1), "ordered"))) ||
        (any(bws$iyuno) &&
         !all(vapply(as.data.frame(ydat[, bws$iyuno]), inherits, logical(1), "factor"))))
      stop(paste("supplied bandwidths do not match", "'ydat'", "in type"))

    if ((any(bws$ixcon) &&
         !all(vapply(as.data.frame(xdat[, bws$ixcon]), inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$ixord) &&
         !all(vapply(as.data.frame(xdat[, bws$ixord]), inherits, logical(1), "ordered"))) ||
        (any(bws$ixuno) &&
         !all(vapply(as.data.frame(xdat[, bws$ixuno]), inherits, logical(1), "factor"))))
      stop(paste("supplied bandwidths do not match", "'xdat'", "in type"))

    ## catch and destroy NA's
    goodrows <- seq_len(nrow(xdat))
    rows.omit <- unclass(na.action(na.omit(data.frame(xdat,ydat))))
    goodrows[rows.omit] <- 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    spec <- npCanonicalConditionalRegSpec(
      regtype = if (is.null(bws$regtype)) "lc" else bws$regtype,
      basis = if (is.null(bws$basis)) "glp" else bws$basis,
      degree = if (is.null(bws$degree)) NULL else bws$degree,
      bernstein.basis = isTRUE(bws$bernstein.basis),
      ncon = bws$xncon,
      where = "npcdensbw"
    )
    .npcdensbw_assert_bounded_cvls_supported(bws, where = "npcdensbw()")
    .npRmpi_require_active_slave_pool(where = "npcdensbw()")
    keep_local_shadow_nn <- bandwidth.compute &&
      identical(spec$regtype.engine, "lp") &&
      identical(bws$method %in% c("cv.ml", "cv.ls"), TRUE) &&
      identical(bws$type %in% c("generalized_nn", "adaptive_nn"), TRUE)
    keep_local_raw_degree1_cvls <- bandwidth.compute &&
      identical(bws$method, "cv.ls") &&
      identical(bws$type, "fixed") &&
      npIsRawDegreeOneConditionalSpec(spec, bws$xncon)
    if (.npRmpi_autodispatch_active() &&
        !keep_local_shadow_nn &&
        !keep_local_raw_degree1_cvls)
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows,,drop = FALSE]


    nrow = nrow(ydat)
    yncol = ncol(ydat)
    xncol = ncol(xdat)

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    ydat = toMatrix(ydat)

    yuno = ydat[, bws$iyuno, drop = FALSE]
    ycon = ydat[, bws$iycon, drop = FALSE]
    yord = ydat[, bws$iyord, drop = FALSE]


    xdat = toMatrix(xdat)

    xuno = xdat[, bws$ixuno, drop = FALSE]
    xcon = xdat[, bws$ixcon, drop = FALSE]
    xord = xdat[, bws$ixord, drop = FALSE]

    tbw <- bws
    spec <- npCanonicalConditionalRegSpec(
      regtype = if (is.null(tbw$regtype)) "lc" else tbw$regtype,
      basis = if (is.null(tbw$basis)) "glp" else tbw$basis,
      degree = if (is.null(tbw$degree)) NULL else tbw$degree,
      bernstein.basis = isTRUE(tbw$bernstein.basis),
      ncon = tbw$xncon,
      where = "npcdensbw"
    )
    tbw$regtype <- spec$regtype
    tbw$pregtype <- switch(spec$regtype,
                           lc = "Local-Constant",
                           ll = "Local-Linear",
                           lp = "Local-Polynomial")
    tbw$basis <- spec$basis
    tbw$degree <- spec$degree
    tbw$bernstein.basis <- spec$bernstein.basis
    tbw$regtype.engine <- spec$regtype.engine
    tbw$basis.engine <- spec$basis.engine
    tbw$degree.engine <- spec$degree.engine
    tbw$bernstein.basis.engine <- spec$bernstein.basis.engine
    reg.code <- if (identical(spec$regtype.engine, "lp")) REGTYPE_LP else REGTYPE_LC
    degree.code <- if (tbw$xncon > 0L) as.integer(spec$degree.engine) else integer(0)
    basis.code <- as.integer(npLpBasisCode(spec$basis.engine))
    bernstein.engine <- isTRUE(spec$bernstein.basis.engine)

    mysd <- EssDee(data.frame(xcon,ycon))
    nconfac <- nrow^(-1.0/(2.0*bws$cxkerorder+bws$ncon))
    ncatfac <- nrow^(-2.0/(2.0*bws$cxkerorder+bws$ncon))

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

    if (bandwidth.compute){
      cont.start <- npContinuousSearchStartControls(
        scale.factor.init.lower,
        scale.factor.init.upper,
        scale.factor.init,
        tbw$scale.factor.search.lower,
        where = "npcdensbw"
      )
      myopti = list(num_obs_train = nrow,
        iMultistart = IMULTI_TRUE,
        iNum_Multistart = nmulti,
        int_use_starting_values = (if (all(bws$ybw==0) && all(bws$xbw==0)) USE_START_NO else USE_START_YES),
        int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
        BANDWIDTH_den_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=(if (remin) RE_MIN_TRUE else RE_MIN_FALSE),
        int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
        bwmethod = switch(bws$method,
          cv.ml = CBWM_CVML,
          cv.ls = CBWM_CVLS),
        xkerneval = switch(bws$cxkertype,
          gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
          uniform = CKER_UNI,
          "truncated gaussian" = CKER_TGAUSS),
        ykerneval = switch(bws$cykertype,
          gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
          uniform = CKER_UNI,
          "truncated gaussian" = CKER_TGAUSS),
        uxkerneval = switch(bws$uxkertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        uykerneval = switch(bws$uykertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        oxkerneval = switch(bws$oxkertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_LR,
        "racineliyan" = OKER_RLY),
        oykerneval = switch(bws$oykertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_NLR,
        "racineliyan" = OKER_RLY),
        ynuno = dim(yuno)[2],
        ynord = dim(yord)[2],
        yncon = dim(ycon)[2],
        xnuno = dim(xuno)[2],
        xnord = dim(xord)[2],
        xncon = dim(xcon)[2],
        old.cdens = FALSE,
        int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
        scale.init.categorical.sample = scale.init.categorical.sample,
        dfc.dir = dfc.dir,
        transform.bounds = transform.bounds,
        cvls.quadrature.grid = .npcdensbw_cvls_quadrature_grid_code(tbw$cvls.quadrature.grid),
        cvls.quadrature.points =
          .npcdensbw_effective_cvls_quadrature_points(tbw$cvls.quadrature.points, tbw$yncon))

      myoptd = list(ftol=ftol, tol=tol, small=small, memfac = memfac,
        lbc.dir = lbc.dir, cfac.dir = cfac.dir, initc.dir = initc.dir,
        lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir,
        lbc.init = cont.start$scale.factor.init.lower,
        hbc.init = cont.start$scale.factor.init.upper,
        cfac.init = cont.start$scale.factor.init,
        lbd.init = lbd.init, hbd.init = hbd.init, dfac.init = dfac.init,
        nconfac = nconfac, ncatfac = ncatfac,
        scale.factor.lower.bound = tbw$scale.factor.search.lower,
        cvls.quadrature.extend.factor = tbw$cvls.quadrature.extend.factor,
        cvls.quadrature.ratios.uniform = tbw$cvls.quadrature.ratios[[1L]],
        cvls.quadrature.ratios.sample = tbw$cvls.quadrature.ratios[[2L]],
        cvls.quadrature.ratios.gl = tbw$cvls.quadrature.ratios[[3L]])

      cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
      cyker.bounds.c <- .npcdensbw_marshal_y_bounds(bws$cykerlb[bws$iycon],
                                                    bws$cykerub[bws$iycon],
                                                    bws$cykerbound)

      if (bws$method != "normal-reference"){
        myout <- npWithLocalLinearRawBasisSearchError(
          if (keep_local_shadow_nn || keep_local_raw_degree1_cvls) {
            .npRmpi_with_local_regression(
              .Call("C_np_density_conditional_bw",
                    as.double(yuno), as.double(yord), as.double(ycon),
                    as.double(xuno), as.double(xord), as.double(xcon),
                    as.double(mysd),
                    as.integer(myopti), as.double(myoptd),
                    as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                                bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                                bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
                    as.integer(nmulti),
                    as.integer(penalty_mode),
                    as.double(penalty.multiplier),
                    as.integer(degree.code),
                    as.integer(bernstein.engine),
                    as.integer(basis.code),
                    as.integer(reg.code),
                    as.double(cxker.bounds.c$lb),
                    as.double(cxker.bounds.c$ub),
                    as.double(cyker.bounds.c$lb),
                    as.double(cyker.bounds.c$ub),
                    PACKAGE="npRmpi")
            )
          } else {
            .Call("C_np_density_conditional_bw",
                  as.double(yuno), as.double(yord), as.double(ycon),
                  as.double(xuno), as.double(xord), as.double(xcon),
                  as.double(mysd),
                  as.integer(myopti), as.double(myoptd),
                  as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                              bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                              bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
                  as.integer(nmulti),
                  as.integer(penalty_mode),
                  as.double(penalty.multiplier),
                  as.integer(degree.code),
                  as.integer(bernstein.engine),
                  as.integer(basis.code),
                  as.integer(reg.code),
                  as.double(cxker.bounds.c$lb),
                  as.double(cxker.bounds.c$ub),
                  as.double(cyker.bounds.c$lb),
                  as.double(cyker.bounds.c$ub),
                  PACKAGE="npRmpi")
          },
          where = "npcdensbw",
          spec = spec,
          bwmethod = bws$method,
          ncon = tbw$xncon
        )
        total.time <- proc.time()[3] - elapsed.start
      } else {
        nbw = double(yncol+xncol)
        gbw = bws$yncon+bws$xncon
        if (gbw > 0){
          gbw_idx <- seq_len(gbw)
          nbw[gbw_idx] = 1.059224
          if(!bws$scaling)
            nbw[gbw_idx]=nbw[gbw_idx]*mysd*nconfac
        }
        myout= list( bw = nbw, fval = c(NA,NA,NA) )
        total.time <- NA
      }

      yr = seq_len(yncol)
      xr = seq_len(xncol)
      rorder = numeric(yncol + xncol)

      ## bandwidths are passed back from the C routine in an unusual order
      ## xc, y[cuo], x[uo]

      rxcon = xr[bws$ixcon]
      rxuno = xr[bws$ixuno]
      rxord = xr[bws$ixord]

      rycon = yr[bws$iycon]
      ryuno = yr[bws$iyuno]
      ryord = yr[bws$iyord]


      ## rorder[c(rxcon,rycon,ryuno,ryord,rxuno,rxord)]=1:(yncol+xncol)

      tbw <- bws
      tbw$ybw[c(rycon,ryuno,ryord)] <- myout$bw[yr+bws$xncon]
      tbw$xbw[c(rxcon,rxuno,rxord)] <- myout$bw[setdiff(seq_len(yncol + xncol), yr + bws$xncon)]

      tbw$fval = myout$fval[1]
      tbw$ifval = myout$fval[2]
      tbw$initial.fval <- if (length(myout$fval) >= 3L) myout$fval[3L] else NA_real_
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
      tbw$fval.history <- myout$fval.history
      tbw$eval.history <- myout$eval.history
      tbw$invalid.history <- myout$invalid.history
      tbw$timing <- myout$timing
      tbw$total.time <- total.time
    }
    tbw <- npSetScaleFactorSearchLower(tbw, npGetScaleFactorSearchLower(bws))

    ## bandwidth metadata
    tbw$sfactor <- tbw$bandwidth <- list(x = tbw$xbw, y = tbw$ybw)

    apply_bw_meta <- function(tl, dfactor){
      for (nm in names(tl)) {
        idx <- tl[[nm]]
        if (length(idx) == 0L)
          next
        if (tbw$scaling) {
          tbw$bandwidth[[nm]][idx] <- tbw$bandwidth[[nm]][idx] * dfactor[[nm]]
        } else {
          tbw$sfactor[[nm]][idx] <- tbw$sfactor[[nm]][idx] / dfactor[[nm]]
        }
      }
    }

    if ((tbw$xnuno+tbw$ynuno) > 0){
      dfactor <- ncatfac
      dfactor <- list(x = dfactor, y = dfactor)

      tl <- list(x = tbw$xdati$iuno, y = tbw$ydati$iuno)

      apply_bw_meta(tl = tl, dfactor = dfactor)
    }

    if ((tbw$xnord+tbw$ynord) > 0){
      dfactor <- ncatfac
      dfactor <- list(x = dfactor, y = dfactor)

      tl <- list(x = tbw$xdati$iord, y = tbw$ydati$iord)

      apply_bw_meta(tl = tl, dfactor = dfactor)
    }


    if (tbw$ncon > 0){
      dfactor <- nconfac
      dfactor <- list(x = EssDee(xcon)*dfactor, y = EssDee(ycon)*dfactor)

      tl <- list(x = tbw$xdati$icon, y = tbw$ydati$icon)

      apply_bw_meta(tl = tl, dfactor = dfactor)
    }

    initial.fval <- tbw$initial.fval
    tbw <- conbandwidth(xbw = tbw$xbw,
                        ybw = tbw$ybw,
                        bwmethod = tbw$method,
                        bwscaling = tbw$scaling,
                        bwtype = tbw$type,
                        cxkertype = tbw$cxkertype,
                        cxkerorder = tbw$cxkerorder,
                        cxkerbound = tbw$cxkerbound,
                        cxkerlb = tbw$cxkerlb[tbw$ixcon],
                        cxkerub = tbw$cxkerub[tbw$ixcon],
                        uxkertype = tbw$uxkertype,
                        oxkertype = tbw$oxkertype,
                        cykertype = tbw$cykertype,
                        cykerorder = tbw$cykerorder,
                        cykerbound = tbw$cykerbound,
                        cykerlb = tbw$cykerlb[tbw$iycon],
                        cykerub = tbw$cykerub[tbw$iycon],
                        uykertype = tbw$uykertype,
                        oykertype = tbw$oykertype,
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
                        total.time = tbw$total.time,
                        regtype = tbw$regtype,
                        pregtype = tbw$pregtype,
                        basis = tbw$basis,
                        degree = tbw$degree,
                        bernstein.basis = tbw$bernstein.basis,
                        regtype.engine = tbw$regtype.engine,
                        basis.engine = tbw$basis.engine,
                        degree.engine = tbw$degree.engine,
                        bernstein.basis.engine = tbw$bernstein.basis.engine)
    tbw <- .np_refresh_xy_bandwidth_metadata(tbw)
    tbw$initial.fval <- if (!is.null(initial.fval)) initial.fval else NA_real_
    tbw <- npSetScaleFactorSearchLower(tbw, npGetScaleFactorSearchLower(bws))
    tbw$cvls.quadrature.grid <- bws$cvls.quadrature.grid
    tbw$cvls.quadrature.extend.factor <- bws$cvls.quadrature.extend.factor
    tbw$cvls.quadrature.points <- bws$cvls.quadrature.points
    tbw$cvls.quadrature.ratios <- bws$cvls.quadrature.ratios
    tbw <- .npcdensbw_restore_explicit_fixed_y_bounds(tbw, bws)

    tbw
  }

.npcdensbw_build_conbandwidth <- function(xdat,
                                          ydat,
                                          bws,
                                          bandwidth.compute,
                                          reg.args) {
  x.info <- untangle(xdat)
  y.info <- untangle(ydat)
  y.idx <- seq_len(ncol(ydat))
  x.idx <- seq_len(ncol(xdat))

  bw.args <- c(
    list(
      xbw = bws[length(y.idx) + x.idx],
      ybw = bws[y.idx],
      nobs = nrow(xdat),
      xdati = x.info,
      ydati = y.info,
      xnames = names(xdat),
      ynames = names(ydat),
      bandwidth.compute = bandwidth.compute
    ),
    reg.args
  )

  tbw <- do.call(conbandwidth, bw.args)
  tbw$scale.factor.search.lower <- reg.args$scale.factor.search.lower
  tbw$cvls.quadrature.grid <- reg.args$cvls.quadrature.grid
  tbw$cvls.quadrature.extend.factor <- reg.args$cvls.quadrature.extend.factor
  tbw$cvls.quadrature.points <- reg.args$cvls.quadrature.points
  tbw$cvls.quadrature.ratios <- reg.args$cvls.quadrature.ratios
  tbw <- .npcdensbw_apply_continuous_search_floor(
    tbw,
    xdat = xdat,
    ydat = ydat,
    scale.factor.search.lower = reg.args$scale.factor.search.lower
  )
  .npcdensbw_restore_explicit_fixed_y_bounds(tbw, reg.args)
}

.npcdensbw_run_fixed_degree <- function(xdat, ydat, bws, reg.args, opt.args) {
  tbw <- .npcdensbw_build_conbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args
  )

  do.call(npcdensbw.conbandwidth, c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args))
}

.npcdensbw_is_explicit_fixed_all_infinite <- function(kerlb, kerub, kerbound) {
  if (is.null(kerbound) || !identical(as.character(kerbound)[1L], "fixed"))
    return(FALSE)
  if (is.null(kerlb) || is.null(kerub))
    return(FALSE)

  lb <- as.double(kerlb)
  ub <- as.double(kerub)

  length(lb) > 0L &&
    length(ub) > 0L &&
    all(is.infinite(lb) & lb < 0) &&
    all(is.infinite(ub) & ub > 0)
}

.npcdensbw_has_explicit_fixed_infinite_y_bound <- function(kerlb, kerub, kerbound) {
  if (is.null(kerbound) || !any(as.character(kerbound) == "fixed", na.rm = TRUE))
    return(FALSE)
  if (is.null(kerlb) || is.null(kerub))
    return(FALSE)

  lb <- as.double(kerlb)
  ub <- as.double(kerub)

  (length(lb) > 0L && any(is.infinite(lb))) ||
    (length(ub) > 0L && any(is.infinite(ub)))
}

.npcdensbw_warn_infinite_response_quadrature <- function(kerlb,
                                                         kerub,
                                                         kerbound,
                                                         points.supplied,
                                                         where = "npcdensbw()") {
  if (isTRUE(points.supplied))
    return(invisible(FALSE))

  if (!.npcdensbw_has_explicit_fixed_infinite_y_bound(kerlb, kerub, kerbound))
    return(invisible(FALSE))

  warning(sprintf(
    "%s with fixed infinite response bounds uses a finite cv.ls quadrature surrogate; consider setting cvls.quadrature.points explicitly for this edge case",
    where
  ), call. = FALSE)
  invisible(TRUE)
}

.npcdensbw_restore_explicit_fixed_y_bounds <- function(tbw, reg.args) {
  if (!.npcdensbw_is_explicit_fixed_all_infinite(reg.args$cykerlb,
                                                 reg.args$cykerub,
                                                 reg.args$cykerbound)) {
    return(tbw)
  }
  if (!identical(tbw$cykerbound, "none") || !any(tbw$iycon))
    return(tbw)

  icon.idx <- which(tbw$iycon)
  pyorder <- switch(tbw$cykerorder / 2,
                    "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order")

  tbw$cykerbound <- "fixed"
  tbw$cykerlb[icon.idx] <- -Inf
  tbw$cykerub[icon.idx] <- Inf
  tbw$pcykertype <- cktToPrint(tbw$cykertype, order = pyorder, kerbound = "fixed")

  if (!is.null(tbw$klist$y)) {
    tbw$klist$y$ckerbound <- "fixed"
    tbw$klist$y$ckerlb <- tbw$cykerlb
    tbw$klist$y$ckerub <- tbw$cykerub
    tbw$klist$y$pckertype <- tbw$pcykertype
  }

  tbw
}

.npcdensbw_marshal_y_bounds <- function(kerlb, kerub, kerbound) {
  if (.npcdensbw_is_explicit_fixed_all_infinite(kerlb, kerub, kerbound)) {
    sentinel <- .Machine$double.xmax / 4
    n <- max(length(kerlb), length(kerub))
    return(list(lb = rep.int(-sentinel, n), ub = rep.int(sentinel, n)))
  }

  npKernelBoundsMarshal(kerlb, kerub)
}

.npcdensbw_eval_only <- function(xdat,
                                 ydat,
                                 bws,
                                 invalid.penalty = c("baseline", "dbmax"),
                                 penalty.multiplier = 10) {
  invalid.penalty <- match.arg(invalid.penalty)

  ydat <- toFrame(ydat)
  xdat <- toFrame(xdat)

  if (length(bws$ybw) != dim(ydat)[2])
    stop("length of bandwidth vector does not match number of columns of 'ydat'")
  if (length(bws$xbw) != dim(xdat)[2])
    stop("length of bandwidth vector does not match number of columns of 'xdat'")
  if (dim(ydat)[1] != dim(xdat)[1])
    stop("number of rows of 'ydat' does not match 'xdat'")

  goodrows <- seq_len(nrow(xdat))
  rows.omit <- unclass(na.action(na.omit(data.frame(xdat, ydat))))
  goodrows[rows.omit] <- 0

  xdat <- xdat[goodrows,, drop = FALSE]
  ydat <- ydat[goodrows,, drop = FALSE]

  ymat <- toMatrix(ydat)
  xmat <- toMatrix(xdat)

  yuno <- ymat[, bws$iyuno, drop = FALSE]
  ycon <- ymat[, bws$iycon, drop = FALSE]
  yord <- ymat[, bws$iyord, drop = FALSE]
  xuno <- xmat[, bws$ixuno, drop = FALSE]
  xcon <- xmat[, bws$ixcon, drop = FALSE]
  xord <- xmat[, bws$ixord, drop = FALSE]

  mysd <- EssDee(data.frame(xcon, ycon))
  nrow <- nrow(ymat)
  nconfac <- nrow^(-1.0 / (2.0 * bws$cxkerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$cxkerorder + bws$ncon))
  penalty_mode <- if (invalid.penalty == "baseline") 1L else 0L
  scale.factor.search.lower <- .npcdensbw_resolve_scale_factor_lower_bound(
    npGetScaleFactorSearchLower(bws),
    fallback = 0.1,
    argname = "bws$scale.factor.search.lower"
  )
  cvls.quadrature.extend.factor <- .npcdensbw_resolve_cvls_quadrature_extend_factor(
    bws$cvls.quadrature.extend.factor,
    fallback = 1,
    argname = "bws$cvls.quadrature.extend.factor"
  )
  cvls.quadrature.points <- .npcdensbw_resolve_cvls_quadrature_points(
    bws$cvls.quadrature.points,
    fallback = c(100L, 50L),
    argname = "bws$cvls.quadrature.points"
  )
  cvls.quadrature.ratios <- .npcdensbw_resolve_cvls_quadrature_ratios(
    bws$cvls.quadrature.ratios,
    fallback = c(0.20, 0.55, 0.25),
    argname = "bws$cvls.quadrature.ratios"
  )
  cvls.quadrature.grid <- .npcdensbw_resolve_cvls_quadrature_grid(
    bws$cvls.quadrature.grid,
    fallback = .npcdensbw_cvls_quadrature_grid_fallback(bws$yncon),
    argname = "bws$cvls.quadrature.grid"
  )
  cvls.quadrature.grid <- .npcdensbw_validate_cvls_quadrature_grid_dimension(
    cvls.quadrature.grid,
    bws$yncon,
    argname = "bws$cvls.quadrature.grid"
  )

  penalty_mode <- if (invalid.penalty == "baseline") 1L else 0L
  reg.code <- if (identical(bws$regtype.engine, "lp")) REGTYPE_LP else REGTYPE_LC
  degree.code <- if (bws$xncon > 0L) as.integer(bws$degree.engine) else integer(0L)
  basis.code <- as.integer(npLpBasisCode(bws$basis.engine))
  bernstein.engine <- isTRUE(bws$bernstein.basis.engine)

  .npcdensbw_assert_bounded_cvls_supported(bws, where = ".npcdensbw_eval_only()")

  myopti <- list(
    num_obs_train = nrow,
    iMultistart = IMULTI_FALSE,
    iNum_Multistart = 0L,
    int_use_starting_values = USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = 0L,
    int_RESTART_FROM_MIN = RE_MIN_FALSE,
    int_MINIMIZE_IO = IO_MIN_TRUE,
    bwmethod = switch(bws$method,
      cv.ml = CBWM_CVML,
      cv.ls = CBWM_CVLS),
    xkerneval = switch(bws$cxkertype,
      gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS),
    ykerneval = switch(bws$cykertype,
      gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS),
    uxkerneval = switch(bws$uxkertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    uykerneval = switch(bws$uykertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    oxkerneval = switch(bws$oxkertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR,
      racineliyan = OKER_RLY),
    oykerneval = switch(bws$oykertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY),
    ynuno = dim(yuno)[2],
    ynord = dim(yord)[2],
    yncon = dim(ycon)[2],
    xnuno = dim(xuno)[2],
    xnord = dim(xord)[2],
    xncon = dim(xcon)[2],
    old.cdens = FALSE,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    scale.init.categorical.sample = FALSE,
    dfc.dir = 0L,
    transform.bounds = FALSE,
    cvls.quadrature.grid = .npcdensbw_cvls_quadrature_grid_code(cvls.quadrature.grid),
    cvls.quadrature.points =
      .npcdensbw_effective_cvls_quadrature_points(cvls.quadrature.points, bws$yncon)
  )

  myoptd <- list(
    ftol = 0,
    tol = 0,
    small = 0,
    memfac = 0,
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
    scale.factor.lower.bound = scale.factor.search.lower,
    cvls.quadrature.extend.factor = cvls.quadrature.extend.factor,
    cvls.quadrature.ratios.uniform = cvls.quadrature.ratios[[1L]],
    cvls.quadrature.ratios.sample = cvls.quadrature.ratios[[2L]],
    cvls.quadrature.ratios.gl = cvls.quadrature.ratios[[3L]]
  )

  cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
  cyker.bounds.c <- .npcdensbw_marshal_y_bounds(bws$cykerlb[bws$iycon],
                                                bws$cykerub[bws$iycon],
                                                bws$cykerbound)

  out <- .Call(
    "C_np_density_conditional_bw_eval",
    as.double(yuno),
    as.double(yord),
    as.double(ycon),
    as.double(xuno),
    as.double(xord),
    as.double(xcon),
    as.double(mysd),
    as.integer(myopti),
    as.double(myoptd),
    as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
    as.integer(1L),
    as.integer(penalty_mode),
    as.double(penalty.multiplier),
    as.integer(degree.code),
    as.integer(bernstein.engine),
    as.integer(basis.code),
    as.integer(reg.code),
    as.double(cxker.bounds.c$lb),
    as.double(cxker.bounds.c$ub),
    as.double(cyker.bounds.c$lb),
    as.double(cyker.bounds.c$ub),
    PACKAGE = "npRmpi"
  )

  list(
    objective = as.numeric(out$fval[1L]),
    num.feval = 1L,
    num.feval.fast = as.numeric(as.numeric(out$fast.history[1L]) > 0)
  )
}

npRmpiNomadShadowPrepareConditionalDensity <- function(c.uno,
                                                       c.ord,
                                                       c.con,
                                                       u.uno,
                                                       u.ord,
                                                       u.con,
                                                       mysd,
                                                       myopti,
                                                       myoptd,
                                                       rbw,
                                                       penalty.mode,
                                                       penalty.multiplier,
                                                       degree,
                                                       bernstein,
                                                       basis,
                                                       regtype,
                                                       cxkerlb,
                                                       cxkerub,
                                                       cykerlb,
                                                       cykerub) {
  if (length(myoptd) <= 23L || length(myopti) <= 28L) {
    rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
    if (isTRUE(rank == 0L))
      stop("resident npcdens NOMAD shadow options are missing quadrature controls", call. = FALSE)
    return(FALSE)
  }

  ok <- .Call(
    "C_np_density_conditional_nomad_shadow_prepare",
    c.uno,
    c.ord,
    c.con,
    u.uno,
    u.ord,
    u.con,
    mysd,
    myopti,
    myoptd,
    rbw,
    penalty.mode,
    penalty.multiplier,
    degree,
    bernstein,
    basis,
    regtype,
    cxkerlb,
    cxkerub,
    cykerlb,
    cykerub,
    PACKAGE = "npRmpi"
  )

  if (isTRUE(ok))
    return(TRUE)

  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
  if (isTRUE(rank == 0L))
    stop("failed to prepare resident npcdens NOMAD shadow state", call. = FALSE)

  FALSE
}

npRmpiNomadShadowEvalConditionalDensity <- function(bw, degree) {
  .Call(
    "C_np_density_conditional_nomad_shadow_eval",
    bw,
    degree,
    PACKAGE = "npRmpi"
  )
}

npRmpiNomadShadowClearConditionalDensity <- function() {
  .Call("C_np_density_conditional_nomad_shadow_clear", PACKAGE = "npRmpi")
}

.npcdensbw_nomad_shadow_prepare_args <- function(xdat,
                                                 ydat,
                                                 bws,
                                                 start.bw = NULL,
                                                 invalid.penalty = c("baseline", "dbmax"),
                                                 penalty.multiplier = 10) {
  invalid.penalty <- match.arg(invalid.penalty)

  ydat <- toFrame(ydat)
  xdat <- toFrame(xdat)

  if (length(bws$ybw) != dim(ydat)[2])
    stop("length of bandwidth vector does not match number of columns of 'ydat'")
  if (length(bws$xbw) != dim(xdat)[2])
    stop("length of bandwidth vector does not match number of columns of 'xdat'")
  if (dim(ydat)[1] != dim(xdat)[1])
    stop("number of rows of 'ydat' does not match 'xdat'")

  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE

  xdat <- xdat[keep.rows,, drop = FALSE]
  ydat <- ydat[keep.rows,, drop = FALSE]

  ymat <- toMatrix(ydat)
  xmat <- toMatrix(xdat)

  yuno <- ymat[, bws$iyuno, drop = FALSE]
  ycon <- ymat[, bws$iycon, drop = FALSE]
  yord <- ymat[, bws$iyord, drop = FALSE]
  xuno <- xmat[, bws$ixuno, drop = FALSE]
  xcon <- xmat[, bws$ixcon, drop = FALSE]
  xord <- xmat[, bws$ixord, drop = FALSE]

  mysd <- EssDee(data.frame(xcon, ycon))
  nrow <- nrow(ymat)
  nconfac <- nrow^(-1.0 / (2.0 * bws$cxkerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$cxkerorder + bws$ncon))
  scale.factor.search.lower <- .npcdensbw_resolve_scale_factor_lower_bound(
    npGetScaleFactorSearchLower(bws),
    fallback = 0.1,
    argname = "bws$scale.factor.search.lower"
  )
  cvls.quadrature.extend.factor <- .npcdensbw_resolve_cvls_quadrature_extend_factor(
    bws$cvls.quadrature.extend.factor,
    fallback = 1,
    argname = "bws$cvls.quadrature.extend.factor"
  )
  cvls.quadrature.points <- .npcdensbw_resolve_cvls_quadrature_points(
    bws$cvls.quadrature.points,
    fallback = c(100L, 50L),
    argname = "bws$cvls.quadrature.points"
  )
  cvls.quadrature.ratios <- .npcdensbw_resolve_cvls_quadrature_ratios(
    bws$cvls.quadrature.ratios,
    fallback = c(0.20, 0.55, 0.25),
    argname = "bws$cvls.quadrature.ratios"
  )
  cvls.quadrature.grid <- .npcdensbw_resolve_cvls_quadrature_grid(
    bws$cvls.quadrature.grid,
    fallback = .npcdensbw_cvls_quadrature_grid_fallback(bws$yncon),
    argname = "bws$cvls.quadrature.grid"
  )
  cvls.quadrature.grid <- .npcdensbw_validate_cvls_quadrature_grid_dimension(
    cvls.quadrature.grid,
    bws$yncon,
    argname = "bws$cvls.quadrature.grid"
  )

  penalty_mode <- if (invalid.penalty == "baseline") 1L else 0L
  reg.code <- if (identical(bws$regtype.engine, "lp")) REGTYPE_LP else REGTYPE_LC
  degree.code <- if (bws$xncon > 0L) as.integer(bws$degree.engine) else integer(0L)
  basis.code <- as.integer(npLpBasisCode(bws$basis.engine))
  bernstein.engine <- isTRUE(bws$bernstein.basis.engine)

  myopti <- list(
    num_obs_train = nrow,
    iMultistart = IMULTI_FALSE,
    iNum_Multistart = 0L,
    int_use_starting_values = USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = 0L,
    int_RESTART_FROM_MIN = RE_MIN_FALSE,
    int_MINIMIZE_IO = IO_MIN_TRUE,
    bwmethod = switch(bws$method,
      cv.ml = CBWM_CVML,
      cv.ls = CBWM_CVLS),
    xkerneval = switch(bws$cxkertype,
      gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS),
    ykerneval = switch(bws$cykertype,
      gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS),
    uxkerneval = switch(bws$uxkertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    uykerneval = switch(bws$uykertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    oxkerneval = switch(bws$oxkertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR,
      racineliyan = OKER_RLY),
    oykerneval = switch(bws$oykertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY),
    ynuno = dim(yuno)[2],
    ynord = dim(yord)[2],
    yncon = dim(ycon)[2],
    xnuno = dim(xuno)[2],
    xnord = dim(xord)[2],
    xncon = dim(xcon)[2],
    old.cdens = FALSE,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    scale.init.categorical.sample = FALSE,
    dfc.dir = 0L,
    transform.bounds = FALSE,
    cvls.quadrature.grid = .npcdensbw_cvls_quadrature_grid_code(cvls.quadrature.grid),
    cvls.quadrature.points =
      .npcdensbw_effective_cvls_quadrature_points(cvls.quadrature.points, bws$yncon)
  )

  myoptd <- list(
    ftol = 0,
    tol = 0,
    small = 0,
    memfac = 0,
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
    scale.factor.lower.bound = scale.factor.search.lower,
    cvls.quadrature.extend.factor = cvls.quadrature.extend.factor,
    cvls.quadrature.ratios.uniform = cvls.quadrature.ratios[[1L]],
    cvls.quadrature.ratios.sample = cvls.quadrature.ratios[[2L]],
    cvls.quadrature.ratios.gl = cvls.quadrature.ratios[[3L]]
  )

  cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
  cyker.bounds.c <- .npcdensbw_marshal_y_bounds(bws$cykerlb[bws$iycon],
                                                bws$cykerub[bws$iycon],
                                                bws$cykerbound)
  if (is.null(start.bw)) {
    start.bw <- c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                  bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                  bws$xbw[bws$ixuno], bws$xbw[bws$ixord])
  } else {
    x.offset <- length(bws$ybw)
    start.bw <- c(start.bw[x.offset + which(bws$ixcon)], start.bw[which(bws$iycon)],
                  start.bw[which(bws$iyuno)], start.bw[which(bws$iyord)],
                  start.bw[x.offset + which(bws$ixuno)], start.bw[x.offset + which(bws$ixord)])
  }

  list(
    c.uno = as.double(yuno),
    c.ord = as.double(yord),
    c.con = as.double(ycon),
    u.uno = as.double(xuno),
    u.ord = as.double(xord),
    u.con = as.double(xcon),
    mysd = as.double(mysd),
    myopti = as.integer(myopti),
    myoptd = as.double(myoptd),
    rbw = as.double(start.bw),
    penalty_mode = as.integer(penalty_mode),
    penalty_multiplier = as.double(penalty.multiplier),
    degree = as.integer(degree.code),
    bernstein = as.integer(bernstein.engine),
    basis = as.integer(basis.code),
    regtype = as.integer(reg.code),
    cxkerlb = as.double(cxker.bounds.c$lb),
    cxkerub = as.double(cxker.bounds.c$ub),
    cykerlb = as.double(cyker.bounds.c$lb),
    cykerub = as.double(cyker.bounds.c$ub)
  )
}

npRmpiNomadShadowSearchConditionalDensity <- function(template,
                                                      setup,
                                                      prep,
                                                      degree.search,
                                                      x0,
                                                      bbin,
                                                      lb,
                                                      ub,
                                                      nomad.nmulti = 1L,
                                                      nomad.inner.nmulti = 0L,
                                                      random.seed = 42L,
                                                      use.runtime.bandwidth.progress = FALSE) {
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

  bwdim <- length(setup$cont_flat) + length(setup$cat_flat)
  ndeg <- length(degree.search$start.degree)

  prepared <- npRmpiNomadShadowPrepareConditionalDensity(
    c.uno = prep$c.uno,
    c.ord = prep$c.ord,
    c.con = prep$c.con,
    u.uno = prep$u.uno,
    u.ord = prep$u.ord,
    u.con = prep$u.con,
    mysd = prep$mysd,
    myopti = prep$myopti,
    myoptd = prep$myoptd,
    rbw = prep$rbw,
    penalty.mode = prep$penalty_mode,
    penalty.multiplier = prep$penalty_multiplier,
    degree = prep$degree,
    bernstein = prep$bernstein,
    basis = prep$basis,
    regtype = prep$regtype,
    cxkerlb = prep$cxkerlb,
    cxkerub = prep$cxkerub,
    cykerlb = prep$cykerlb,
    cykerub = prep$cykerub
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
    npRmpiNomadShadowClearConditionalDensity()
  }, add = TRUE)
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[bwdim + seq_len(ndeg)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npcdensbw_nomad_point_to_bw(point[seq_len(bwdim)], template = template, setup = setup)
    x.offset <- length(template$ybw)
    flat.bw <- c(bw_vec[x.offset + which(template$ixcon)], bw_vec[which(template$iycon)],
                 bw_vec[which(template$iyuno)], bw_vec[which(template$iyord)],
                 bw_vec[x.offset + which(template$ixuno)], bw_vec[x.offset + which(template$ixord)])
    out <- npRmpiNomadShadowEvalConditionalDensity(
      bw = as.double(flat.bw),
      degree = as.integer(degree)
    )
    nomad.num.feval.total <<- nomad.num.feval.total + as.numeric(out[2L])
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out[3L])

    list(
      objective = as.numeric(out[1L]),
      degree = degree,
      num.feval = as.numeric(out[2L]),
      num.feval.fast = as.numeric(out[3L])
    )
  }

  external.progress <- if (isTRUE(use.runtime.bandwidth.progress) &&
                           isTRUE(rank == 0L) &&
                           !is.null(.np_progress_runtime$bandwidth_state)) {
    .np_progress_runtime$bandwidth_state
  } else {
    NULL
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
    direction = "max",
    objective_name = "fval",
    nmulti = nomad.nmulti,
    nomad.inner.nmulti = nomad.inner.nmulti,
    random.seed = random.seed,
    progress_state = external.progress,
    manage_progress_lifecycle = is.null(external.progress),
    bind_bandwidth_runtime = !is.null(external.progress),
    handoff_before_build = identical(degree.search$engine, "nomad+powell"),
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

.npcdensbw_nomad_bw_setup <- function(xdat,
                                      ydat,
                                      template,
                                      bandwidth.scale.categorical = 1e4) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  xmat <- toMatrix(xdat)
  ymat <- toMatrix(ydat)
  xcon <- xmat[, template$ixcon, drop = FALSE]
  ycon <- ymat[, template$iycon, drop = FALSE]
  nrow <- nrow(xmat)
  nconfac <- nrow^(-1.0 / (2.0 * template$cxkerorder + template$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * template$cxkerorder + template$ncon))

  x_offset <- length(template$ybw)
  y_cont_flat <- which(template$iycon)
  x_cont_flat <- x_offset + which(template$ixcon)
  y_uno_flat <- which(template$iyuno)
  y_ord_flat <- which(template$iyord)
  x_uno_flat <- x_offset + which(template$ixuno)
  x_ord_flat <- x_offset + which(template$ixord)

  cat_upper_one <- function(values, kernel) {
    if (identical(kernel, "aitchisonaitken")) {
      nlev <- length(unique(values))
      return((nlev - 1) / nlev)
    }
    1
  }

  cat_upper <- c(
    if (length(y_uno_flat)) vapply(which(template$iyuno), function(i) cat_upper_one(ydat[[i]], template$uykertype), numeric(1L)) else numeric(0L),
    rep.int(1, length(y_ord_flat)),
    if (length(x_uno_flat)) vapply(which(template$ixuno), function(i) cat_upper_one(xdat[[i]], template$uxkertype), numeric(1L)) else numeric(0L),
    rep.int(1, length(x_ord_flat))
  )

  list(
    cont_flat = c(y_cont_flat, x_cont_flat),
    cont_scale = c(EssDee(ycon), EssDee(xcon)) * nconfac,
    cat_flat = c(y_uno_flat, y_ord_flat, x_uno_flat, x_ord_flat),
    ncatfac = ncatfac,
    bandwidth.scale.categorical = bandwidth.scale.categorical,
    cat_upper = cat_upper
  )
}

.npcdensbw_nomad_point_to_bw <- function(point, template, setup) {
  point <- as.numeric(point)
  ncont <- length(setup$cont_flat)
  ncat <- length(setup$cat_flat)
  bws <- numeric(length(template$ybw) + length(template$xbw))

  if (ncont > 0L) {
    gamma <- point[seq_len(ncont)]
    ext_bw <- gamma * setup$cont_scale
    bws[setup$cont_flat] <- if (isTRUE(template$scaling)) gamma else ext_bw
  }

  if (ncat > 0L) {
    lambda_scaled <- point[ncont + seq_len(ncat)]
    ext_bw <- lambda_scaled / setup$bandwidth.scale.categorical
    bws[setup$cat_flat] <- if (isTRUE(template$scaling)) ext_bw / setup$ncatfac else ext_bw
  }

  bws
}

.npcdensbw_nomad_bw_to_point <- function(bws, template, setup) {
  point <- numeric(length(setup$cont_flat) + length(setup$cat_flat))

  if (length(setup$cont_flat) > 0L) {
    raw <- bws[setup$cont_flat]
    point[seq_along(setup$cont_flat)] <- if (isTRUE(template$scaling)) {
      raw
    } else {
      raw / setup$cont_scale
    }
  }

  if (length(setup$cat_flat) > 0L) {
    raw <- bws[setup$cat_flat]
    ext_bw <- if (isTRUE(template$scaling)) raw * setup$ncatfac else raw
    point[length(setup$cont_flat) + seq_along(setup$cat_flat)] <- ext_bw * setup$bandwidth.scale.categorical
  }

  point
}

.npcdensbw_nomad_continuous_lower_bound <- function(template) {
  npGetScaleFactorSearchLower(
    template,
    fallback = 0.1,
    argname = "template$scale.factor.search.lower"
  )
}

.npcdensbw_powell_progress_fields <- function(state,
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

.npcdensbw_with_powell_refinement_progress <- function(degree, expr) {
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
    active.state$unknown_total_fields <- .npcdensbw_powell_progress_fields
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

.npcdensbw_nomad_search <- function(xdat,
                                    ydat,
                                    bws,
                                    reg.args,
                                    opt.args,
                                    degree.search,
                                    nomad.inner.nmulti = 0L,
                                    random.seed = 42L) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$pregtype <- "Local-Polynomial"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template.reg.args$regtype.engine <- "lp"
  template.reg.args$degree.engine <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

  template <- .npcdensbw_build_conbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args
  )

  if (!identical(template$type, "fixed"))
    stop("automatic degree search with search.engine='nomad' currently requires bwtype='fixed'")
  setup <- .npcdensbw_nomad_bw_setup(xdat = xdat, ydat = ydat, template = template)
  bwdim <- length(setup$cont_flat) + length(setup$cat_flat)
  ndeg <- length(degree.search$start.degree)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(dim(ydat)[2]+dim(xdat)[2]) else npValidateNmulti(opt.args$nmulti[1L])
  objective.direction <- "max"
  cont_lower <- .npcdensbw_nomad_continuous_lower_bound(template)
  bw_lower <- c(rep.int(cont_lower, length(setup$cont_flat)), rep.int(0, length(setup$cat_flat)))
  bw_upper <- c(rep.int(1e6, length(setup$cont_flat)), setup$cat_upper * setup$bandwidth.scale.categorical)

  x0 <- c(
    .np_nomad_complete_start_point(
      point = {
        raw <- c(template$ybw, template$xbw)
        if (all(raw == 0)) NULL else .npcdensbw_nomad_bw_to_point(raw, template = template, setup = setup)
      },
      lower = bw_lower,
      upper = bw_upper,
      ncont = length(setup$cont_flat)
    ),
    as.integer(degree.search$start.degree)
  )
  lb <- c(bw_lower, degree.search$lower)
  ub <- c(bw_upper, degree.search$upper)
  bbin <- c(rep.int(0L, bwdim), rep.int(1L, ndeg))
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0
  direct.meta.context <- if (is.recursive(bws) &&
                             !is.null(bws$nconfac) &&
                             !is.null(bws$ncatfac) &&
                             !is.null(bws$sdev) &&
                             !identical(bws$nconfac, NA) &&
                             !identical(bws$ncatfac, NA) &&
                             !identical(bws$sdev, NA)) {
    list(
      nconfac = bws$nconfac,
      ncatfac = bws$ncatfac,
      sdev = bws$sdev
    )
  } else {
    txmat <- toMatrix(xdat)
    tymat <- toMatrix(ydat)
    xcon <- txmat[, template$ixcon, drop = FALSE]
    ycon <- tymat[, template$iycon, drop = FALSE]
    list(
      nconfac = nrow(xdat)^(-1.0 / (2.0 * template$cxkerorder + template$ncon)),
      ncatfac = nrow(xdat)^(-2.0 / (2.0 * template$cxkerorder + template$ncon)),
      sdev = EssDee(data.frame(xcon, ycon))
    )
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    degree <- as.integer(best_record$degree)
    bw_vec <- .npcdensbw_nomad_point_to_bw(point[seq_len(bwdim)], template = template, setup = setup)
    powell.elapsed <- NA_real_

    build_direct_payload <- function() {
      final.reg.args <- reg.args
      final.reg.args$regtype <- "lp"
      final.reg.args$pregtype <- "Local-Polynomial"
      final.reg.args$degree <- degree
      final.reg.args$bernstein.basis <- degree.search$bernstein.basis
      final.reg.args$regtype.engine <- "lp"
      final.reg.args$degree.engine <- degree
      final.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

      tbw <- .npcdensbw_build_conbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bw_vec,
        bandwidth.compute = FALSE,
        reg.args = final.reg.args
      )
      payload <- tbw
      payload$nconfac <- direct.meta.context$nconfac
      payload$ncatfac <- direct.meta.context$ncatfac
      payload$sdev <- direct.meta.context$sdev
      payload <- .np_refresh_xy_bandwidth_metadata(payload)
      payload$method <- if (!is.null(payload$method) && length(payload$method)) {
        as.character(payload$method[1L])
      } else if (!is.null(reg.args$bwmethod) && length(reg.args$bwmethod)) {
        as.character(reg.args$bwmethod[1L])
      } else {
        "cv.ml"
      }
      payload$pmethod <- bwmToPrint(payload$method)
      payload$fval <- as.numeric(best_record$objective)
      payload$ifval <- NA_real_
      payload$num.feval <- as.numeric(nomad.num.feval.total)
      payload$num.feval.fast <- as.numeric(nomad.num.feval.fast.total)
      payload$fval.history <- NA_real_
      payload$eval.history <- NA_real_
      payload$invalid.history <- NA_real_
      payload$timing <- NA_real_
      payload$total.time <- NA_real_
      payload
    }

    direct.payload <- build_direct_payload()
    if (is.null(direct.payload$timing.profile) && is.list(best_record$timing.profile))
      direct.payload$timing.profile <- best_record$timing.profile
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$pregtype <- "Local-Polynomial"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.reg.args$regtype.engine <- "lp"
      hot.reg.args$degree.engine <- degree
      hot.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
      hot.opt.args <- opt.args
      hot.opt.args$nmulti <- .np_nomad_powell_hotstart_nmulti("disable_multistart")
      powell.start <- proc.time()[3L]
      hot.payload <- .npcdensbw_with_powell_refinement_progress(
        degree,
        .npcdensbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          bws = bw_vec,
          reg.args = hot.reg.args,
          opt.args = hot.opt.args
        )
      )
      powell.elapsed <- proc.time()[3L] - powell.start
      direct.payload$num.feval <- as.numeric(direct.payload$num.feval[1L]) + as.numeric(hot.payload$num.feval[1L])
      direct.payload$num.feval.fast <- as.numeric(direct.payload$num.feval.fast[1L]) + as.numeric(hot.payload$num.feval.fast[1L])
      hot.payload$num.feval <- direct.payload$num.feval
      hot.payload$num.feval.fast <- direct.payload$num.feval.fast
      if (!is.null(hot.payload$method) && length(hot.payload$method))
        hot.payload$pmethod <- bwmToPrint(as.character(hot.payload$method[1L]))
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = objective.direction)) {
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  if (.npRmpi_has_active_slave_pool(comm = 1L) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    start.bw <- .npcdensbw_nomad_point_to_bw(x0[seq_len(bwdim)], template = template, setup = setup)
    prep <- .npcdensbw_nomad_shadow_prepare_args(
      xdat = xdat,
      ydat = ydat,
      bws = template,
      start.bw = start.bw,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier
    )
    # Workers only need the fields consumed by the shadow NOMAD search.
    search.template <- list(
      scaling = isTRUE(template$scaling),
      ybw = template$ybw,
      xbw = template$xbw,
      iycon = template$iycon,
      ixcon = template$ixcon,
      iyuno = template$iyuno,
      iyord = template$iyord,
      ixuno = template$ixuno,
      ixord = template$ixord
    )
    search.setup <- list(
      cont_flat = setup$cont_flat,
      cont_scale = setup$cont_scale,
      cat_flat = setup$cat_flat,
      ncatfac = setup$ncatfac,
      bandwidth.scale.categorical = setup$bandwidth.scale.categorical
    )
    search.degree <- list(
      engine = degree.search$engine,
      start.degree = degree.search$start.degree,
      candidates = degree.search$candidates,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      start.user = degree.search$start.user
    )

    mc <- substitute(
      get("npRmpiNomadShadowSearchConditionalDensity", envir = asNamespace("npRmpi"), inherits = FALSE)(
        TEMPLATE,
        SETUP,
        PREP,
        DEGREESEARCH,
        X0,
        BBIN,
        LB,
        UB,
        NOMADNMULTI,
        INNERNMULTI,
        RSEED,
        RPROGRESS
      ),
      list(
        TEMPLATE = search.template,
        SETUP = search.setup,
        PREP = prep,
        DEGREESEARCH = search.degree,
        X0 = x0,
        BBIN = bbin,
        LB = lb,
        UB = ub,
        NOMADNMULTI = nomad.nmulti,
        INNERNMULTI = nomad.inner.nmulti,
        RSEED = random.seed,
        RPROGRESS = TRUE
      )
    )

    if (isTRUE(.npRmpi_autodispatch_called_from_bcast())) {
      search.result <- eval(mc, envir = environment())
    } else {
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
          .np_degree_better(payload_result$objective, search.result$best$objective, direction = objective.direction)) {
        search.result$best$objective <- as.numeric(payload_result$objective[1L])
      }
    } else {
      search.result$best_payload <- payload_result
    }

    return(.npRmpi_reconcile_nomad_search_timing(search.result))
  }

  .np_nomad_baseline_note(degree.search$start.degree)

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[bwdim + seq_len(ndeg)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npcdensbw_nomad_point_to_bw(point[seq_len(bwdim)], template = template, setup = setup)

    eval.reg.args <- reg.args
    eval.reg.args$regtype <- "lp"
    eval.reg.args$pregtype <- "Local-Polynomial"
    eval.reg.args$degree <- degree
    eval.reg.args$bernstein.basis <- degree.search$bernstein.basis
    eval.reg.args$regtype.engine <- "lp"
    eval.reg.args$degree.engine <- degree
    eval.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

    tbw <- .npcdensbw_build_conbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = eval.reg.args
    )

    out <- .npcdensbw_eval_only(
      xdat = xdat,
      ydat = ydat,
      bws = tbw,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier
    )
    nomad.num.feval.total <<- nomad.num.feval.total + as.numeric(out$num.feval[1L])
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])

    list(
      objective = out$objective,
      degree = degree,
      num.feval = out$num.feval,
      num.feval.fast = out$num.feval.fast
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
    direction = objective.direction,
    objective_name = "fval",
    nmulti = nomad.nmulti,
    nomad.inner.nmulti = nomad.inner.nmulti,
    random.seed = random.seed,
    handoff_before_build = identical(degree.search$engine, "nomad+powell"),
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

.npcdensbw_degree_search_controls <- function(regtype,
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
  search.engine <- .np_degree_search_engine_controls(search.engine)

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (ncon < 1L)
    stop("automatic degree search requires at least one continuous conditioning predictor")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = ncon,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  baseline.degree <- rep.int(0L, ncon)
  # Density/distribution degree search should anchor on the local-constant
  # baseline unless the user explicitly requests another start.
  default.start.degree <- baseline.degree
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, default.start.degree))
  } else {
    start.raw <- npValidateGlpDegree(regtype = "lp", degree = degree.start, ncon = ncon, argname = "degree.start")
    out.of.range <- vapply(seq_len(ncon), function(j) !(start.raw[j] %in% bounds$candidates[[j]]), logical(1))
    if (any(out.of.range))
      stop("degree.start must lie within the searched degree candidates for every continuous conditioning predictor")
    start.raw
  }

  list(
    method = if (identical(search.engine, "cell")) degree.select else search.engine,
    engine = search.engine,
    candidates = bounds$candidates,
    lower = bounds$lower,
    upper = bounds$upper,
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

.npcdensbw_attach_degree_search <- function(bws, search_result) {
  metadata <- list(
    mode = search_result$method,
    direction = search_result$direction,
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

npcdensbw.NULL <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           bws, ...){
    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain y names and 'toFrame'
    ydat <- toFrame(ydat)

    ## do bandwidths

    bws = double(ncol(ydat)+ncol(xdat))

    tbw <- npcdensbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw
  }

npcdensbw.default <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           bws,
           bandwidth.compute = TRUE,
           bwmethod,
           bwscaling,
           bwtype,
           cfac.dir,
           scale.factor.init,
           cxkerbound,
           cxkerlb,
           cxkerorder,
           cxkertype,
           cxkerub,
           cykerbound,
           cykerlb,
           cykerorder,
           cykertype,
           cykerub,
           dfac.dir,
           dfac.init,
           dfc.dir,
           ftol,
           scale.factor.init.upper,
           hbd.dir,
           hbd.init,
           initc.dir,
           initd.dir,
           invalid.penalty,
           itmax,
           lbc.dir,
           scale.factor.init.lower,
           lbd.dir,
           lbd.init,
           memfac,
           nmulti,
           oxkertype,
           oykertype,
           penalty.multiplier,
           remin,
           scale.init.categorical.sample,
           scale.factor.search.lower = NULL,
           cvls.quadrature.grid = c("hybrid", "uniform", "sample"),
           cvls.quadrature.extend.factor = 1,
           cvls.quadrature.points = c(100L, 50L),
           cvls.quadrature.ratios = c(0.20, 0.55, 0.25),
           small,
           tol,
           transform.bounds,
           uxkertype,
           uykertype,
           regtype = c("lc", "ll", "lp"),
           basis = c("glp", "additive", "tensor"),
           degree = NULL,
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
           bernstein.basis = FALSE,
           ## dummy arguments for conbandwidth() function call
           ...){
    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain y names and 'toFrame'
    ydat <- toFrame(ydat)

    x.info <- untangle(xdat)
    y.info <- untangle(ydat)

    mc <- match.call(expand.dots = FALSE)
    mc.names <- names(mc)
    mc.expanded <- match.call(expand.dots = TRUE)
    if ("cvls.i1.rescue" %in% names(mc.expanded))
      stop("cvls.i1.rescue has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive.tol" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive.tol has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive.grid.hy.ratio" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive.grid.hy.ratio has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    if ("cvls.quadrature.adaptive.floor.tol" %in% names(mc.expanded))
      stop("cvls.quadrature.adaptive.floor.tol has been removed; use cvls.quadrature.grid",
           call. = FALSE)
    cvls.quadrature.grid <- .npcdensbw_resolve_cvls_quadrature_grid(
      if ("cvls.quadrature.grid" %in% mc.names) cvls.quadrature.grid else NULL,
      fallback = .npcdensbw_cvls_quadrature_grid_fallback(sum(y.info$icon)),
      argname = "cvls.quadrature.grid"
    )
    cvls.quadrature.grid <- .npcdensbw_validate_cvls_quadrature_grid_dimension(
      cvls.quadrature.grid,
      sum(y.info$icon),
      argname = "cvls.quadrature.grid"
    )
    scale.factor.search.lower <- .npcdensbw_resolve_scale_factor_lower_bound(
      scale.factor.search.lower,
      fallback = 0.1,
      argname = "scale.factor.search.lower"
    )
    cvls.quadrature.extend.factor <- .npcdensbw_resolve_cvls_quadrature_extend_factor(
      cvls.quadrature.extend.factor,
      fallback = 1,
      argname = "cvls.quadrature.extend.factor"
    )
    cvls.quadrature.points <- .npcdensbw_resolve_cvls_quadrature_points(
      cvls.quadrature.points,
      fallback = c(100L, 50L),
      argname = "cvls.quadrature.points"
    )
    cvls.quadrature.ratios <- .npcdensbw_resolve_cvls_quadrature_ratios(
      cvls.quadrature.ratios,
      fallback = c(0.20, 0.55, 0.25),
      argname = "cvls.quadrature.ratios"
    )
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
      where = "npcdensbw"
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

    regtype.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "regtype")
    basis.named <- any(mc.names == "basis")
    degree.named <- any(mc.names == "degree")
    bernstein.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "bernstein.basis")

    regtype <- if (!is.null(nomad.shortcut$values$regtype)) {
      match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp"))
    } else {
      "lc"
    }
    if (identical(regtype, "lc") && (basis.named || degree.named || bernstein.named))
      stop("regtype='lc' does not accept basis/degree/bernstein.basis; use regtype='lp' for local-polynomial controls")
    if (identical(regtype, "ll")) {
      if (degree.named)
        stop("regtype='ll' uses canonical LP(degree=1, basis='glp'); remove 'degree' or use regtype='lp'")
      if (basis.named && !identical(match.arg(basis), "glp"))
        stop("regtype='ll' uses canonical basis='glp'; use regtype='lp' for alternate LP bases")
      if (bernstein.named && isTRUE(bernstein.basis))
        stop("regtype='ll' uses canonical bernstein.basis=FALSE; use regtype='lp' for Bernstein LP")
    }

    bernstein.value <- if (!is.null(nomad.shortcut$values$bernstein.basis)) {
      nomad.shortcut$values$bernstein.basis
    } else {
      bernstein.basis
    }
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    degree.setup <- npSetupGlpDegree(
      regtype = regtype,
      degree = degree,
      ncon = sum(x.info$icon),
      degree.select = degree.select.value
    )
    spec <- npCanonicalConditionalRegSpec(
      regtype = regtype,
      basis = basis,
      degree = degree.setup,
      bernstein.basis = bernstein.value,
      ncon = sum(x.info$icon),
      where = "npcdensbw"
    )
    pregtype <- switch(spec$regtype,
                       lc = "Local-Constant",
                       ll = "Local-Linear",
                       lp = "Local-Polynomial")

    search.mc.names <- names(mc)
    lp.dot.args <- list(...)
    random.seed.value <- .np_degree_extract_random_seed(lp.dot.args)
    search.engine.value <- if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell"
    degree.min.value <- nomad.shortcut$values$degree.min
    degree.max.value <- nomad.shortcut$values$degree.max
    degree.start.value <- if ("degree.start" %in% search.mc.names) degree.start else NULL
    degree.restarts.value <- if ("degree.restarts" %in% search.mc.names) degree.restarts else 0L
    degree.max.cycles.value <- if ("degree.max.cycles" %in% search.mc.names) degree.max.cycles else 20L
    degree.verify.value <- if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE
    degree.search <- .npcdensbw_degree_search_controls(
      regtype = regtype,
      regtype.named = regtype.named,
      ncon = sum(x.info$icon),
      nobs = NROW(xdat),
      basis = if (basis.named) basis else "glp",
      degree.select = degree.select.value,
      search.engine = search.engine.value,
      degree.min = degree.min.value,
      degree.max = degree.max.value,
      degree.start = degree.start.value,
      degree.restarts = degree.restarts.value,
      degree.max.cycles = degree.max.cycles.value,
      degree.verify = degree.verify.value,
      bernstein.basis = bernstein.value,
      bernstein.named = bernstein.named
    )
    nomad.inner <- .np_nomad_validate_inner_multistart(
      call_names = mc.names,
      dot.args = lp.dot.args,
      nomad.nmulti = nomad.nmulti,
      regtype = regtype,
      automatic.degree.search = !is.null(degree.search),
      search.engine = if (is.null(degree.search)) "" else degree.search$engine
    )
    nomad.inner.named <- nomad.inner$named
    nomad.inner.nmulti <- nomad.inner$nmulti

    if (!is.null(degree.search)) {
      spec$bernstein.basis <- degree.search$bernstein.basis
      spec$bernstein.basis.engine <- degree.search$bernstein.basis
    }

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    mc.names <- names(mc)
    margs <- c("bwmethod", "bwscaling", "bwtype", "cxkertype", "cxkerorder",
               "cxkerbound", "cxkerlb", "cxkerub",
               "cykertype", "cykerorder", "cykerbound", "cykerlb", "cykerub",
               "uxkertype", "uykertype", "oxkertype", "oykertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    y.idx <- seq_len(length(ydat))
    x.idx <- seq_len(length(xdat))
    bw.args <- list(
      xbw = bws[length(ydat) + x.idx],
      ybw = bws[y.idx],
      nobs = nrow(xdat),
      xdati = x.info,
      ydati = y.info,
      xnames = names(xdat),
      ynames = names(ydat),
      bandwidth.compute = bandwidth.compute,
      regtype = spec$regtype,
      pregtype = pregtype,
      basis = spec$basis,
      degree = spec$degree,
      bernstein.basis = spec$bernstein.basis,
      regtype.engine = spec$regtype.engine,
      basis.engine = spec$basis.engine,
      degree.engine = spec$degree.engine,
      bernstein.basis.engine = spec$bernstein.basis.engine
    )
    if (any.m) {
      nms <- mc.names[m]
      bw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    reg.args <- bw.args[setdiff(names(bw.args), c("xbw", "ybw", "nobs", "xdati", "ydati", "xnames", "ynames", "bandwidth.compute"))]
    reg.args$scale.factor.search.lower <- scale.factor.search.lower
    reg.args$cvls.quadrature.grid <- cvls.quadrature.grid
    reg.args$cvls.quadrature.extend.factor <- cvls.quadrature.extend.factor
    reg.args$cvls.quadrature.points <- cvls.quadrature.points
    reg.args$cvls.quadrature.ratios <- cvls.quadrature.ratios
    reg.bwmethod <- if (is.null(reg.args$bwmethod)) "cv.ls" else reg.args$bwmethod
    if (isTRUE(bandwidth.compute) &&
        identical(as.character(reg.bwmethod)[1L], "cv.ls")) {
      .npcdensbw_warn_infinite_response_quadrature(
        reg.args$cykerlb,
        reg.args$cykerub,
        reg.args$cykerbound,
        points.supplied = "cvls.quadrature.points" %in% mc.names
      )
    }
    tbw <- do.call(conbandwidth, bw.args)
    .npRmpi_require_active_slave_pool(where = "npcdensbw()")
    keep_local_shadow_nn <- bandwidth.compute &&
      identical(tbw$regtype.engine, "lp") &&
      identical(tbw$method %in% c("cv.ml", "cv.ls"), TRUE) &&
      identical(tbw$type %in% c("generalized_nn", "adaptive_nn"), TRUE)
    keep_local_raw_degree1_cvls <- bandwidth.compute &&
      identical(tbw$method, "cv.ls") &&
      identical(tbw$type, "fixed") &&
      npIsRawDegreeOneConditionalSpec(spec, tbw$xncon)
    if (.npRmpi_autodispatch_active() &&
        !keep_local_shadow_nn &&
        !keep_local_raw_degree1_cvls &&
        is.null(degree.search))
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))
    ## next grab dummies for actual bandwidth selection and perform call

    mc.names <- names(mc)
    margs <- c("nmulti", "remin", "itmax", "ftol",
               "tol", "small", "memfac",
               "lbc.dir", "dfc.dir", "cfac.dir","initc.dir",
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.init.categorical.sample",
               "transform.bounds",
               "cvls.quadrature.grid",
               "cvls.quadrature.extend.factor",
               "cvls.quadrature.points",
               "cvls.quadrature.ratios",
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
    opt.args <- c(list(bandwidth.compute = bandwidth.compute), opt.args)

    if (!is.null(degree.search)) {
      if (identical(degree.search$engine, "cell")) {
        eval_fun <- function(degree.vec) {
          cell.reg.args <- reg.args
          cell.reg.args$regtype <- "lp"
          cell.reg.args$pregtype <- "Local-Polynomial"
          cell.reg.args$degree <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
          cell.reg.args$regtype.engine <- "lp"
          cell.reg.args$degree.engine <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
          cell.bws <- .npcdensbw_run_fixed_degree(
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            reg.args = cell.reg.args,
            opt.args = opt.args
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
          direction = "max",
          trace_level = "full",
          objective_name = "fval"
        )
      } else {
        search.result <- .npcdensbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          reg.args = reg.args,
          opt.args = opt.args,
          degree.search = degree.search,
          nomad.inner.nmulti = nomad.inner.nmulti,
          random.seed = random.seed.value
        )
      }
      tbw <- .npcdensbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
    } else {
      tbw <- .npcdensbw_build_conbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bws,
        bandwidth.compute = bandwidth.compute,
        reg.args = reg.args
      )
      bwsel.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting conditional density bandwidth",
        do.call(npcdensbw.conbandwidth, bwsel.args)
      )
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)
  }
