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
    
    tbw <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat), list(...)))

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
    dots <- list(...)
    .np_nomad_native_reject_unsupported_options_from_dots(
      dots,
      "native npreg NOMAD route"
    )

    xdat <- toFrame(xdat)

    bws = double(dim(xdat)[2])
    
    tbw <- do.call(npregbw.default,
                   c(list(xdat = xdat, ydat = ydat, bws = bws), dots))

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)
    
    tbw
  }

.npregbw_nomad_native_target <- function(template, bwsolver) {
  method <- if (!is.null(template$method) && length(template$method)) {
    as.character(template$method[1L])
  } else {
    "cv.ls"
  }
  bwtype <- if (!is.null(template$type) && length(template$type)) {
    as.character(template$type[1L])
  } else {
    ""
  }
  regtype <- if (!is.null(template$regtype) && length(template$regtype)) {
    as.character(template$regtype[1L])
  } else {
    "lc"
  }

  bwtype %in% c("fixed", "generalized_nn", "adaptive_nn") &&
    method %in% c("cv.ls", "cv.aic") &&
    regtype %in% c("lc", "ll", "lp") &&
    bwsolver %in% c("mads", "mads+powell")
}

.npregbw_nomad_degree_native_target <- function(template, degree.search) {
  method <- if (!is.null(template$method) && length(template$method)) {
    as.character(template$method[1L])
  } else {
    "cv.ls"
  }
  bwtype <- if (!is.null(template$type) && length(template$type)) {
    as.character(template$type[1L])
  } else {
    ""
  }
  regtype <- if (!is.null(template$regtype) && length(template$regtype)) {
    as.character(template$regtype[1L])
  } else {
    ""
  }
  engine <- if (!is.null(degree.search$engine) && length(degree.search$engine)) {
    as.character(degree.search$engine[1L])
  } else {
    ""
  }

  bwtype %in% c("fixed", "generalized_nn", "adaptive_nn") &&
    method %in% c("cv.ls", "cv.aic") &&
    identical(regtype, "lp") &&
    engine %in% c("nomad", "nomad+powell")
}

.npregbw_nomad_native_require_crs <- function() {
  if (!requireNamespace("crs", quietly = TRUE))
    stop("native npreg NOMAD route requires crs >= 0.15-44", call. = FALSE)
  if (utils::packageVersion("crs") < "0.15.44")
    stop("native npreg NOMAD route requires crs >= 0.15-44", call. = FALSE)
  invisible(TRUE)
}

.npregbw_nomad_native_option_vectors <- function(opts) {
  if (is.null(opts) || !length(opts))
    return(list(names = character(), values = character()))

  .np_nomad_native_reject_unsupported_options(opts, "native npreg NOMAD route")

  option.names <- names(opts)
  if (is.null(option.names) || any(!nzchar(option.names)))
    stop("native npreg NOMAD route received unnamed NOMAD options", call. = FALSE)

  option.values <- vapply(opts, function(value) {
    if (is.logical(value)) {
      if (isTRUE(value[1L])) "true" else "false"
    } else if (length(value) > 1L) {
      paste0("( ", paste(as.character(value), collapse = " "), " )")
    } else {
      as.character(value[1L])
    }
  }, character(1L))

  list(names = as.character(option.names), values = option.values)
}

.npregbw_nomad_native_prepare_args <- function(xdat,
                                               ydat,
                                               bws,
                                               invalid.penalty = c("baseline", "dbmax"),
                                               penalty.multiplier = 10,
                                               itmax = 10000L,
                                               ftol = 1.490116e-07,
                                               tol = 1.490116e-04,
                                               small = 1.490116e-05,
                                               lbc.dir = 0.5,
                                               cfac.dir = 2.5 * (3.0 - sqrt(5)),
                                               initc.dir = 1.0,
                                               lbd.dir = 0.1,
                                               hbd.dir = 1,
                                               dfac.dir = 0.25 * (3.0 - sqrt(5)),
                                               initd.dir = 1.0,
                                               scale.factor.init.lower = 0.1,
                                               scale.factor.init.upper = 2.0,
                                               scale.factor.init = 0.5,
                                               lbd.init = 0.1,
                                               hbd.init = 0.9,
                                               dfac.init = 0.375,
                                               scale.factor.search.lower = NULL,
                                               scale.init.categorical.sample = FALSE,
                                               transform.bounds = FALSE) {
  invalid.penalty <- match.arg(invalid.penalty)
  scale.factor.search.lower <- npResolveScaleFactorLowerBound(
    if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
  )
  cont.start <- npContinuousSearchStartControls(
    scale.factor.init.lower,
    scale.factor.init.upper,
    scale.factor.init,
    scale.factor.search.lower,
    where = "npregbw"
  )

  xdat <- toFrame(xdat)
  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector")

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
  cont_scale <- mysd * nconfac
  bandwidth.scale.categorical <- 1e4

  reg.spec <- npCanonicalConditionalRegSpec(
    regtype = bws$regtype,
    basis = bws$basis,
    degree = bws$degree,
    bernstein.basis = bws$bernstein.basis,
    ncon = bws$ncon,
    where = "npregbw"
  )
  reg.c <- npRegtypeToC(regtype = reg.spec$regtype.engine,
                        degree = reg.spec$degree.engine,
                        ncon = bws$ncon,
                        context = "npregbw")
  degree.c <- if (bws$ncon > 0) as.integer(reg.spec$degree.engine) else integer(1L)

  myopti <- list(
    num_obs_train = nrow,
    iMultistart = IMULTI_TRUE,
    iNum_Multistart = 1L,
    int_use_starting_values = USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_reg_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = itmax,
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
    int_do_tree = .npregbw_tree_code(bws, ncon = bws$ncon, ncat = bws$nuno + bws$nord),
    scale.init.categorical.sample = scale.init.categorical.sample,
    dfc.dir = 3L,
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
    lbc.init = cont.start$scale.factor.init.lower,
    hbc.init = cont.start$scale.factor.init.upper,
    cfac.init = cont.start$scale.factor.init,
    lbd.init = lbd.init,
    hbd.init = hbd.init,
    dfac.init = dfac.init,
    nconfac = nconfac,
    ncatfac = ncatfac,
    scale.factor.lower.bound = scale.factor.search.lower
  )

  cker.bounds <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])
  decode.scale <- numeric(length(bws$bw))
  if (length(bws$icon) > 0L)
    decode.scale[bws$icon] <- if (identical(bws$type, "fixed") && !isTRUE(bws$scaling)) cont_scale else 1
  if (length(bws$iuno) > 0L && any(bws$iuno)) {
    decode.scale[bws$iuno] <- if (isTRUE(bws$scaling)) {
      1 / (bandwidth.scale.categorical * ncatfac)
    } else {
      1 / bandwidth.scale.categorical
    }
  }
  if (length(bws$iord) > 0L && any(bws$iord)) {
    decode.scale[bws$iord] <- if (isTRUE(bws$scaling)) {
      1 / (bandwidth.scale.categorical * ncatfac)
    } else {
      1 / bandwidth.scale.categorical
    }
  }
  decode.scale <- as.double(c(decode.scale[bws$icon], decode.scale[bws$iuno], decode.scale[bws$iord]))
  list(
    runo = as.double(runo),
    rord = as.double(rord),
    rcon = as.double(rcon),
    y = as.double(ydat),
    mysd = as.double(mysd),
    myopti = as.integer(myopti),
    myoptd = as.double(myoptd),
    degree = as.integer(degree.c),
    bernstein = as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
    basis = as.integer(npLpBasisCode(reg.spec$basis.engine)),
    penalty_mode = as.integer(if (invalid.penalty == "baseline") 1L else 0L),
    penalty_multiplier = as.double(penalty.multiplier),
    ckerlb = as.double(cker.bounds$lb),
    ckerub = as.double(cker.bounds$ub),
    decode_scale = decode.scale
  )
}

npNomadNativeSearchRegression <- function(prep,
                                          x0,
                                          bbin,
                                          lb,
                                          ub,
                                          max.eval = 0L,
                                          random.seed = 42L,
                                          inner.start.count = 0L,
                                          option.names = character(),
                                          option.values = character()) {
  native.call <- .np_nomad_capture_solver_output(.Call(
    "C_np_regression_nomad_native_search",
    as.double(prep$runo),
    as.double(prep$rord),
    as.double(prep$rcon),
    as.double(prep$y),
    as.double(prep$mysd),
    as.integer(prep$myopti),
    as.double(prep$myoptd),
    as.integer(prep$degree),
    as.integer(prep$bernstein),
    as.integer(prep$basis),
    as.double(x0),
    as.integer(bbin),
    as.double(lb),
    as.double(ub),
    as.integer(max.eval),
    as.integer(random.seed),
    as.integer(inner.start.count),
    as.character(option.names),
    as.character(option.values),
    as.integer(prep$penalty_mode),
    as.double(prep$penalty_multiplier),
    as.double(prep$ckerlb),
    as.double(prep$ckerub),
    as.double(prep$decode_scale),
    PACKAGE = "np"
  ), capture.output = TRUE)
  .np_nomad_native_call_value(native.call)
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
           powell.remin = TRUE,
           bwsolver = c("powell", "mads", "mads+powell"),
           scale.init.categorical.sample = FALSE,
           scale.factor.search.lower = NULL,
           small = 1.490116e-05,
           tol = 1.490116e-04,
           transform.bounds = FALSE,
           ...,
           nomad.opts = list()){
    nomad.opts <- .np_nomad_normalize_user_opts(nomad.opts, "npregbw")
    dots <- list(...)
    if (length(nomad.opts))
      .np_nomad_native_reject_unsupported_options_for_route(
        opts = nomad.opts,
        route = "native npreg NOMAD route",
        bwsolver = bwsolver
      )
    elapsed.start <- proc.time()[3]

    xdat <- toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(dim(xdat)[2])
    }
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    bwsolver <- npValidateBwsolver(bwsolver)
    remin <- npValidateScalarLogical(powell.remin, "powell.remin")
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

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (length(bws$bw) != dim(xdat)[2])
      stop("length of bandwidth vector does not match number of columns of 'xdat'")

    npValidateRegressionNnLowerBound(
      bws,
      where = "npregbw",
      allow.zero.placeholder = TRUE
    )
    npValidateRegressionExtendedNn(
      bws,
      where = "npregbw",
      bandwidth.compute = bandwidth.compute
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
    xdat.frame <- xdat

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
    reg.spec <- npCanonicalConditionalRegSpec(
      regtype = tbw$regtype,
      basis = tbw$basis,
      degree = tbw$degree,
      bernstein.basis = tbw$bernstein.basis,
      ncon = tbw$ncon,
      where = "npregbw"
    )

    mysd <- EssDee(rcon)
    nconfac <- nrow^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- nrow^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

    reg.c <- npRegtypeToC(regtype = reg.spec$regtype.engine,
                          degree = reg.spec$degree.engine,
                          ncon = tbw$ncon,
                          context = "npregbw")
    npCheckRegressionDesignCondition(reg.code = reg.c$code,
                                     xcon = rcon,
                                     basis = reg.spec$basis.engine,
                                     degree = reg.spec$degree.engine,
                                     bernstein.basis = reg.spec$bernstein.basis.engine,
                                     where = "npregbw")
    degree.c <- if (tbw$ncon > 0) {
      as.integer(if (is.null(reg.c$degree)) rep.int(0L, tbw$ncon) else reg.c$degree)
    } else {
      integer(1)
    }

    if (bandwidth.compute){
      if (npBwsolverUsesMads(bwsolver)) {
        opt.args <- list(
          nmulti = nmulti,
          itmax = itmax,
          powell.remin = remin,
          scale.init.categorical.sample = scale.init.categorical.sample,
          ftol = ftol,
          tol = tol,
          small = small,
          lbc.dir = lbc.dir,
          dfc.dir = dfc.dir,
          cfac.dir = cfac.dir,
          initc.dir = initc.dir,
          lbd.dir = lbd.dir,
          hbd.dir = hbd.dir,
          dfac.dir = dfac.dir,
          initd.dir = initd.dir,
          scale.factor.init.lower = scale.factor.init.lower,
          scale.factor.init.upper = scale.factor.init.upper,
          scale.factor.init = scale.factor.init,
          lbd.init = lbd.init,
          hbd.init = hbd.init,
          dfac.init = dfac.init,
          scale.factor.search.lower = scale.factor.search.lower,
          invalid.penalty = invalid.penalty,
          penalty.multiplier = penalty.multiplier,
          transform.bounds = transform.bounds,
          bandwidth.compute = TRUE,
          bwsolver = bwsolver,
          nomad.opts = nomad.opts
        )
        return(.npregbw_run_fixed_degree_mads(
          xdat = xdat.frame,
          ydat = ydat,
          bws = tbw$bw,
          reg.args = list(
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
            scale.factor.search.lower = scale.factor.search.lower
          ),
          opt.args = opt.args,
          yname = tbw$ynames,
          bwsolver = bwsolver
        ))
      }
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
        int_do_tree = .npregbw_tree_code(bws, ncon = bws$ncon, ncat = bws$nuno + bws$nord),
        scale.init.categorical.sample = scale.init.categorical.sample,
        dfc.dir = dfc.dir,
        transform.bounds = transform.bounds)
      
      myoptd = list(ftol=ftol, tol=tol, small=small,
        lbc.dir = lbc.dir, cfac.dir = cfac.dir, initc.dir = initc.dir, 
        lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir, 
        lbc.init = cont.start$scale.factor.init.lower, hbc.init = cont.start$scale.factor.init.upper, cfac.init = cont.start$scale.factor.init, 
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
                as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
                as.integer(npLpBasisCode(reg.spec$basis.engine)),
                as.double(cker.bounds.c$lb),
                as.double(cker.bounds.c$ub),
                PACKAGE = "np"))[1]
      

      rorder = numeric(ncol)
      ord_idx <- seq_len(ncol)
      rorder[c(ord_idx[bws$icon], ord_idx[bws$iuno], ord_idx[bws$iord])] <- ord_idx

      tbw$bw <- myout$bw[rorder]
      tbw$fval <- myout$fval[1]
      tbw$ifval <- myout$fval[2]
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
      tbw$nn.cache <- .np_nn_cache_stats(myout$nn.cache)
      tbw$fval.history <- myout$fval.history
      tbw$eval.history <- myout$eval.history
      tbw$invalid.history <- myout$invalid.history
      tbw$timing <- myout$timing
    }

    tbw$total.time <- proc.time()[3] - elapsed.start

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
                      nn.cache = tbw$nn.cache,
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

.npregbw_tree_code <- function(bws, ncon, ncat) {
  code <- npDoTreeOrCategoricalCompress(ncon = ncon, ncat = ncat, bws = bws)

  if (!identical(code, DO_TREE_YES))
    return(code)

  method <- if (!is.null(bws$method) && length(bws$method)) {
    as.character(bws$method[1L])
  } else {
    "cv.ls"
  }
  bwtype <- if (!is.null(bws$type) && length(bws$type)) {
    as.character(bws$type[1L])
  } else {
    "fixed"
  }
  regtype <- if (!is.null(bws$regtype.engine) && length(bws$regtype.engine)) {
    as.character(bws$regtype.engine[1L])
  } else if (!is.null(bws$regtype) && length(bws$regtype)) {
    as.character(bws$regtype[1L])
  } else {
    "lc"
  }

  if (ncon > 0L &&
      bwtype %in% c("generalized_nn", "adaptive_nn") &&
      !identical(regtype, "lc")) {
    return(DO_TREE_NO)
  }

  code
}

.npregbw_eval_only <- function(xdat,
                               ydat,
                               bws,
                               invalid.penalty = c("baseline", "dbmax"),
                               penalty.multiplier = 10,
                               objective = c("ls", "ks")) {
  objective <- match.arg(objective)
  out <- .npregbw_call_fixed_degree_core(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    invalid.penalty = invalid.penalty,
    penalty.multiplier = penalty.multiplier,
    eval.only = TRUE,
    objective = objective
  )

  list(
    objective = out$objective,
    num.feval = out$num.feval,
    num.feval.fast = out$num.feval.fast
  )
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
                                            eval.only = FALSE,
                                            objective = c("ls", "ks")) {
  invalid.penalty <- match.arg(invalid.penalty)
  objective <- match.arg(objective)
  if (identical(objective, "ks") && !isTRUE(eval.only))
    stop("internal Klein-Spady regression objective is eval-only", call. = FALSE)
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
  reg.spec <- npCanonicalConditionalRegSpec(
    regtype = bws$regtype,
    basis = bws$basis,
    degree = bws$degree,
    bernstein.basis = bws$bernstein.basis,
    ncon = bws$ncon,
    where = "npregbw"
  )
  reg.c <- npRegtypeToC(regtype = reg.spec$regtype.engine,
                        degree = reg.spec$degree.engine,
                        ncon = bws$ncon,
                        context = "npregbw")
  degree.c <- if (bws$ncon > 0) as.integer(reg.spec$degree.engine) else integer(1L)
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
    bwmethod = if (identical(objective, "ks")) RBWM_CVKS else switch(bws$method,
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
    int_do_tree = .npregbw_tree_code(bws, ncon = bws$ncon, ncat = bws$nuno + bws$nord),
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
    out <- .Call(
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
      as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
      as.integer(npLpBasisCode(reg.spec$basis.engine)),
      as.double(cker.bounds.c$lb),
      as.double(cker.bounds.c$ub),
      PACKAGE = "np"
    )
  } else {
    out <- .Call(
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
      as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
      as.integer(npLpBasisCode(reg.spec$basis.engine)),
      as.double(cker.bounds.c$lb),
      as.double(cker.bounds.c$ub),
      PACKAGE = "np"
    )
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
    nn.cache = .np_nn_cache_stats(out$nn.cache),
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
  tbw$nn.cache <- core$nn.cache
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
    remin = opt.value("powell.remin", TRUE),
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

.npregbw_run_fixed_degree_mads <- function(xdat,
                                           ydat,
                                           bws,
                                           reg.args,
                                           opt.args,
                                           yname,
                                           bwsolver = c("mads", "mads+powell")) {
  bwsolver <- npValidateBwsolver(bwsolver)
  template <- .npregbw_build_rbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    reg.args = reg.args,
    yname = yname
  )
  if (!(template$type %in% c("fixed", "generalized_nn", "adaptive_nn")))
    stop("bwsolver='mads' requires bwtype='fixed', 'generalized_nn', or 'adaptive_nn'")

  setup <- .npregbw_nomad_bw_setup(xdat = xdat, template = template, allow.extended.nn = TRUE)
  bounds <- .npregbw_nomad_bw_bounds(template = template, setup = setup)
  point.start <- if (all(template$bw == 0)) {
    NULL
  } else {
    .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup)
  }
  x0 <- .npregbw_nomad_complete_bw_start_point(
    point = point.start,
    bounds = bounds,
    setup = setup
  )

  opt.value <- function(name, default) {
    if (is.null(opt.args[[name]])) default else opt.args[[name]]
  }
  mads.nmulti <- opt.value("nmulti", npDefaultNmulti(dim(toFrame(xdat))[2L]))
  mads.inner.nmulti <- opt.value("mads.nmulti", opt.value("nomad.nmulti", 0L))
  random.seed <- opt.value("random.seed", 42L)
  nomad.opts <- opt.value("nomad.opts", list())
  source <- "explicit"
  reason <- NULL
  mads.num.feval.total <- 0
  mads.num.feval.fast.total <- 0

  eval_fun <- function(point) {
    bw_vec <- .npregbw_nomad_point_to_bw(point, template = template, setup = setup)
    tbw <- .npregbw_build_rbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = reg.args,
      yname = yname
    )
    out <- .npregbw_eval_only(
      xdat = xdat,
      ydat = ydat,
      bws = tbw,
      invalid.penalty = opt.value("invalid.penalty", "baseline"),
      penalty.multiplier = opt.value("penalty.multiplier", 10)
    )
    mads.num.feval.total <<- mads.num.feval.total + as.numeric(out$num.feval[1L])
    mads.num.feval.fast.total <<- mads.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])

    list(
      objective = out$objective,
      degree = integer(0L),
      num.feval = out$num.feval,
      num.feval.fast = out$num.feval.fast
    )
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    bw_vec <- .npregbw_nomad_point_to_bw(point, template = template, setup = setup)
    final.tbw <- .npregbw_build_rbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = reg.args,
      yname = yname
    )
    powell.elapsed <- NA_real_
    direct.payload <- .npregbw_finalize_fixed_degree_payload(
      xdat = xdat,
      ydat = ydat,
      bws = final.tbw,
      core = list(
        bw = as.numeric(final.tbw$bw),
        objective = as.numeric(best_record$objective),
        ifval = as.numeric(best_record$objective),
        num.feval = as.numeric(mads.num.feval.total),
        num.feval.fast = as.numeric(mads.num.feval.fast.total),
        fval.history = as.numeric(best_record$objective),
        eval.history = if (!is.null(solution$bbe)) rep(1, max(1L, as.integer(solution$bbe))) else 1,
        invalid.history = 0,
        timing = NA_real_
      ),
      total.time = NA_real_
    )
    direct.objective <- as.numeric(best_record$objective)

    if (identical(bwsolver, "mads+powell")) {
      hot.opt.args <- .np_nomad_powell_hotstart_opt_args(
        opt.args,
        strategy = "disable_multistart",
        remin = isTRUE(opt.args$powell.remin)
      )
      powell.start <- proc.time()[3L]
      hot.payload <- .npregbw_with_powell_refinement_progress(integer(0L), local({
        .npregbw_run_fixed_degree_source_of_truth(
          xdat = xdat,
          ydat = ydat,
          bws = final.tbw,
          opt.args = hot.opt.args
        )
      }))
      powell.elapsed <- proc.time()[3L] - powell.start
      hot.payload$num.feval <- as.numeric(direct.payload$num.feval[1L]) + as.numeric(hot.payload$num.feval[1L])
      hot.payload$num.feval.fast <- as.numeric(direct.payload$num.feval.fast[1L]) + as.numeric(hot.payload$num.feval.fast[1L])
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min")) {
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  native.start.bounds <- .np_nomad_bw_restart_start_bounds(
    bounds = bounds,
    setup = setup,
    opt.value = opt.value,
    where = "npregbw"
  )
  if (is.null(point.start)) {
    x0 <- .npregbw_nomad_complete_bw_start_point(
      point = NULL,
      bounds = bounds,
      setup = setup,
      initial = native.start.bounds$initial,
      where = "npregbw"
    )
  }

  if (.npregbw_nomad_native_target(template, bwsolver)) {
    .npregbw_nomad_native_require_crs()
    native.nmulti <- npValidateNmulti(mads.nmulti)
    native.inner.nmulti <- npValidateNonNegativeInteger(mads.inner.nmulti, "nomad.nmulti")
    native.inner.nmulti <- as.integer(native.inner.nmulti[1L])
    if (isTRUE(opt.args$nomad.remin))
      stop("native npreg NOMAD route does not support NOMAD remin", call. = FALSE)

    native.nomad.opts <- .np_nomad_prepare_solver_opts(
      random.seed = random.seed,
      nomad.opts = nomad.opts,
      geometry.policy = "user-only",
      where = "npregbw native NOMAD source geometry"
    )
    native.option.vectors <- .npregbw_nomad_native_option_vectors(native.nomad.opts)
    native.start.matrix <- .np_nomad_build_starts(
      x0 = x0,
      bbin = bounds$bbin,
      lb = bounds$lower,
      ub = bounds$upper,
      nmulti = native.nmulti,
      random.seed = random.seed,
      degree_spec = NULL,
      start.lower = native.start.bounds$lower,
      start.upper = native.start.bounds$upper
    )
    native.prep <- .npregbw_nomad_native_prepare_args(
      xdat = xdat,
      ydat = ydat,
      bws = template,
      invalid.penalty = opt.value("invalid.penalty", "baseline"),
      penalty.multiplier = opt.value("penalty.multiplier", 10),
      itmax = opt.value("itmax", 10000L),
      ftol = opt.value("ftol", 1.490116e-07),
      tol = opt.value("tol", 1.490116e-04),
      small = opt.value("small", 1.490116e-05),
      scale.factor.search.lower = opt.value("scale.factor.search.lower", NULL),
      scale.init.categorical.sample = opt.value("scale.init.categorical.sample", FALSE),
      transform.bounds = opt.value("transform.bounds", FALSE)
    )

    native.results <- vector("list", nrow(native.start.matrix))
    native.best.index <- NA_integer_
    native.best.objective <- Inf
    native.nomad.elapsed <- 0
    native.num.feval.total <- 0
    native.num.feval.fast.total <- 0
    native.num.feval.invalid.total <- 0
    for (i in seq_len(nrow(native.start.matrix))) {
      native.start <- proc.time()[3L]
      native.i <- npNomadNativeSearchRegression(
        prep = native.prep,
        x0 = as.numeric(native.start.matrix[i, ]),
        bbin = bounds$bbin,
        lb = bounds$lower,
        ub = bounds$upper,
        max.eval = 0L,
        random.seed = random.seed,
        inner.start.count = native.inner.nmulti,
        option.names = native.option.vectors$names,
        option.values = native.option.vectors$values
      )
      native.elapsed <- proc.time()[3L] - native.start
      native.nomad.elapsed <- native.nomad.elapsed + native.elapsed
      if (!identical(as.integer(native.i$status[1L]), 0L) ||
          !identical(as.integer(native.i$result_status[1L]), 0L)) {
        stop(sprintf(
          "native npreg NOMAD route failed (status=%s, result_status=%s): %s",
          as.integer(native.i$status[1L]),
          as.integer(native.i$result_status[1L]),
          as.character(native.i$message[1L])
        ), call. = FALSE)
      }
      official.objective.i <- as.numeric(native.i$official_objective[1L])
      objective.i <- as.numeric(native.i$objective[1L])
      native.results[[i]] <- list(
        restart = i,
        start = as.numeric(native.start.matrix[i, ]),
        elapsed = native.elapsed,
        status = "ok",
        message = as.character(native.i$message[1L]),
        objective = official.objective.i,
        bbe = as.numeric(native.i$blackbox_evaluations[1L]),
        iterations = as.numeric(native.i$iterations[1L]),
        solution = as.numeric(native.i$solution),
        best_point = as.numeric(native.i$best_point),
        best_objective = objective.i,
        native = native.i
      )
      native.num.feval.total <- native.num.feval.total + as.numeric(native.i$total_num.feval[1L])
      native.num.feval.fast.total <- native.num.feval.fast.total + as.numeric(native.i$total_num.feval.fast[1L])
      native.num.feval.invalid.total <- native.num.feval.invalid.total + as.numeric(native.i$total_num.feval.invalid[1L])
      if (is.finite(objective.i) && objective.i < native.best.objective) {
        native.best.objective <- objective.i
        native.best.index <- i
      }
    }
    if (!is.finite(native.best.index))
      stop("native npreg NOMAD route did not return a finite solution", call. = FALSE)

    native.best <- native.results[[native.best.index]]
    native.handoff.point <- as.numeric(native.best$best_point)
    if (any(!is.finite(native.handoff.point)))
      stop("native npreg NOMAD route did not return a finite best point", call. = FALSE)
    native.bw <- .npregbw_nomad_point_to_bw(native.handoff.point, template = template, setup = setup)
    native.record <- list(
      eval_id = as.integer(native.best$native$compiled_callback_calls[1L]),
      degree = integer(0L),
      objective = native.best.objective,
      status = "ok",
      cached = FALSE,
      message = native.best$message,
      elapsed = native.best$elapsed,
      num.feval = as.numeric(native.best$native$best_num.feval[1L]),
      num.feval.fast = as.numeric(native.best$native$best_num.feval.fast[1L]),
      num.feval.invalid = as.numeric(native.best$native$best_num.feval.invalid[1L])
    )
    mads.num.feval.total <- native.num.feval.total
    mads.num.feval.fast.total <- native.num.feval.fast.total
    payload.result <- build_payload(
      point = native.handoff.point,
      best_record = native.record,
      solution = native.best,
      interrupted = FALSE
    )
    search.result <- list(
      best = native.record,
      best_point = native.handoff.point,
      best_payload = payload.result$payload,
      completed = TRUE,
      method = "nomad",
      source = source,
      reason = reason,
      restart.results = native.results,
      best.restart = native.best.index,
      nomad.time = native.nomad.elapsed,
      powell.time = payload.result$powell.time,
      optim.time = native.nomad.elapsed + as.numeric(payload.result$powell.time[1L]),
      num.feval.total = native.num.feval.total,
      num.feval.fast.total = native.num.feval.fast.total,
      num.feval.invalid.total = native.num.feval.invalid.total,
      native.diagnostics = list(
        raw.point = native.handoff.point,
        bandwidth = native.bw,
        objective = native.best.objective,
        official.solution = as.numeric(native.best$solution),
        official.objective = as.numeric(native.best$objective[1L]),
        compiled.callback.count = as.integer(native.best$native$compiled_callback_calls[1L]),
        compiled.callback.failures = as.integer(native.best$native$compiled_callback_failures[1L]),
        crs.callback.evaluations = as.integer(native.best$native$crs_callback_evaluations[1L]),
        blackbox.evaluations = as.integer(native.best$native$blackbox_evaluations[1L]),
        cache.hits = as.integer(native.best$native$cache_hits[1L]),
        cache.size = as.integer(native.best$native$cache_size[1L]),
        total.evaluations = as.integer(native.best$native$total_evaluations[1L]),
        iterations = as.integer(native.best$native$iterations[1L])
      )
    )
    if (isTRUE(getOption("np.developer.native.nomad.diagnostics", FALSE)) &&
        !is.null(search.result$best_payload))
      attr(search.result$best_payload, "native.nomad.diagnostics") <- search.result$native.diagnostics
    if (!is.null(payload.result$objective) &&
        .np_degree_better(payload.result$objective, search.result$best$objective, direction = "min"))
      search.result$best$objective <- as.numeric(payload.result$objective[1L])
  } else {
    search.result <- .np_nomad_search(
      engine = "nomad",
      baseline_record = NULL,
      start_degree = integer(0L),
      x0 = x0,
      bbin = bounds$bbin,
      lb = bounds$lower,
      ub = bounds$upper,
      eval_fun = eval_fun,
      build_payload = build_payload,
      direction = "min",
      objective_name = "fval",
      nmulti = mads.nmulti,
      nomad.inner.nmulti = mads.inner.nmulti,
      random.seed = random.seed,
      handoff_before_build = identical(bwsolver, "mads+powell"),
      remin = isTRUE(opt.args$nomad.remin),
      nomad.opts = nomad.opts,
      start.lower = native.start.bounds$lower,
      start.upper = native.start.bounds$upper
    )
  }
  search.result$method <- bwsolver
  payload <- search.result$best_payload
  payload$bwsolver <- bwsolver
  payload$search.engine <- bwsolver
  if (!is.null(search.result$nomad.time))
    payload$nomad.time <- as.numeric(search.result$nomad.time[1L])
  if (!is.null(search.result$powell.time))
    payload$powell.time <- as.numeric(search.result$powell.time[1L])
  payload$total.time <- as.numeric(search.result$optim.time[1L])
  payload <- .np_attach_nomad_restart_summary(payload, search.result)
  payload
}

.npregbw_nomad_controls <- function(search.engine) {
  .np_degree_search_engine_controls(search.engine)
}

.npregbw_nomad_bw_setup <- function(xdat,
                                    template,
                                    bandwidth.scale.categorical = 1e4,
                                    allow.extended.nn = FALSE,
                                    evaldat = NULL) {
  xdat <- toFrame(xdat)
  evaldat <- if (is.null(evaldat)) xdat else toFrame(evaldat)
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

  cont_extendednn_upper <- if (isTRUE(allow.extended.nn)) {
    npContinuousExtendedNnNomadUpper(
      traindat = xdat,
      evaldat = evaldat,
      bwtype = template$type,
      ckertype = template$ckertype,
      cont.idx = cont_idx
    )
  } else {
    NULL
  }

  list(
    type = as.character(template$type[1L]),
    nobs = nrow,
    cont_scale = cont_scale,
    cont_idx = cont_idx,
    cat_idx = cat_idx,
    ncatfac = ncatfac,
    bandwidth.scale.categorical = bandwidth.scale.categorical,
    cat_upper = cat_upper,
    cont_extendednn_upper = cont_extendednn_upper
  )
}

.npregbw_nomad_bw_bounds <- function(template, setup) {
  .np_nomad_bw_bounds(
    template = template,
    setup = setup,
    fixed.lower = npGetScaleFactorSearchLower(
      template,
      argname = "template$scale.factor.search.lower"
    ),
    nn.lower = npRegressionNnLowerBound(template),
    where = "npregbw"
  )
}

.np_nomad_bw_restart_start_bounds <- function(bounds,
                                             setup,
                                             opt.value,
                                             where = "NOMAD bandwidth search") {
  start.lower <- as.numeric(bounds$lower)
  start.upper <- as.numeric(bounds$upper)
  start.initial <- rep(NA_real_, length(start.lower))
  n <- length(start.lower)
  if (length(start.upper) != n)
    stop(sprintf("%s: search bounds must have matching lengths", where), call. = FALSE)

  ncon <- if (!is.null(bounds$ncon)) {
    as.integer(bounds$ncon[1L])
  } else if (!is.null(setup$cont_idx)) {
    length(setup$cont_idx)
  } else if (!is.null(setup$cont_flat)) {
    length(setup$cont_flat)
  } else {
    0L
  }
  ncat <- if (!is.null(bounds$ncat)) {
    as.integer(bounds$ncat[1L])
  } else if (!is.null(setup$cat_idx)) {
    length(setup$cat_idx)
  } else if (!is.null(setup$cat_flat)) {
    length(setup$cat_flat)
  } else {
    max(0L, n - ncon)
  }

  type <- if (!is.null(setup$type)) as.character(setup$type[1L]) else ""
  if (ncon > 0L && identical(type, "fixed")) {
    cont.idx <- seq_len(ncon)
    search.lower <- max(as.numeric(bounds$lower[cont.idx]), na.rm = TRUE)
    cont.start <- npContinuousSearchStartControls(
      opt.value("scale.factor.init.lower", 0.1),
      opt.value("scale.factor.init.upper", 2.0),
      opt.value("scale.factor.init", 0.5),
      search.lower,
      where = where
    )
    start.lower[cont.idx] <- cont.start$scale.factor.init.lower
    start.upper[cont.idx] <- cont.start$scale.factor.init.upper
    start.initial[cont.idx] <- cont.start$scale.factor.init
  }

  if (ncat > 0L) {
    cat.idx <- ncon + seq_len(ncat)
    lbd.init <- npValidatePositiveFiniteNumeric(opt.value("lbd.init", 0.1), "lbd.init")
    hbd.init <- npValidatePositiveFiniteNumeric(opt.value("hbd.init", 0.9), "hbd.init")
    dfac.init <- npValidatePositiveFiniteNumeric(opt.value("dfac.init", 0.375), "dfac.init")
    if (hbd.init < lbd.init)
      stop(sprintf("%s: 'hbd.init' must be greater than or equal to 'lbd.init'", where),
           call. = FALSE)
    cat.scale <- if (!is.null(setup$bandwidth.scale.categorical)) {
      as.double(setup$bandwidth.scale.categorical[1L])
    } else {
      1e4
    }
    start.lower[cat.idx] <- as.double(lbd.init) * cat.scale
    start.upper[cat.idx] <- as.double(hbd.init) * cat.scale
    start.initial[cat.idx] <- as.double(dfac.init) * cat.scale
  }

  finite.lower <- is.finite(start.lower) & is.finite(bounds$lower)
  finite.upper <- is.finite(start.upper) & is.finite(bounds$upper)
  start.lower[finite.lower] <- pmax(start.lower[finite.lower], bounds$lower[finite.lower])
  start.upper[finite.upper] <- pmin(start.upper[finite.upper], bounds$upper[finite.upper])

  both.finite <- is.finite(start.lower) & is.finite(start.upper)
  inverted <- both.finite & start.upper < start.lower
  if (any(inverted)) {
    stop(sprintf(
      "%s: effective NOMAD random-start interval is empty after applying search bounds",
      where
    ), call. = FALSE)
  }

  finite.initial <- is.finite(start.initial)
  start.initial[finite.initial] <- pmax(bounds$lower[finite.initial],
                                        pmin(bounds$upper[finite.initial],
                                             start.initial[finite.initial]))

  list(lower = start.lower, upper = start.upper, initial = start.initial)
}

.np_nomad_explicit_or_initial_start <- function(point,
                                                initial,
                                                n,
                                                where = "NOMAD start") {
  if (!is.null(point))
    return(point)
  if (is.null(initial))
    return(NULL)

  initial <- as.numeric(initial)
  if (length(initial) != n)
    stop(sprintf("%s: initial start vector must have length %d", where, n),
         call. = FALSE)
  initial
}

.npregbw_nomad_complete_bw_start_point <- function(point,
                                                   bounds,
                                                   setup,
                                                   initial = NULL,
                                                   where = "npregbw") {
  .np_nomad_bw_complete_start_point(
    point = point,
    bounds = bounds,
    template = list(type = setup$type),
    setup = setup,
    initial = initial,
    where = where
  )
}

.npregbw_nomad_point_to_bw <- function(point, template, setup) {
  .np_nomad_bw_point_to_storage(
    point = point,
    template = template,
    setup = setup,
    storage.length = length(template$bw),
    clamp.nn = FALSE
  )
}

.npregbw_nomad_bw_to_point <- function(bws, template, setup) {
  .np_nomad_bw_storage_to_point(bws = bws, template = template, setup = setup)
}

.npregbw_powell_progress_fields <- function(state,
                                            done = NULL,
                                            detail = NULL,
                                            now = .np_progress_now()) {
  fields <- .np_degree_progress_context_fields()
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
  reuse.active <- !is.null(old.state) &&
    isTRUE(old.state$enabled) &&
    isTRUE(old.state$visible)
  active.state <- if (isTRUE(reuse.active)) {
    old.state
  } else {
    .np_progress_begin(
      label = .np_nomad_powell_progress_label(),
      domain = "general",
      surface = "bandwidth"
    )
  }

  on.exit({
    if (!isTRUE(reuse.active)) {
      current.state <- .np_progress_runtime$bandwidth_state
      if (!is.null(current.state)) {
        .np_progress_end(current.state)
      }
    }
    .np_progress_runtime$bandwidth_state <- old.state
  }, add = TRUE)

  active.state$label <- .np_nomad_powell_progress_label()
  active.state$unknown_total_fields <- .npregbw_powell_progress_fields
  active.state$nomad_nmulti <- 1L
  active.state$nomad_current_degree <- as.integer(degree)
  active.state$started <- .np_progress_now()
  active.state$last_done <- NULL
  active.state <- .np_progress_show_now(active.state)
  .np_progress_runtime$bandwidth_state <- active.state

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

.npregbw_nomad_search <- function(xdat,
                                  ydat,
                                  bws,
                                  reg.args,
                                  opt.args,
                                  yname,
                                  degree.search,
                                  nomad.inner.nmulti = 0L,
                                  random.seed = 42L,
                                  nomad.opts = list(),
                                  source = "explicit",
                                  reason = NULL,
                                  progress_label = NULL) {
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
  if (!(template$type %in% c("fixed", "generalized_nn", "adaptive_nn")))
    stop("automatic degree search with search.engine='nomad' requires bwtype='fixed', 'generalized_nn', or 'adaptive_nn'")

  setup <- .npregbw_nomad_bw_setup(xdat = xdat, template = template, allow.extended.nn = TRUE)
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)
  opt.value <- function(name, default) {
    if (is.null(opt.args[[name]])) default else opt.args[[name]]
  }
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(dim(xdat)[2]) else npValidateNmulti(opt.args$nmulti[1L])

  bw_bounds <- .npregbw_nomad_bw_bounds(template = template, setup = setup)
  bw_start_bounds <- .np_nomad_bw_restart_start_bounds(
    bounds = bw_bounds,
    setup = setup,
    opt.value = opt.value,
    where = "npregbw"
  )

  point.start <- if (all(template$bw == 0)) {
    NULL
  } else {
    .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup)
  }
  x0 <- c(
    .npregbw_nomad_complete_bw_start_point(
      point = point.start,
      bounds = bw_bounds,
      setup = setup,
      initial = bw_start_bounds$initial,
      where = "npregbw"
    ),
    as.integer(degree.search$start.degree)
  )
  lb <- c(bw_bounds$lower, degree.search$lower)
  ub <- c(bw_bounds$upper, degree.search$upper)
  bbin <- c(bw_bounds$bbin, rep.int(1L, ncon))
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  .np_nomad_baseline_note(degree.search$start.degree)

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[ncon + ncat + seq_len(ncon)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)

    eval.reg.args <- reg.args
    eval.reg.args$regtype <- "lp"
    eval.reg.args$degree <- degree
    eval.reg.args$bernstein.basis <- degree.search$bernstein.basis

    tbw <- .npregbw_build_rbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = eval.reg.args,
      yname = yname
    )

    out <- .npregbw_eval_only(
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
        eval.history = if (!is.null(solution$bbe)) rep(1, max(1L, as.integer(solution$bbe))) else 1,
        invalid.history = 0,
        timing = NA_real_
      ),
      total.time = NA_real_
    )
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.opt.args <- .np_nomad_powell_hotstart_opt_args(
        opt.args,
        strategy = "disable_multistart",
        remin = isTRUE(opt.args$powell.remin)
      )

      powell.start <- proc.time()[3L]
      hot.payload <- .npregbw_with_powell_refinement_progress(degree, local({
        .npregbw_run_fixed_degree_source_of_truth(
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
      if (!is.null(hot.payload$method) && length(hot.payload$method))
        hot.payload$pmethod <- bwmToPrint(as.character(hot.payload$method[1L]))

      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min")) {
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  search.engine.used <- if (identical(degree.search$engine, "nomad+powell")) {
    "nomad"
  } else {
    degree.search$engine
  }

  if (.npregbw_nomad_degree_native_target(template, degree.search)) {
    .npregbw_nomad_native_require_crs()
    native.nmulti <- npValidateNmulti(nomad.nmulti)
    native.inner.nmulti <- npValidateNonNegativeInteger(nomad.inner.nmulti, "nomad.inner.nmulti")
    native.inner.nmulti <- as.integer(native.inner.nmulti[1L])

    native.nomad.opts <- .np_nomad_prepare_solver_opts(
      random.seed = random.seed,
      nomad.opts = nomad.opts,
      coordinate.roles = .np_nomad_coordinate_roles(bw_bounds, degree.search),
      expected.length = length(lb),
      geometry.policy = "generate-central",
      where = "npregbw native NOMAD degree source geometry"
    )
    native.option.vectors <- .npregbw_nomad_native_option_vectors(native.nomad.opts)
    native.start.matrix <- .np_nomad_build_starts(
      x0 = x0,
      bbin = bbin,
      lb = lb,
      ub = ub,
      nmulti = native.nmulti,
      random.seed = random.seed,
      start.lower = c(bw_start_bounds$lower, degree.search$lower),
      start.upper = c(bw_start_bounds$upper, degree.search$upper),
      degree_spec = list(
        initial = degree.search$start.degree,
        lower = degree.search$lower,
        upper = degree.search$upper,
        basis = degree.search$basis,
        nobs = degree.search$nobs,
        user_supplied = degree.search$start.user
      )
    )
    native.prep <- .npregbw_nomad_native_prepare_args(
      xdat = xdat,
      ydat = ydat,
      bws = template,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier,
      itmax = if (is.null(opt.args$itmax)) 10000L else opt.args$itmax,
      ftol = if (is.null(opt.args$ftol)) 1.490116e-07 else opt.args$ftol,
      tol = if (is.null(opt.args$tol)) 1.490116e-04 else opt.args$tol,
      small = if (is.null(opt.args$small)) 1.490116e-05 else opt.args$small,
      scale.factor.search.lower = if (is.null(opt.args$scale.factor.search.lower)) NULL else opt.args$scale.factor.search.lower,
      scale.init.categorical.sample = if (is.null(opt.args$scale.init.categorical.sample)) FALSE else opt.args$scale.init.categorical.sample,
      transform.bounds = if (is.null(opt.args$transform.bounds)) FALSE else opt.args$transform.bounds
    )

    ndeg <- length(degree.search$start.degree)
    degree.idx <- (ncol(native.start.matrix) - ndeg + 1L):ncol(native.start.matrix)
    native.results <- vector("list", nrow(native.start.matrix))
    native.best.index <- NA_integer_
    native.best.objective <- Inf
    native.nomad.elapsed <- 0
    native.num.feval.total <- 0
    native.num.feval.fast.total <- 0
    native.num.feval.invalid.total <- 0
    native.callback.total <- 0L
    native.baseline.record <- NULL
    native.progress <- .np_nomad_native_progress_begin(
      nmulti = native.nmulti,
      baseline_degree = degree.search$start.degree,
      best_record = native.baseline.record,
      label = progress_label
    )
    on.exit(.np_nomad_native_progress_abort(native.progress), add = TRUE)

    make_native_record <- function(native, objective, degree, elapsed) {
      list(
        eval_id = as.integer(native$compiled_callback_calls[1L]),
        degree = as.integer(degree),
        objective = as.numeric(objective[1L]),
        status = "ok",
        cached = FALSE,
        message = as.character(native$message[1L]),
        elapsed = as.numeric(elapsed[1L]),
        num.feval = as.numeric(native$best_num.feval[1L]),
        num.feval.fast = as.numeric(native$best_num.feval.fast[1L])
      )
    }

    run_native_restart <- function(start, restart.index, remin = FALSE) {
      native.restart.degree <- if (ndeg > 0L) {
        as.integer(round(start[degree.idx]))
      } else {
        integer(0L)
      }
      .np_nomad_native_progress_restart(
        handle = native.progress,
        restart_index = restart.index,
        degree = native.restart.degree,
        best_record = native.baseline.record,
        restart_durations = native.nomad.elapsed,
        eval_offset = native.callback.total
      )
      native.start <- proc.time()[3L]
      native <- npNomadNativeSearchRegression(
        prep = native.prep,
        x0 = as.numeric(start),
        bbin = bbin,
        lb = lb,
        ub = ub,
        max.eval = 0L,
        random.seed = random.seed,
        inner.start.count = native.inner.nmulti,
        option.names = native.option.vectors$names,
        option.values = native.option.vectors$values
      )
      native.elapsed <- proc.time()[3L] - native.start
      if (!identical(as.integer(native$status[1L]), 0L) ||
          !identical(as.integer(native$result_status[1L]), 0L)) {
        stop(sprintf(
          "native npreg NOMAD degree-search route failed (status=%s, result_status=%s): %s",
          as.integer(native$status[1L]),
          as.integer(native$result_status[1L]),
          as.character(native$message[1L])
        ), call. = FALSE)
      }
      if (is.null(native$best_point) || any(!is.finite(native$best_point)))
        stop("native npreg NOMAD degree-search route did not return a finite best point", call. = FALSE)

      native.degree <- if (!is.null(native$best_degree) && length(native$best_degree)) {
        as.integer(native$best_degree)
      } else {
        as.integer(round(native$best_point[degree.idx]))
      }
      native.objective <- as.numeric(native$objective[1L])
      list(
        restart = as.integer(restart.index),
        remin = isTRUE(remin),
        start = as.numeric(start),
        degree.start = native.restart.degree,
        elapsed = native.elapsed,
        status = "ok",
        message = as.character(native$message[1L]),
        objective = native.objective,
        bbe = as.numeric(native$blackbox_evaluations[1L]),
        iterations = as.numeric(native$iterations[1L]),
        solution = as.numeric(native$solution),
        best_point = as.numeric(native$best_point),
        best_degree = native.degree,
        first_degree = if (!is.null(native$first_degree)) as.integer(native$first_degree) else integer(0L),
        first_objective = as.numeric(native$first_objective[1L]),
        native = native
      )
    }

    for (i in seq_len(nrow(native.start.matrix))) {
      native.i <- run_native_restart(
        start = as.numeric(native.start.matrix[i, ]),
        restart.index = i
      )
      native.results[[i]] <- native.i
      native.nomad.elapsed <- native.nomad.elapsed + as.numeric(native.i$elapsed[1L])
      native.num.feval.total <- native.num.feval.total + as.numeric(native.i$native$total_num.feval[1L])
      native.num.feval.fast.total <- native.num.feval.fast.total + as.numeric(native.i$native$total_num.feval.fast[1L])
      native.num.feval.invalid.total <- native.num.feval.invalid.total + as.numeric(native.i$native$total_num.feval.invalid[1L])
      native.callback.total <- native.callback.total + as.integer(native.i$native$compiled_callback_calls[1L])
      if (is.null(native.baseline.record) && length(native.i$first_degree)) {
        native.baseline.record <- list(
          eval_id = 1L,
          degree = as.integer(native.i$first_degree),
          objective = as.numeric(native.i$first_objective[1L]),
          status = "ok",
          cached = FALSE,
          message = native.i$message,
          elapsed = native.i$elapsed,
          num.feval = NA_real_
        )
      }
      if (is.finite(native.i$objective) &&
          .np_degree_better(native.i$objective, native.best.objective, direction = "min")) {
        native.best.objective <- native.i$objective
        native.best.index <- i
      }
    }
    if (!is.finite(native.best.index))
      stop("native npreg NOMAD degree-search route did not return a finite solution", call. = FALSE)

    if (isTRUE(opt.args$nomad.remin)) {
      remin.index <- length(native.results) + 1L
      remin.start <- as.numeric(native.results[[native.best.index]]$best_point)
      native.remin <- run_native_restart(
        start = remin.start,
        restart.index = remin.index,
        remin = TRUE
      )
      native.results[[remin.index]] <- native.remin
      native.nomad.elapsed <- native.nomad.elapsed + as.numeric(native.remin$elapsed[1L])
      native.num.feval.total <- native.num.feval.total + as.numeric(native.remin$native$total_num.feval[1L])
      native.num.feval.fast.total <- native.num.feval.fast.total + as.numeric(native.remin$native$total_num.feval.fast[1L])
      native.num.feval.invalid.total <- native.num.feval.invalid.total + as.numeric(native.remin$native$total_num.feval.invalid[1L])
      native.callback.total <- native.callback.total + as.integer(native.remin$native$compiled_callback_calls[1L])
      if (is.finite(native.remin$objective) &&
          .np_degree_better(native.remin$objective, native.best.objective, direction = "min")) {
        native.best.objective <- native.remin$objective
        native.best.index <- remin.index
      }
    }

    native.best <- native.results[[native.best.index]]
    native.record <- make_native_record(
      native = native.best$native,
      objective = native.best$objective,
      degree = native.best$best_degree,
      elapsed = native.best$elapsed
    )
    if (is.null(native.baseline.record))
      native.baseline.record <- native.record
    nomad.num.feval.total <- native.num.feval.total
    nomad.num.feval.fast.total <- native.num.feval.fast.total
    payload.result <- build_payload(
      point = native.best$best_point,
      best_record = native.record,
      solution = native.best,
      interrupted = FALSE
    )
    .np_nomad_native_progress_end(
      handle = native.progress,
      degree = native.record$degree,
      best_record = native.record
    )
    search.result <- list(
      method = degree.search$engine,
      source = source,
      reason = reason,
      direction = "min",
      verify = FALSE,
      completed = TRUE,
      certified = FALSE,
      interrupted = FALSE,
      baseline = native.baseline.record,
      best = native.record,
      best_payload = payload.result$payload,
      best_point = native.best$best_point,
      n.unique = as.integer(native.callback.total),
      n.visits = as.integer(native.callback.total),
      n.cached = 0L,
      nomad.time = native.nomad.elapsed,
      powell.time = payload.result$powell.time,
      optim.time = sum(c(native.nomad.elapsed, payload.result$powell.time), na.rm = TRUE),
      grid.size = NA_integer_,
      best.restart = native.best.index,
      nomad.remin = isTRUE(opt.args$nomad.remin),
      nomad.remin.index = if (any(vapply(native.results, function(x) isTRUE(x$remin), logical(1)))) {
        which(vapply(native.results, function(x) isTRUE(x$remin), logical(1)))[1L]
      } else {
        NA_integer_
      },
      nomad.remin.roundtrip = NULL,
      restart.starts = lapply(seq_len(nrow(native.start.matrix)), function(i) as.numeric(native.start.matrix[i, ])),
      restart.degree.starts = lapply(seq_len(nrow(native.start.matrix)), function(i) as.integer(native.start.matrix[i, degree.idx])),
      restart.bandwidth.starts = lapply(seq_len(nrow(native.start.matrix)), function(i) as.numeric(native.start.matrix[i, seq_len(degree.idx[1L] - 1L)])),
      restart.start.info = list(
        basis = if (is.null(degree.search$basis)) "glp" else degree.search$basis,
        degree.start.policy = .np_lp_nomad_degree_start_policy(),
        lower = as.integer(degree.search$lower),
        upper = as.integer(degree.search$upper),
        user_supplied_start = isTRUE(degree.search$start.user)
      ),
      restart.results = native.results,
      trace = data.frame(
        trace_id = seq_along(native.results),
        eval_id = vapply(native.results, function(x) as.integer(x$native$compiled_callback_calls[1L]), integer(1L)),
        degree = vapply(native.results, function(x) paste(as.integer(x$best_degree), collapse = ","), character(1L)),
        fval = vapply(native.results, function(x) as.numeric(x$objective[1L]), numeric(1L)),
        status = vapply(native.results, `[[`, character(1L), "status"),
        cached = rep(FALSE, length(native.results)),
        message = vapply(native.results, function(x) if (is.null(x$message)) "" else as.character(x$message[1L]), character(1L)),
        elapsed = vapply(native.results, function(x) as.numeric(x$elapsed[1L]), numeric(1L)),
        num.feval = vapply(native.results, function(x) as.numeric(x$native$best_num.feval[1L]), numeric(1L)),
        stringsAsFactors = FALSE
      ),
      native.diagnostics = list(
        raw.point = as.numeric(native.best$best_point),
        degree = as.integer(native.best$best_degree),
        objective = as.numeric(native.best$objective[1L]),
        official.solution = as.numeric(native.best$solution),
        official.objective = as.numeric(native.best$native$official_objective[1L]),
        compiled.callback.count = as.integer(native.best$native$compiled_callback_calls[1L]),
        compiled.callback.failures = as.integer(native.best$native$compiled_callback_failures[1L]),
        crs.callback.evaluations = as.integer(native.best$native$crs_callback_evaluations[1L]),
        blackbox.evaluations = as.integer(native.best$native$blackbox_evaluations[1L]),
        cache.hits = as.integer(native.best$native$cache_hits[1L]),
        cache.size = as.integer(native.best$native$cache_size[1L]),
        total.evaluations = as.integer(native.best$native$total_evaluations[1L]),
        iterations = as.integer(native.best$native$iterations[1L])
      )
    )
    if (!is.null(payload.result$objective) &&
        .np_degree_better(payload.result$objective, search.result$best$objective, direction = "min"))
      search.result$best$objective <- as.numeric(payload.result$objective[1L])
    if (isTRUE(getOption("np.developer.native.nomad.diagnostics", FALSE)) &&
        !is.null(search.result$best_payload))
      attr(search.result$best_payload, "native.nomad.diagnostics") <- search.result$native.diagnostics

    return(search.result)
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
    remin = isTRUE(opt.args$nomad.remin),
    nomad.opts = nomad.opts,
    source = source,
    reason = reason,
    progress_label = progress_label,
    start.lower = c(bw_start_bounds$lower, degree.search$lower),
    start.upper = c(bw_start_bounds$upper, degree.search$upper),
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
                                            bernstein.named,
                                            nomad.source = "explicit",
                                            nomad.auto.filled = character()) {
  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  if (identical(degree.select, "manual"))
    return(NULL)
  resolved <- .np_degree_resolve_auto_engine(
    search.engine = search.engine,
    degree.select = degree.select,
    ncon = ncon,
    source = nomad.source,
    auto.filled = nomad.auto.filled
  )
  search.engine <- .npregbw_nomad_controls(resolved$search.engine)
  degree.select <- resolved$degree.select

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
    candidates = bounds$candidates,
    lower = bounds$lower,
    upper = bounds$upper,
    grid.size = bounds$grid.size,
    singleton = bounds$singleton,
    fixed.degree = bounds$fixed.degree,
    baseline.degree = baseline.degree,
    start.degree = start.degree,
    start.user = !is.null(degree.start),
    basis = if (missing(basis) || is.null(basis)) "glp" else as.character(basis[1L]),
    nobs = as.integer(nobs[1L]),
    restarts = npValidateNonNegativeInteger(degree.restarts, "degree.restarts"),
    max.cycles = npValidatePositiveInteger(degree.max.cycles, "degree.max.cycles"),
    verify = npValidateScalarLogical(degree.verify, "degree.verify"),
    bernstein.basis = bern.auto,
    source = resolved$source,
    reason = resolved$reason
  )
}

.npregbw_attach_degree_search <- function(bws, search_result) {
  cumulative_count <- function(name, total.name) {
    if (!is.null(search_result[[total.name]])) {
      total <- suppressWarnings(as.numeric(search_result[[total.name]][1L]))
      if (is.finite(total))
        return(total)
    }

    trace <- search_result$trace
    trace.summable <- !(search_result$method %in% c("nomad", "nomad+powell"))
    if (!is.null(search_result$engine))
      trace.summable <- trace.summable && !(search_result$engine %in% c("nomad", "nomad+powell"))
    if (!isTRUE(trace.summable))
      return(NA_real_)
    if (!is.data.frame(trace) || !nrow(trace) || !(name %in% names(trace)))
      return(NA_real_)

    actual <- rep(TRUE, nrow(trace))
    if ("status" %in% names(trace))
      actual <- actual & !is.na(trace$status) & as.character(trace$status) == "ok"
    if ("cached" %in% names(trace)) {
      cached <- suppressWarnings(as.logical(trace$cached))
      cached[is.na(cached)] <- FALSE
      actual <- actual & !cached
    }

    values <- suppressWarnings(as.numeric(trace[[name]]))
    values <- values[actual & is.finite(values)]
    if (!length(values))
      return(NA_real_)
    sum(values)
  }

  metadata <- .np_degree_search_metadata(search_result, default_direction = "min")

  if (!is.null(search_result$nomad.time))
    bws$nomad.time <- as.numeric(search_result$nomad.time[1L])
  if (!is.null(search_result$powell.time))
    bws$powell.time <- as.numeric(search_result$powell.time[1L])
  if (!is.null(search_result$optim.time) && is.finite(search_result$optim.time))
    bws$total.time <- as.numeric(search_result$optim.time[1L])
  bws <- .np_attach_nomad_restart_summary(bws, search_result)
  total.num.feval <- cumulative_count("num.feval", "num.feval.total")
  total.num.feval.fast <- cumulative_count("num.feval.fast", "num.feval.fast.total")
  if (is.finite(total.num.feval))
    bws$num.feval <- total.num.feval
  if (is.finite(total.num.feval.fast))
    bws$num.feval.fast <- total.num.feval.fast
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
           nomad.remin = FALSE,
           powell.remin,
           bwsolver = c("powell", "mads", "mads+powell"),
           scale.init.categorical.sample,
           scale.factor.search.lower = NULL,
           small,
           tol,
           transform.bounds = FALSE,
           ukertype,
           ...,
           nomad.opts = list()){

    xdat <- toFrame(xdat)
    yname <- deparse(substitute(ydat))
    nomad.opts <- .np_nomad_normalize_user_opts(nomad.opts, "npregbw")
    lp.dot.args <- list(...)
    if (length(nomad.opts))
      .np_nomad_native_reject_unsupported_options_for_route(
        opts = nomad.opts,
        route = "native npreg NOMAD route",
        nomad = nomad,
        degree.select = degree.select,
        search.engine = search.engine,
        bwsolver = bwsolver
      )
    if ("remin" %in% names(lp.dot.args)) {
      legacy.remin <- npValidateScalarLogical(lp.dot.args$remin, "remin")
      warning("npregbw: argument 'remin' is deprecated; use 'powell.remin' and 'nomad.remin'",
              call. = FALSE)
      if (missing(powell.remin))
        powell.remin <- legacy.remin
      if (missing(nomad.remin))
        nomad.remin <- legacy.remin
      lp.dot.args$remin <- NULL
    }
    npRejectLegacyLpArgs(names(lp.dot.args), where = "npregbw")
    random.seed.value <- .np_degree_extract_random_seed(lp.dot.args)
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

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
      if (sum(untangle(xdat)$icon) == 0L)
        stop("nomad=TRUE requires at least one continuous predictor for degree search",
             call. = FALSE)
      if ("degree" %in% mc.names)
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$bwtype, c("fixed", "generalized_nn", "adaptive_nn")))[1L] %in%
              c("fixed", "generalized_nn", "adaptive_nn")))
        stop("nomad=TRUE requires bwtype='fixed', 'generalized_nn', or 'adaptive_nn'")
      if ("degree.select" %in% mc.names &&
          identical(as.character(match.arg(nomad.shortcut$values$degree.select, c("manual", "coordinate", "exhaustive")))[1L], "manual"))
        stop("nomad=TRUE requires automatic degree search; use degree.select='coordinate' or 'exhaustive'")
      if (!identical(nomad.shortcut$metadata$source, "auto") &&
          "search.engine" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$search.engine, c("nomad+powell", "cell", "nomad")))[1L] %in%
              c("nomad", "nomad+powell")))
        stop("nomad=TRUE requires search.engine='nomad' or 'nomad+powell'")
      if ("degree.verify" %in% mc.names &&
          isTRUE(npValidateScalarLogical(nomad.shortcut$values$degree.verify, "degree.verify")))
        stop("nomad=TRUE currently requires degree.verify=FALSE")
    }

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

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("nmulti", "nomad.remin", "powell.remin", "bwsolver", "itmax", "ftol", "tol",
               "nomad.opts",
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
    search.mc.names <- mc.names
    regtype.value <- if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else "lc"
    bernstein.value <- if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else TRUE
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    search.engine.value <- if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell"
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
      bernstein.named = isTRUE(nomad.shortcut$enabled) || ("bernstein.basis" %in% search.mc.names),
      nomad.source = nomad.shortcut$metadata$source,
      nomad.auto.filled = nomad.shortcut$metadata$auto.filled
    )
    if (!is.null(degree.search) &&
        "bwsolver" %in% search.mc.names &&
        npBwsolverUsesMads(bwsolver)) {
      stop("bwsolver is for fixed-degree bandwidth searches; use search.engine for automatic degree search")
    }
    nomad.inner.named <- "nomad.nmulti" %in% search.mc.names
    nomad.inner.nmulti <- if (nomad.inner.named) {
      npValidateNonNegativeInteger(nomad.nmulti, "nomad.nmulti")
    } else {
      0L
    }
    if (nomad.inner.named &&
        (is.null(degree.search) || !(degree.search$engine %in% c("nomad", "nomad+powell")))) {
      stop("nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad+powell'")
    }
    if (!is.null(degree.search) && is.null(reg.args$degree)) {
      reg.args$degree <- npSetupGlpDegree(
        regtype = regtype.value,
        degree = NULL,
        ncon = ncon,
        degree.select = degree.select.value
      )
    }

    if (!is.null(degree.search)) {
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
          num.feval = if (!is.null(cell.bws$num.feval)) as.numeric(cell.bws$num.feval[1L]) else NA_real_,
          num.feval.fast = if (!is.null(cell.bws$num.feval.fast)) as.numeric(cell.bws$num.feval.fast[1L]) else NA_real_,
          nn.cache = cell.bws$nn.cache
        )
      }

      if (isTRUE(degree.search$singleton)) {
        search.result <- .np_degree_singleton_search_result(
          degree.search = degree.search,
          eval_result = eval_fun(degree.search$fixed.degree),
          direction = "min",
          objective_name = "fval"
        )
      } else if (identical(degree.search$engine, "cell")) {
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
          source = degree.search$source,
          reason = degree.search$reason,
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
          random.seed = random.seed.value,
          nomad.opts = if (is.null(opt.args$nomad.opts)) list() else opt.args$nomad.opts,
          source = degree.search$source,
          reason = degree.search$reason,
          progress_label = .np_degree_search_label(degree.search$engine, degree.search$source)
        )
      }
      tbw <- .npregbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
      tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)
      mc <- match.call(expand.dots = FALSE)
      environment(mc) <- parent.frame()
      tbw$call <- mc
      return(tbw)
    }

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
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
    
  }
