npudistbw <- function(...){
  mc <- match.call(expand.dots = FALSE)
  npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npudistbw")
  npRejectUnsupportedLpDegreeSearchArgs(names(mc$...), where = "npudistbw")
  target <- .np_bw_dispatch_target(dots = mc$...,
                                   data_arg_names = c("dat", "gdat"),
                                   eval_env = parent.frame())
  UseMethod("npudistbw", target)
}

.npudistbw_method_name <- function(bws, where = "npudistbw") {
  method <- bws[["method"]]
  method <- tryCatch(as.character(method)[1L], error = function(e) NA_character_)
  if (is.na(method) || !nzchar(method))
    stop(where, " requires a valid bandwidth method")
  switch(method,
         cv.cdf = method,
         "normal-reference" = method,
         stop(where, " does not support bwmethod = '", method, "'"))
}

.npudistbw_method_code <- function(bws, where = "npudistbw") {
  switch(.npudistbw_method_name(bws, where = where),
         cv.cdf = DBWM_CVLS,
         "normal-reference" = NA_integer_)
}

.npudistbw_prepare_cdf_grid <- function(dat,
                                        bws,
                                        gdat = NULL,
                                        do.full.integral = FALSE,
                                        ngrid = 100L) {
  if (!is.null(gdat)) {
    gdat <- toFrame(gdat)
    if (any(is.na(gdat)))
      stop("na's not allowed to be present in cdf gdata")
    gdat.matrix <- toMatrix(gdat)
    return(list(
      guno = gdat.matrix[, bws$iuno, drop = FALSE],
      gord = gdat.matrix[, bws$iord, drop = FALSE],
      gcon = gdat.matrix[, bws$icon, drop = FALSE],
      cdf_on_train = FALSE,
      nog = nrow(gdat.matrix)
    ))
  }

  if (isTRUE(do.full.integral)) {
    return(list(
      guno = data.frame(),
      gord = data.frame(),
      gcon = data.frame(),
      cdf_on_train = TRUE,
      nog = 0L
    ))
  }

  nog <- as.integer(ngrid)
  probs <- seq(0, 1, length.out = nog)
  ev <- dat[rep(1L, nog), , drop = FALSE]
  for (i in seq_len(ncol(ev)))
    ev[, i] <- cast(uocquantile(dat[, i], probs), dat[, i])
  ev <- toMatrix(ev)
  list(
    guno = ev[, bws$iuno, drop = FALSE],
    gord = ev[, bws$iord, drop = FALSE],
    gcon = ev[, bws$icon, drop = FALSE],
    cdf_on_train = FALSE,
    nog = nog
  )
}

npudistbw.formula <-
  function(formula, data, subset, na.action, call, gdata = NULL, ...){
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

    if (attr(attr(mf, "terms"), "response") != 0)
      stop("invalid distribution formula")

    dat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    has.gval <- !is.null(gdata)
    if (has.gval) {
      gmf <- match.call(expand.dots = FALSE)
      gm <- match(c("formula", "gdata"), names(gmf), nomatch = 0)
      gmf <- gmf[c(1,gm)]

      gmf[[1]] <- as.name("model.frame")
      names(gmf)[3] <- "data"
      gmf.args <- as.list(gmf[-1L])
      gmf <- do.call(stats::model.frame, gmf.args, envir = parent.frame())

      gdat <- gmf[, attr(attr(gmf, "terms"),"term.labels"), drop = FALSE]

    }

    dots <- list(...)
    seed.args <- c(list(dat = dat), if (has.gval) list(gdat = gdat) else list(), dots)
    tbw <- do.call(npudistbw, seed.args)
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$terms <- attr(mf,"terms")
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw
  }


npudistbw.NULL <-
  function(dat = stop("invoked without input data 'dat'"),
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npudistbw()")
    dots <- list(...)
    .np_nomad_native_reject_unsupported_options_from_dots(
      dots,
      "native npudist NOMAD route"
    )
    if (.npRmpi_autodispatch_active())
      {
        dat.preflight <- toFrame(dat)
        if (anyNA(dat.preflight) && !any(stats::complete.cases(dat.preflight)))
          stop("Data has no rows without NAs")
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_expand_dots_call(match.call(expand.dots = FALSE)),
        parent.frame()))
      }

    t.names <- NULL
    if(!is.data.frame(dat) && !is.matrix(dat))
      t.names <- paste(deparse(substitute(dat)), collapse = "")

    dat = toFrame(dat)

    if(!is.null(t.names))
      names(dat) <- t.names

    if (anyNA(dat) && !any(stats::complete.cases(dat)))
      stop("Data has no rows without NAs")

    bws = double(dim(dat)[2])

    tbw <- do.call(npudistbw.default, c(list(dat = dat, bws = bws), dots))

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw
  }

.npudistbw_nomad_native_target <- function(template, bwsolver) {
  method <- if (!is.null(template$method) && length(template$method)) {
    as.character(template$method[1L])
  } else {
    "cv.cdf"
  }
  bwtype <- if (!is.null(template$type) && length(template$type)) {
    as.character(template$type[1L])
  } else {
    ""
  }

  method %in% c("cv.cdf") &&
    bwtype %in% c("fixed", "generalized_nn", "adaptive_nn") &&
    bwsolver %in% c("mads", "mads+powell")
}

.npudistbw_nomad_native_require_crs <- function() {
  if (!requireNamespace("crs", quietly = TRUE))
    stop("native npudist NOMAD route requires crs >= 0.15-44", call. = FALSE)
  if (utils::packageVersion("crs") < "0.15.44")
    stop("native npudist NOMAD route requires crs >= 0.15-44", call. = FALSE)
  invisible(TRUE)
}

.npudistbw_nomad_native_option_vectors <- .npudensbw_nomad_native_option_vectors

.npudistbw_nomad_native_prepare_args <- function(dat,
                                                 bws,
                                                 gdat = NULL,
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
                                                 transform.bounds = FALSE,
                                                 do.full.integral = FALSE,
                                                 ngrid = 100L,
                                                 memfac = 500.0) {
  invalid.penalty <- match.arg(invalid.penalty)
  dat <- toFrame(dat)
  dat.matrix <- toMatrix(dat)

  duno <- dat.matrix[, bws$iuno, drop = FALSE]
  dcon <- dat.matrix[, bws$icon, drop = FALSE]
  dord <- dat.matrix[, bws$iord, drop = FALSE]
  mysd <- EssDee(dcon)
  nrow <- dim(dat.matrix)[1L]
  nconfac <- nrow^(-1.0 / (2.0 * bws$ckerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$ckerorder + bws$ncon))
  sfloor <- npResolveScaleFactorLowerBound(
    if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
  )
  cont.start <- npContinuousSearchStartControls(
    scale.factor.init.lower,
    scale.factor.init.upper,
    scale.factor.init,
    sfloor,
    where = "npudistbw"
  )

  cdf.grid <- .npudistbw_prepare_cdf_grid(
    dat = dat,
    bws = bws,
    gdat = gdat,
    do.full.integral = do.full.integral,
    ngrid = ngrid
  )
  guno <- cdf.grid$guno
  gord <- cdf.grid$gord
  gcon <- cdf.grid$gcon
  cdf_on_train <- cdf.grid$cdf_on_train
  nog <- cdf.grid$nog

  myopti <- list(
    num_obs_train = nrow,
    num_obs_eval = nog,
    iMultistart = IMULTI_TRUE,
    iNum_Multistart = 1L,
    int_use_starting_values = USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = itmax,
    int_RESTART_FROM_MIN = RE_MIN_FALSE,
    int_MINIMIZE_IO = IO_MIN_TRUE,
    bwmethod = .npudistbw_method_code(bws, where = "npudistbw"),
    ckerneval = switch(bws$ckertype,
      gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS),
    ukerneval = switch(bws$ukertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    okerneval = switch(bws$okertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY),
    cdf_on_train = cdf_on_train,
    nuno = dim(duno)[2L],
    nord = dim(dord)[2L],
    ncon = dim(dcon)[2L],
    int_do_tree = npDoTreeOrCategoricalCompress(
      ncon = dim(dcon)[2L],
      ncat = dim(duno)[2L] + dim(dord)[2L],
      bws = bws),
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
    memfac = memfac,
    scale.factor.lower.bound = sfloor
  )

  cker.bounds <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])
  list(
    duno = as.double(duno),
    dord = as.double(dord),
    dcon = as.double(dcon),
    guno = as.double(guno),
    gord = as.double(gord),
    gcon = as.double(gcon),
    mysd = as.double(mysd),
    myopti = as.integer(myopti),
    myoptd = as.double(myoptd),
    penalty_mode = as.integer(if (invalid.penalty == "baseline") 1L else 0L),
    penalty_multiplier = as.double(penalty.multiplier),
    ckerlb = as.double(cker.bounds$lb),
    ckerub = as.double(cker.bounds$ub)
  )
}

npNomadNativeSearchDistribution <- function(prep,
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
    "C_np_distribution_nomad_native_search",
    as.double(prep$duno),
    as.double(prep$dord),
    as.double(prep$dcon),
    as.double(prep$guno),
    as.double(prep$gord),
    as.double(prep$gcon),
    as.double(prep$mysd),
    as.integer(prep$myopti),
    as.double(prep$myoptd),
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
    PACKAGE = "npRmpi"
  ), capture.output = TRUE)
  .np_nomad_native_call_value(native.call)
}

.npudistbw_run_mads <- function(dat,
                                bws,
                                gdat = NULL,
                                opt.args,
                                bwsolver = c("mads", "mads+powell")) {
  bwsolver <- npValidateBwsolver(bwsolver)
  template <- bws
  if (!(template$type %in% c("fixed", "generalized_nn", "adaptive_nn")))
    stop("bwsolver='mads' requires bwtype='fixed', 'generalized_nn', or 'adaptive_nn'")

  setup <- .npregbw_nomad_bw_setup(
    xdat = dat,
    template = template,
    allow.extended.nn = TRUE,
    evaldat = if (is.null(gdat)) dat else gdat
  )
  bounds <- .npregbw_nomad_bw_bounds(template = template, setup = setup)
  point.start <- if (all(template$bw == 0)) NULL else .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup)
  x0 <- .npregbw_nomad_complete_bw_start_point(point = point.start, bounds = bounds, setup = setup)
  opt.value <- function(name, default) {
    if (is.null(opt.args[[name]])) default else opt.args[[name]]
  }
  mads.num.feval.total <- 0
  mads.num.feval.fast.total <- 0

  eval_fun <- function(point) {
    bw_vec <- .npregbw_nomad_point_to_bw(point, template = template, setup = setup)
    tbw <- bws
    tbw$bw <- bw_vec
    out <- npudistbw.dbandwidth(
      dat = dat,
      bws = tbw,
      gdat = gdat,
      bandwidth.compute = TRUE,
      nmulti = 1L,
      powell.remin = FALSE,
      invalid.penalty = opt.value("invalid.penalty", "baseline"),
      penalty.multiplier = opt.value("penalty.multiplier", 10),
      scale.factor.search.lower = opt.value("scale.factor.search.lower", NULL),
      do.full.integral = opt.value("do.full.integral", FALSE),
      ngrid = opt.value("ngrid", 100L),
      memfac = opt.value("memfac", 500.0),
      bwsolver = "powell",
      eval.only = TRUE
    )
    mads.num.feval.total <<- mads.num.feval.total + as.numeric(out$num.feval[1L])
    mads.num.feval.fast.total <<- mads.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])
    list(objective = out$fval, degree = integer(0L), num.feval = out$num.feval)
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    bw_vec <- .npregbw_nomad_point_to_bw(point, template = template, setup = setup)
    final.tbw <- bws
    final.tbw$bw <- bw_vec
    final.tbw$fval <- final.tbw$ifval <- as.numeric(best_record$objective)
    final.tbw$num.feval <- as.numeric(mads.num.feval.total)
    final.tbw$num.feval.fast <- as.numeric(mads.num.feval.fast.total)
    final.tbw$fval.history <- as.numeric(best_record$objective)
    final.tbw$eval.history <- if (!is.null(solution$bbe)) rep(1, max(1L, as.integer(solution$bbe))) else 1
    final.tbw$invalid.history <- 0
    final.tbw$timing <- NA_real_
    direct.payload <- npudistbw.dbandwidth(dat = dat, bws = final.tbw, gdat = gdat, bandwidth.compute = FALSE)
    direct.objective <- as.numeric(best_record$objective)
    powell.elapsed <- NA_real_

    if (identical(bwsolver, "mads+powell")) {
      hot.start <- proc.time()[3L]
      hot.payload <- npudistbw.dbandwidth(
        dat = dat,
        bws = final.tbw,
        gdat = gdat,
        bandwidth.compute = TRUE,
        nmulti = 1L,
        powell.remin = isTRUE(opt.args$powell.remin),
        invalid.penalty = opt.value("invalid.penalty", "baseline"),
        penalty.multiplier = opt.value("penalty.multiplier", 10),
        scale.factor.search.lower = opt.value("scale.factor.search.lower", NULL),
        do.full.integral = opt.value("do.full.integral", FALSE),
        ngrid = opt.value("ngrid", 100L),
        memfac = opt.value("memfac", 500.0),
        bwsolver = "powell"
      )
      powell.elapsed <- proc.time()[3L] - hot.start
      hot.payload$num.feval <- as.numeric(direct.payload$num.feval[1L]) + as.numeric(hot.payload$num.feval[1L])
      hot.payload$num.feval.fast <- as.numeric(direct.payload$num.feval.fast[1L]) + as.numeric(hot.payload$num.feval.fast[1L])
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min"))
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  native.start.bounds <- .np_nomad_bw_restart_start_bounds(
    bounds = bounds,
    setup = setup,
    opt.value = opt.value,
    where = "npudistbw"
  )
  if (is.null(point.start)) {
    x0 <- .npregbw_nomad_complete_bw_start_point(
      point = NULL,
      bounds = bounds,
      setup = setup,
      initial = native.start.bounds$initial,
      where = "npudistbw"
    )
  }

  if (.npudistbw_nomad_native_target(template, bwsolver)) {
    .npudistbw_nomad_native_require_crs()
    native.nmulti <- npValidateNmulti(opt.value("nmulti", npDefaultNmulti(dim(toFrame(dat))[2L])))
    native.inner.nmulti <- npValidateNonNegativeInteger(
      opt.value("mads.nmulti", opt.value("nomad.nmulti", 0L)),
      "nomad.nmulti"
    )
    native.inner.nmulti <- as.integer(native.inner.nmulti[1L])
    if (isTRUE(opt.args$nomad.remin))
      stop("native npudist NOMAD route does not support NOMAD remin", call. = FALSE)

    native.random.seed <- opt.value("random.seed", 42L)
    native.nomad.opts <- .np_nomad_prepare_solver_opts(
      random.seed = native.random.seed,
      nomad.opts = opt.value("nomad.opts", list()),
      coordinate.roles = .np_nomad_coordinate_roles(bounds),
      expected.length = length(bounds$lower),
      geometry.policy = "generate-central",
      where = "npudistbw native NOMAD source geometry"
    )
    native.option.vectors <- .npudistbw_nomad_native_option_vectors(native.nomad.opts)
    native.start.matrix <- .np_nomad_build_starts(
      x0 = x0,
      bbin = bounds$bbin,
      lb = bounds$lower,
      ub = bounds$upper,
      nmulti = native.nmulti,
      random.seed = native.random.seed,
      degree_spec = NULL,
      start.lower = native.start.bounds$lower,
      start.upper = native.start.bounds$upper
    )
    native.prep <- .npudistbw_nomad_native_prepare_args(
      dat = dat,
      bws = bws,
      gdat = gdat,
      invalid.penalty = opt.value("invalid.penalty", "baseline"),
      penalty.multiplier = opt.value("penalty.multiplier", 10),
      itmax = opt.value("itmax", 10000L),
      ftol = opt.value("ftol", 1.490116e-07),
      tol = opt.value("tol", 1.490116e-04),
      small = opt.value("small", 1.490116e-05),
      scale.factor.init.lower = opt.value("scale.factor.init.lower", 0.1),
      scale.factor.init.upper = opt.value("scale.factor.init.upper", 2.0),
      scale.factor.init = opt.value("scale.factor.init", 0.5),
      lbd.init = opt.value("lbd.init", 0.1),
      hbd.init = opt.value("hbd.init", 0.9),
      dfac.init = opt.value("dfac.init", 0.375),
      scale.factor.search.lower = opt.value("scale.factor.search.lower", NULL),
      do.full.integral = opt.value("do.full.integral", FALSE),
      ngrid = opt.value("ngrid", 100L),
      memfac = opt.value("memfac", 500.0)
    )

    native.results <- vector("list", nrow(native.start.matrix))
    native.best.index <- NA_integer_
    native.best.objective <- Inf
    native.nomad.elapsed <- 0
    native.num.feval.total <- 0
    native.num.feval.fast.total <- 0
    native.num.feval.guarded.total <- 0
    for (i in seq_len(nrow(native.start.matrix))) {
      native.start <- proc.time()[3L]
      native.i <- npNomadNativeSearchDistribution(
        prep = native.prep,
        x0 = as.numeric(native.start.matrix[i, ]),
        bbin = bounds$bbin,
        lb = bounds$lower,
        ub = bounds$upper,
        max.eval = 0L,
        random.seed = native.random.seed,
        inner.start.count = native.inner.nmulti,
        option.names = native.option.vectors$names,
        option.values = native.option.vectors$values
      )
      native.elapsed <- proc.time()[3L] - native.start
      native.nomad.elapsed <- native.nomad.elapsed + native.elapsed
      if (!identical(as.integer(native.i$status[1L]), 0L) ||
          !identical(as.integer(native.i$result_status[1L]), 0L)) {
        stop(sprintf(
          "native npudist NOMAD route failed (status=%s, result_status=%s): %s",
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
      native.num.feval.guarded.total <- native.num.feval.guarded.total + as.numeric(native.i$total_num.feval.guarded[1L])
      if (is.finite(objective.i) && objective.i < native.best.objective) {
        native.best.objective <- objective.i
        native.best.index <- i
      }
    }
    if (!is.finite(native.best.index))
      stop("native npudist NOMAD route did not return a finite solution", call. = FALSE)

    native.best <- native.results[[native.best.index]]
    native.handoff.point <- as.numeric(native.best$best_point)
    if (any(!is.finite(native.handoff.point)))
      stop("native npudist NOMAD route did not return a finite best point", call. = FALSE)
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
      num.feval.guarded = as.numeric(native.best$native$best_num.feval.guarded[1L])
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
      restart.results = native.results,
      restart.starts = lapply(seq_len(nrow(native.start.matrix)), function(i) as.numeric(native.start.matrix[i, ])),
      best.restart = native.best.index,
      nomad.time = native.nomad.elapsed,
      powell.time = payload.result$powell.time,
      optim.time = native.nomad.elapsed + as.numeric(payload.result$powell.time[1L]),
      num.feval.total = native.num.feval.total,
      num.feval.fast.total = native.num.feval.fast.total,
      num.feval.guarded.total = native.num.feval.guarded.total,
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
      nmulti = opt.value("nmulti", npDefaultNmulti(dim(toFrame(dat))[2L])),
      nomad.inner.nmulti = opt.value("mads.nmulti", opt.value("nomad.nmulti", 0L)),
      random.seed = opt.value("random.seed", 42L),
      handoff_before_build = identical(bwsolver, "mads+powell"),
      remin = isTRUE(opt.args$nomad.remin),
      nomad.opts = opt.value("nomad.opts", list()),
      start.lower = native.start.bounds$lower,
      start.upper = native.start.bounds$upper
    )
  }
  search.result$method <- bwsolver
  out <- search.result$best_payload
  out$bwsolver <- bwsolver
  out$search.engine <- bwsolver
  out$nomad.time <- as.numeric(search.result$nomad.time[1L])
  out$powell.time <- as.numeric(search.result$powell.time[1L])
  out$total.time <- as.numeric(search.result$optim.time[1L])
  .np_attach_nomad_restart_summary(out, search.result)
}

npudistbw.dbandwidth <- 
  function(dat = stop("invoked without input data 'dat'"),
           bws,
           gdat = NULL,
           bandwidth.compute = TRUE,
           cfac.dir = 2.5*(3.0-sqrt(5)),
           scale.factor.init = 0.5,
           dfac.dir = 0.25*(3.0-sqrt(5)),
           dfac.init = 0.375,
           dfc.dir = 3,
           do.full.integral = FALSE,
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
           ngrid = 100,
           nmulti,
           penalty.multiplier = 10,
           powell.remin = TRUE,
           bwsolver = c("powell", "mads", "mads+powell"),
           scale.init.categorical.sample = FALSE,
           scale.factor.search.lower = NULL,
           small = 1.490116e-05,
           tol = 1.490116e-04,
           transform.bounds = FALSE,
           eval.only = FALSE,
           ...){
    dot.args <- list(...)
    elapsed.start <- proc.time()[3]
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    bwsolver <- npValidateBwsolver(bwsolver)
    eval.only <- npValidateScalarLogical(eval.only, "eval.only")
    remin <- npValidateScalarLogical(powell.remin, "powell.remin")
    do.full.integral <- npValidateScalarLogical(do.full.integral, "do.full.integral")
    scale.init.categorical.sample <-
      npValidateScalarLogical(scale.init.categorical.sample, "scale.init.categorical.sample")
    transform.bounds <- npValidateScalarLogical(transform.bounds, "transform.bounds")
    itmax <- npValidatePositiveInteger(itmax, "itmax")
    ngrid <- npValidatePositiveInteger(ngrid, "ngrid")
    ftol <- npValidatePositiveFiniteNumeric(ftol, "ftol")
    tol <- npValidatePositiveFiniteNumeric(tol, "tol")
    small <- npValidatePositiveFiniteNumeric(small, "small")
    memfac <- npValidatePositiveFiniteNumeric(memfac, "memfac")
    penalty.multiplier <- npValidatePositiveFiniteNumeric(penalty.multiplier, "penalty.multiplier")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    if (!missing(nmulti))
      nmulti <- npValidateNmulti(nmulti)
    .npRmpi_require_active_slave_pool(where = "npudistbw()")
    if (.npRmpi_autodispatch_active()) {
      if (bandwidth.compute)
        .npudistbw_method_name(bws, where = "npudistbw")
      dat.preflight <- toFrame(dat)
      if (anyNA(dat.preflight) && !any(stats::complete.cases(dat.preflight)))
        stop("Data has no rows without NAs")
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_expand_dots_call(match.call(expand.dots = FALSE)),
        parent.frame()))
    }

    dat = toFrame(dat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(dim(dat)[2])
    }
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)

    if (length(bws$bw) != dim(dat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of",
           "'dat'"))

    if ((any(bws$icon) &&
         !all(vapply(as.data.frame(dat[, bws$icon]), inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iord) &&
         !all(vapply(as.data.frame(dat[, bws$iord]), inherits, logical(1), "ordered"))) ||
        (any(bws$iuno) &&
         !all(vapply(as.data.frame(dat[, bws$iuno]), inherits, logical(1), "factor"))))
      stop(paste("supplied bandwidths do not match", "'dat'", "in type"))

    npValidateExtendedNnContinuousBandwidth(bws, where = "npudistbw")

    if(any(bws$iuno))
      stop("distribution bandwidth selection does not support unordered data types")

    dat <- na.omit(dat)
    rows.omit <- unclass(na.action(dat))
    if (nrow(dat) == 0L)
      stop("Data has no rows without NAs")

    nrow = dim(dat)[1]
    ncol = dim(dat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    odat <- dat
    dat.frame <- dat
    dat = toMatrix(dat)

    duno = dat[, bws$iuno, drop = FALSE]
    dcon = dat[, bws$icon, drop = FALSE]
    dord = dat[, bws$iord, drop = FALSE]

    tbw <- bws

    mysd <- EssDee(dcon)
    nconfac <- nrow^(-1.0/(1.0+bws$ckerorder))
    ncatfac <- nrow^(-2.0/(1.0+bws$ckerorder))

    cdf.grid <- .npudistbw_prepare_cdf_grid(
      dat = odat,
      bws = bws,
      gdat = gdat,
      do.full.integral = do.full.integral,
      ngrid = ngrid
    )
    guno = cdf.grid$guno
    gord = cdf.grid$gord
    gcon = cdf.grid$gcon
    cdf_on_train = cdf.grid$cdf_on_train
    nog = cdf.grid$nog

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

    if (bandwidth.compute && !eval.only && npBwsolverUsesMads(bwsolver)) {
      return(.npudistbw_run_mads(
        dat = dat.frame,
        bws = bws,
        gdat = gdat,
        opt.args = list(
          nmulti = nmulti,
          mads.nmulti = dot.args$mads.nmulti,
          nomad.nmulti = dot.args$nomad.nmulti,
          nomad.remin = if (is.null(dot.args$nomad.remin)) FALSE else dot.args$nomad.remin,
          powell.remin = remin,
          bwsolver = bwsolver,
          itmax = itmax,
          ftol = ftol,
          tol = tol,
          small = small,
          invalid.penalty = invalid.penalty,
          penalty.multiplier = penalty.multiplier,
          scale.factor.init.lower = scale.factor.init.lower,
          scale.factor.init.upper = scale.factor.init.upper,
          scale.factor.init = scale.factor.init,
          lbd.init = lbd.init,
          hbd.init = hbd.init,
          dfac.init = dfac.init,
          scale.init.categorical.sample = scale.init.categorical.sample,
          transform.bounds = transform.bounds,
          scale.factor.search.lower = scale.factor.search.lower,
          do.full.integral = do.full.integral,
          ngrid = ngrid,
          memfac = memfac,
          nomad.opts = dot.args$nomad.opts
        ),
        bwsolver = bwsolver
      ))
    }

    if (bandwidth.compute){
      cont.start <- npContinuousSearchStartControls(
        scale.factor.init.lower,
        scale.factor.init.upper,
        scale.factor.init,
        scale.factor.search.lower,
        where = "npudistbw"
      )
      myopti = list(num_obs_train = dim(dat)[1],
        num_obs_eval = nog,
        iMultistart = IMULTI_TRUE,
        iNum_Multistart = nmulti,
        int_use_starting_values = (if (all(bws$bw==0)) USE_START_NO else USE_START_YES),
        int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
        BANDWIDTH_den_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=(if (remin) RE_MIN_TRUE else RE_MIN_FALSE),
        int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
        bwmethod = .npudistbw_method_code(bws, where = "npudistbw"),
        ckerneval = switch(bws$ckertype,
          gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
          uniform = CKER_UNI,
          "truncated gaussian" = CKER_TGAUSS),
        ukerneval = switch(bws$ukertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        okerneval = switch(bws$okertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_NLR,
        "racineliyan" = OKER_RLY),
        cdf_on_train = cdf_on_train,
        nuno = dim(duno)[2],
        nord = dim(dord)[2],
        ncon = dim(dcon)[2],
        int_do_tree = npDoTreeOrCategoricalCompress(
          ncon = dim(dcon)[2],
          ncat = dim(duno)[2] + dim(dord)[2],
          bws = bws),
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
        nconfac = nconfac, ncatfac = ncatfac, memfac = memfac,
        scale.factor.lower.bound = scale.factor.search.lower)
      cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

      if (bws$method != "normal-reference"){
        if (isTRUE(eval.only)) {
          myout <-
            .Call("C_np_distribution_bw_eval",
                  as.double(duno), as.double(dord), as.double(dcon),
                  as.double(guno), as.double(gord), as.double(gcon), as.double(mysd),
                  as.integer(myopti), as.double(myoptd),
                  as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
                  as.integer(nmulti),
                  as.integer(penalty_mode),
                  as.double(penalty.multiplier),
                  as.double(cker.bounds.c$lb),
                  as.double(cker.bounds.c$ub),
                  PACKAGE="npRmpi")
        } else {
          myout <-
            .Call("C_np_distribution_bw",
                as.double(duno), as.double(dord), as.double(dcon),
                as.double(guno), as.double(gord), as.double(gcon), as.double(mysd),
                as.integer(myopti), as.double(myoptd),
                as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
                as.integer(nmulti),
                as.integer(penalty_mode),
                as.double(penalty.multiplier),
                as.double(cker.bounds.c$lb),
                as.double(cker.bounds.c$ub),
                PACKAGE="npRmpi")
        }
        total.time <- proc.time()[3] - elapsed.start
      } else {
        nbw = double(ncol)
        gbw = bws$ncon
        if (gbw > 0){
          gbw_idx <- seq_len(gbw)
          nbw[gbw_idx] = 1.587
          if(!bws$scaling)
            nbw[gbw_idx]=nbw[gbw_idx]*mysd*nconfac
        }
        myout= list( bw = nbw, fval = c(NA,NA) )
        total.time <- NA
      }

      rorder = numeric(ncol)
      ord_idx <- seq_len(ncol)
      rorder[c(ord_idx[bws$icon], ord_idx[bws$iuno], ord_idx[bws$iord])] <- ord_idx

      tbw$bw <- myout$bw[rorder]

      tbw$fval = myout$fval[1]
      tbw$ifval = myout$fval[2]
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
      tbw$nn.cache <- .np_nn_cache_stats(myout$nn.cache)
      tbw$fval.history <- myout$fval.history
      tbw$eval.history <- myout$eval.history
      tbw$invalid.history <- myout$invalid.history
      tbw$timing <- myout$timing
      tbw$total.time <- total.time
    }

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

    tbw <- dbandwidth(bw = tbw$bw,
                      bwmethod = tbw$method,
                      bwscaling = tbw$scaling,
                      bwtype = tbw$type,
                      ckertype = tbw$ckertype,
                      ckerorder = tbw$ckerorder,
                      ckerbound = tbw$ckerbound,
                      ckerlb = tbw$ckerlb,
                      ckerub = tbw$ckerub,
                      ukertype = c("aitchisonaitken"),
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
                      xnames = tbw$xnames,
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

npudistbw.default <-
  function(dat = stop("invoked without input data 'dat'"),
           bws,
           gdat,
           bandwidth.compute = TRUE,
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
           dfac.dir,
           dfac.init,
           dfc.dir,
           do.full.integral,
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
           ngrid,
           nmulti,
           okertype,
           penalty.multiplier,
           powell.remin,
           bwsolver = c("powell", "mads", "mads+powell"),
           scale.init.categorical.sample,
           scale.factor.search.lower = NULL,
           small,
           tol,
           transform.bounds,
           ...){
    .npRmpi_require_active_slave_pool(where = "npudistbw()")
    if (.npRmpi_autodispatch_active()) {
      dat.preflight <- toFrame(dat)
      if (anyNA(dat.preflight) && !any(stats::complete.cases(dat.preflight)))
        stop("Data has no rows without NAs")
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_expand_dots_call(match.call(expand.dots = FALSE)),
        parent.frame()))
    }

    t.names <- NULL
    if(!is.data.frame(dat) && !is.matrix(dat))
      t.names <- paste(deparse(substitute(dat)), collapse = "")

    dat <- toFrame(dat)

    if(!is.null(t.names))
      names(dat) <- t.names

    if (anyNA(dat) && !any(stats::complete.cases(dat)))
      stop("Data has no rows without NAs")

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    bw.args <- list(
      bw = bws,
      nobs = dim(dat)[1],
      xdati = untangle(dat),
      xnames = names(dat),
      bandwidth.compute = bandwidth.compute
    )
    if (!missing(bwmethod)) bw.args$bwmethod <- bwmethod
    if (!missing(bwscaling)) bw.args$bwscaling <- bwscaling
    if (!missing(bwtype)) bw.args$bwtype <- bwtype
    if (!missing(ckertype)) bw.args$ckertype <- ckertype
    if (!missing(ckerorder)) bw.args$ckerorder <- ckerorder
    if (!missing(ckerbound)) bw.args$ckerbound <- ckerbound
    if (!missing(ckerlb)) bw.args$ckerlb <- ckerlb
    if (!missing(ckerub)) bw.args$ckerub <- ckerub
    if (!missing(okertype)) bw.args$okertype <- okertype
    tbw <- do.call(dbandwidth, bw.args)


    ## next grab dummies for actual bandwidth selection and perform call

    dots <- list(...)
    dot.names <- names(dots)
    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("gdat","bandwidth.compute", "nmulti", "powell.remin", "bwsolver", "itmax",
               "do.full.integral", "ngrid", "ftol", "tol",
               "small", "lbc.dir", "dfc.dir", "cfac.dir", "initc.dir",
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.init.categorical.sample", "scale.factor.search.lower", "memfac",
               "transform.bounds",
               "invalid.penalty",
               "penalty.multiplier",
               "mads.nmulti", "nomad.nmulti", "nomad.remin", "nomad.opts")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    bwsel.args <- list(dat = dat, bws = tbw)
    if (any.m) {
      nms <- mc.names[m]
      bwsel.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    dotted.arg.names <- intersect(margs, dot.names)
    if (length(dotted.arg.names)) {
      bwsel.args[dotted.arg.names] <- dots[dotted.arg.names]
    }
    tbw <- .np_progress_select_bandwidth_enhanced(
      "Selecting distribution bandwidth",
      do.call(npudistbw.dbandwidth, bwsel.args)
    )

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }
