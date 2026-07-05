
npindexbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npindexbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat"),
                                     eval_env = parent.frame())
    UseMethod("npindexbw", target)
  }

npindexbw.formula <-
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

    mf[[1]] <- as.name("model.frame")

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

    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())

    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    tbw <- do.call(npindexbw, c(list(xdat = xdat, ydat = ydat), list(...)))

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

npindexbw.NULL <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npindexbw()")
    mc <- match.call(expand.dots = FALSE)
    dots <- list(...)
    dot.names <- names(dots)
    nomad.requested <- "nomad" %in% dot.names &&
      (npValidateNomadControl(dots$nomad, "nomad") %in% c("true", "auto"))
    degree.select.value <- if ("degree.select" %in% dot.names) {
      match.arg(as.character(dots$degree.select[[1L]]),
                c("manual", "coordinate", "exhaustive"))
    } else {
      "manual"
    }
    automatic.degree.search <- isTRUE(nomad.requested) ||
      !identical(degree.select.value, "manual")
    method.value <- if ("method" %in% dot.names) {
      match.arg(as.character(dots$method[[1L]]), c("ichimura", "kleinspady"))
    } else {
      "ichimura"
    }
    collective.degree.search <- isTRUE(automatic.degree.search) &&
      identical(method.value, "kleinspady")
    regtype.value <- if ("regtype" %in% dot.names) {
      match.arg(as.character(dots$regtype[[1L]]), c("lc", "ll", "lp"))
    } else if (isTRUE(nomad.requested)) {
      "lp"
    } else {
      "lc"
    }
    search.engine.value <- if ("search.engine" %in% dot.names) {
      match.arg(as.character(dots$search.engine[[1L]]), c("nomad+powell", "cell", "nomad"))
    } else if (isTRUE(nomad.requested)) {
      "nomad+powell"
    } else {
      "nomad+powell"
    }
    ichimura.lp.nomad.degree.search <- isTRUE(automatic.degree.search) &&
      identical(method.value, "ichimura") &&
      identical(regtype.value, "lp") &&
      search.engine.value %in% c("nomad", "nomad+powell")
    .np_nomad_validate_inner_multistart(
      call_names = names(mc),
      dot.args = dots,
      regtype = regtype.value,
      automatic.degree.search = automatic.degree.search,
      search.engine = search.engine.value
    )
    if (.npRmpi_autodispatch_active() &&
        (!isTRUE(automatic.degree.search) ||
           isTRUE(collective.degree.search) ||
           isTRUE(ichimura.lp.nomad.degree.search) ||
           .npRmpi_safe_int(mpi.comm.size(1L)) > 2L))
      return(.npRmpi_autodispatch_call(mc, parent.frame()))

    xdat <- toFrame(xdat)

    bws <- double(ncol(xdat)+1)

    tbw <- npindexbw.default(xdat = xdat,
                             ydat = ydat,
                             bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)

    tbw
  }

.npindex_resolve_spec <- function(source, where = "npindex") {
  if (!is.null(source$regtype.engine)) {
    return(list(
      regtype = if (is.null(source$regtype)) "lc" else as.character(source$regtype),
      basis = if (is.null(source$basis)) "glp" else as.character(source$basis),
      degree = if (is.null(source$degree)) integer(0) else as.integer(source$degree),
      bernstein.basis = isTRUE(source$bernstein.basis),
      regtype.engine = as.character(source$regtype.engine),
      basis.engine = if (is.null(source$basis.engine)) "glp" else as.character(source$basis.engine),
      degree.engine = if (is.null(source$degree.engine)) integer(0) else as.integer(source$degree.engine),
      bernstein.basis.engine = isTRUE(source$bernstein.basis.engine)
    ))
  }

  npCanonicalConditionalRegSpec(
    regtype = if (is.null(source$regtype)) "lc" else as.character(source$regtype),
    basis = if (is.null(source$basis)) "glp" else as.character(source$basis),
    degree = source$degree,
    bernstein.basis = isTRUE(source$bernstein.basis),
    ncon = 1L,
    where = where
  )
}

.npindex_objective_policy <- function(bws,
                                      spec,
                                      bandwidth.compute = TRUE,
                                      where = "npindexbw") {
  method <- as.character(bws$method[1L])
  bwtype <- as.character(bws$type[1L])
  public.regtype <- as.character(spec$regtype[1L])
  engine.regtype <- as.character(spec$regtype.engine[1L])
  objective.spec <- spec
  route <- if (identical(public.regtype, "lc")) {
    "lc"
  } else if (identical(public.regtype, "ll")) {
    "ll"
  } else {
    paste0("lp", paste(as.integer(spec$degree), collapse = ","))
  }

  canonical.degree <- if (identical(public.regtype, "lc")) {
    0L
  } else if (identical(public.regtype, "ll")) {
    1L
  } else {
    as.integer(spec$degree.engine)
  }

  canonical.degree0 <- length(canonical.degree) > 0L && all(canonical.degree == 0L)
  if (!(identical(method, "ichimura") || identical(method, "kleinspady"))) {
    stop(
      sprintf("%s received unsupported npindex method '%s'", where, method),
      call. = FALSE
    )
  }
  executor <- "npreg_loo"

  if (canonical.degree0) {
    objective.spec$regtype.engine <- "lc"
    objective.spec$basis.engine <- "glp"
    objective.spec$degree.engine <- 0L
    objective.spec$bernstein.basis.engine <- FALSE
  }

  list(
    method = method,
    bwtype = bwtype,
    bandwidth.compute = isTRUE(bandwidth.compute),
    public.regtype = public.regtype,
    public.degree = as.integer(spec$degree),
    route = route,
    canonical.regtype = if (identical(public.regtype, "lc") || identical(public.regtype, "ll")) "lp" else engine.regtype,
    canonical.degree = canonical.degree,
    executor = executor,
    objective.spec = objective.spec,
    support = "supported",
    where = where
  )
}

.npindex_nn_candidate_bandwidth <- function(h, bwtype, nobs) {
  if (identical(bwtype, "fixed")) {
    return(list(ok = is.finite(h) && (h > 0), value = as.double(h)))
  }

  if (!is.finite(h)) {
    return(list(ok = FALSE, value = NA_real_))
  }

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  hard.upper <- .Machine$integer.max / 2
  k <- .np_round_half_to_even(h)
  upper.ok <- (k <= upper) ||
    (npExtendedNnEnabled() && (k <= hard.upper) &&
       (as.character(bwtype)[1L] %in% c("generalized_nn", "adaptive_nn")))
  list(ok = (k >= lower) && upper.ok, value = as.double(k))
}

.npindexbw_h_start_controls <- function(scale.factor.init.lower = 0.1,
                                        scale.factor.init.upper = 2.0,
                                        scale.factor.init = 0.5,
                                        scale.factor.search.lower = 0,
                                        where = "npindexbw") {
  cont.start <- npContinuousSearchStartControls(
    scale.factor.init.lower,
    scale.factor.init.upper,
    scale.factor.init,
    scale.factor.search.lower,
    where = where
  )
  cont.start$scale.factor.search.lower <- as.double(scale.factor.search.lower)
  cont.start
}

.npindex_start_bandwidth_scale <- function(fit, nobs) {
  EssDee(fit) * nobs^(-1 / 5)
}

.npindex_default_start_bandwidth <- function(fit,
                                             bwtype,
                                             nobs,
                                             start.controls = .npindexbw_h_start_controls()) {
  if (identical(bwtype, "fixed")) {
    return(start.controls$scale.factor.init * .npindex_start_bandwidth_scale(fit = fit, nobs = nobs))
  }

  lower <- 2L
  max(lower, min(max(1L, as.integer(nobs) - 1L), .np_round_half_to_even(sqrt(nobs))))
}

.npindex_random_start_bandwidth <- function(fit,
                                            bwtype,
                                            nobs,
                                            start.controls = .npindexbw_h_start_controls()) {
  if (identical(bwtype, "fixed")) {
    return(runif(1, min = start.controls$scale.factor.init.lower, max = start.controls$scale.factor.init.upper) *
             .npindex_start_bandwidth_scale(fit = fit, nobs = nobs))
  }

  upper <- max(1L, as.integer(nobs) - 1L)
  runif(1, min = 2, max = max(2L, upper))
}

.npindex_ols_beta_tail <- function(ols.fit) {
  slopes <- as.double(coef(ols.fit)[-1L])
  if (length(slopes) <= 1L)
    return(numeric(0))

  anchor <- slopes[1L]
  tail <- slopes[-1L]
  finite.slopes <- slopes[is.finite(slopes)]
  scale <- if (length(finite.slopes)) max(1, max(abs(finite.slopes))) else 1
  if (!is.finite(anchor) || abs(anchor) <= sqrt(.Machine$double.eps) * scale) {
    tail[!is.finite(tail)] <- 0
    return(tail)
  }

  out <- tail / anchor
  out[!is.finite(out)] <- 0
  out
}

.npindex_index_from_beta_tail <- function(xmat, beta.tail) {
  xmat <- toMatrix(xmat)
  beta.tail <- as.double(beta.tail)
  beta <- c(1, beta.tail)
  if (length(beta) != ncol(xmat))
    stop("npindexbw: beta/index geometry length mismatch", call. = FALSE)
  index <- as.double(xmat %*% beta)
  index[!is.finite(index)] <- 0
  index
}

.npindex_beta_coordinate_setup <- function(xmat) {
  xmat <- toMatrix(xmat)
  p <- ncol(xmat)
  if (p <= 1L) {
    factor <- numeric(0)
  } else {
    scales <- apply(xmat, 2L, stats::sd)
    scales <- as.double(scales)
    scales[!is.finite(scales) | scales <= 0] <- 1
    anchor <- scales[1L]
    if (!is.finite(anchor) || anchor <= 0)
      anchor <- 1
    factor <- scales[-1L] / anchor
    factor[!is.finite(factor) | factor <= 0] <- 1
  }
  list(
    factor = factor,
    to_search = function(beta.tail) {
      beta.tail <- as.double(beta.tail)
      if (!length(beta.tail))
        return(numeric(0))
      beta.tail * factor
    },
    to_public = function(beta.search) {
      beta.search <- as.double(beta.search)
      if (!length(beta.search))
        return(numeric(0))
      beta.search / factor
    }
  )
}

.npindex_finalize_bandwidth <- function(h,
                                        bwtype,
                                        nobs,
                                        lower = NULL,
                                        where = "npindexbw") {
  candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bwtype, nobs = nobs)
  if (!candidate$ok) {
    if (identical(bwtype, "fixed")) {
      stop(sprintf("%s: bandwidth must be positive and finite", where), call. = FALSE)
    }
    upper <- max(2L, as.integer(nobs) - 1L)
    if (!identical(bwtype, "fixed") && is.finite(h) && h > upper &&
        !npExtendedNnEnabled()) {
      stop(
        sprintf(
          "%s: nearest-neighbor bandwidth exceeds n-1; set options(np.extendednn = TRUE) to allow extended generalized_nn/adaptive_nn bandwidths",
          where
        ),
        call. = FALSE
      )
    }
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth candidate must map to an integer in [2, %d]",
        where,
        max(2L, as.integer(nobs) - 1L)
      ),
      call. = FALSE
    )
  }

  if (identical(bwtype, "fixed") && !is.null(lower) && candidate$value < lower) {
    stop(sprintf("%s: bandwidth is below the continuous scale-factor lower bound", where),
         call. = FALSE)
  }

  candidate$value
}

.npindexbw_nomad_fixed_h_scale <- function(fit,
                                           h.start.raw,
                                           nobs,
                                           start.controls = .npindexbw_h_start_controls()) {
  scale <- .npindex_start_bandwidth_scale(fit = fit, nobs = nobs)
  if (!is.finite(scale) || scale <= 0) {
    scale <- max(
      if (is.finite(h.start.raw)) abs(as.double(h.start.raw)) else 0,
      1e-3
    )
  }
  as.double(scale)
}

.npindexbw_nomad_fixed_start_setup <- function(xmat,
                                               ydat,
                                               baseline.bws,
                                               degree.search,
                                               nmulti,
                                               random.seed,
                                               h.start.controls = .npindexbw_h_start_controls()) {
  p <- ncol(xmat)
  nobs <- nrow(xmat)
  beta.free <- if (p > 1L) seq_len(p - 1L) else integer(0)
  h.col <- length(beta.free) + 1L
  degree.col <- h.col + 1L

  ols.fit <- lm(ydat ~ xmat, x = TRUE)

  ols.beta <- if (length(beta.free)) {
    .npindex_ols_beta_tail(ols.fit)
  } else {
    numeric(0)
  }
  if (length(ols.beta))
    ols.beta[!is.finite(ols.beta)] <- 0

  beta.start.raw <- if (length(beta.free)) {
    beta.user <- as.double(baseline.bws$beta[beta.free + 1L])
    if (setequal(beta.user, c(0))) ols.beta else beta.user
  } else {
    numeric(0)
  }
  if (length(beta.start.raw))
    beta.start.raw[!is.finite(beta.start.raw)] <- 0
  fit.proxy <- .npindex_index_from_beta_tail(xmat, beta.start.raw)

  if (isTRUE(all.equal(as.double(baseline.bws$bw[1L]), 0))) {
    h.start.raw <- .npindex_default_start_bandwidth(
      fit = fit.proxy,
      bwtype = "fixed",
      nobs = nobs,
      start.controls = h.start.controls
    )
  } else {
    h.lower.raw <- h.start.controls$scale.factor.search.lower *
      .npindex_start_bandwidth_scale(fit = fit.proxy, nobs = nobs)
    h.start.raw <- tryCatch(
      .npindex_finalize_bandwidth(
        h = baseline.bws$bw[1L],
        bwtype = "fixed",
        nobs = nobs,
        lower = h.lower.raw,
        where = "npindexbw"
      ),
      error = function(e) .npindex_default_start_bandwidth(
        fit = fit.proxy,
        bwtype = "fixed",
        nobs = nobs,
        start.controls = h.start.controls
      )
    )
  }
  if (!is.finite(h.start.raw) || h.start.raw <= 0)
    h.start.raw <- 1e-3

  h.scale <- .npindexbw_nomad_fixed_h_scale(
    fit = fit.proxy,
    h.start.raw = h.start.raw,
    nobs = nobs,
    start.controls = h.start.controls
  )
  degree.starts <- .np_lp_nomad_build_degree_starts(
    initial = degree.search$start.degree,
    lower = degree.search$lower,
    upper = degree.search$upper,
    basis = degree.search$basis,
    nobs = degree.search$nobs,
    nmulti = nmulti,
    random.seed = random.seed,
    user_supplied = isTRUE(degree.search$start.user)
  )

  start_matrix.raw <- matrix(0, nrow = nmulti, ncol = degree.col)
  if (length(beta.free))
    start_matrix.raw[1L, beta.free] <- beta.start.raw
  start_matrix.raw[1L, h.col] <- h.start.raw
  start_matrix.raw[, degree.col] <- as.integer(degree.starts[, 1L])

  if (nmulti > 1L) {
    seed.state <- .np_seed_enter(random.seed)
    on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

    for (j in 2:nmulti) {
      if (length(beta.free)) {
        beta.rand <- runif(length(ols.beta), min = 0.5, max = 1.5) * ols.beta
        beta.rand[!is.finite(beta.rand)] <- 0
        start_matrix.raw[j, beta.free] <- beta.rand
      }

      h.rand <- .npindex_random_start_bandwidth(
        fit = fit.proxy,
        bwtype = "fixed",
        nobs = nobs,
        start.controls = h.start.controls
      )
      if (!is.finite(h.rand) || h.rand <= 0)
        h.rand <- h.start.raw
      start_matrix.raw[j, h.col] <- h.rand
    }
  }

  start_matrix.point <- start_matrix.raw
  start_matrix.point[, h.col] <- start_matrix.point[, h.col] / h.scale

  list(
    beta.free = beta.free,
    ols.beta = ols.beta,
    h.col = h.col,
    degree.col = degree.col,
    h.scale = h.scale,
    start_matrix.raw = start_matrix.raw,
    start_matrix.point = start_matrix.point
  )
}

.npindexbw_fast_eligible <- function(h, bws, eval.index) {
  if (!npLogicalOption("np.largeh", TRUE))
    return(FALSE)

  if (!identical(bws$type, "fixed"))
    return(FALSE)

  ckerbound <- if (is.null(bws$ckerbound) || !length(bws$ckerbound)) {
    "none"
  } else {
    as.character(bws$ckerbound)[1L]
  }
  if (!identical(ckerbound, "none"))
    return(FALSE)

  ckertype <- as.character(bws$ckertype)[1L]
  ckerorder <- as.integer(bws$ckerorder)[1L]
  if (identical(ckertype, "truncated gaussian"))
    return(FALSE)
  if (!identical(ckertype, "uniform") &&
      (!is.finite(ckerorder) || ckerorder != 2L))
    return(FALSE)

  fast_largeh_tol <- npLargehRelTol()

  cont_utol <- switch(
    ckertype,
    gaussian = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
    epanechnikov = sqrt(fast_largeh_tol),
    uniform = 1.0 - 32.0 * .Machine$double.eps,
    0.0
  )

  h <- as.double(h)
  if (length(h) != 1L)
    return(FALSE)
  if (!is.finite(cont_utol) || cont_utol <= 0 || !is.finite(h) || h <= 0)
    return(FALSE)

  vals <- as.double(eval.index)
  if (!length(vals) || any(!is.finite(vals)))
    return(FALSE)

  h >= (diff(range(vals)) / cont_utol)
}

.npindexbw_is_degree0_policy <- function(policy) {
  degree <- as.integer(policy$canonical.degree)
  length(degree) > 0L && all(degree == 0L)
}

.npindexbw_check_index_bound_contract <- function(bws,
                                                  policy,
                                                  where = "npindexbw") {
  if (!identical(as.character(policy$executor[1L]), "npreg_loo"))
    stop("internal error: npindex objective policy selected an unsupported executor", call. = FALSE)

  ckerbound <- if (is.null(bws$ckerbound) || !length(bws$ckerbound)) {
    "none"
  } else {
    as.character(bws$ckerbound[1L])
  }

  if (identical(ckerbound, "fixed")) {
    stop(
      sprintf(
        "%s does not support ckerbound='fixed' for single-index objective evaluation; use ckerbound='range' to compute bounds on the scalar index or ckerbound='none'",
        where
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.npindexbw_build_lp_regression_leaf <- function(index,
                                                ydat,
                                                h,
                                                bws,
                                                spec) {
  index.df <- data.frame(index = as.double(index))
  engine.regtype <- if (is.null(spec$regtype.engine) ||
                        !length(spec$regtype.engine)) {
    "lc"
  } else {
    as.character(spec$regtype.engine[1L])
  }
  engine.basis <- if (identical(engine.regtype, "lp")) {
    as.character(spec$basis.engine)
  } else {
    "glp"
  }
  engine.degree <- if (identical(engine.regtype, "lp")) {
    as.integer(spec$degree.engine)
  } else {
    NULL
  }
  engine.bernstein <- identical(engine.regtype, "lp") &&
    isTRUE(spec$bernstein.basis.engine)
  reg.args <- list(
    regtype = engine.regtype,
    basis = engine.basis,
    degree = engine.degree,
    bernstein.basis = engine.bernstein,
    bwmethod = "cv.ls",
    bwtype = bws$type,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder,
    ckerbound = bws$ckerbound,
    ckerlb = bws$ckerlb,
    ckerub = bws$ckerub,
    ukertype = if (is.null(bws$ukertype)) "aitchisonaitken" else bws$ukertype,
    okertype = if (is.null(bws$okertype)) "liracine" else bws$okertype
  )

  list(
    xdat = index.df,
    bws = .npregbw_build_rbandwidth(
      xdat = index.df,
      ydat = ydat,
      bws = c(h),
      bandwidth.compute = FALSE,
      reg.args = reg.args,
      yname = if (is.null(bws$ynames)) "y" else as.character(bws$ynames[1L])
    )
  )
}

.npindexbw_with_inner_bandwidth_progress_suppressed <- function(expr) {
  old.state <- .np_progress_runtime$bandwidth_state
  .np_progress_runtime$bandwidth_state <- NULL
  on.exit(.np_progress_runtime$bandwidth_state <- old.state, add = TRUE)
  force(expr)
}

.npindexbw_eval_ichimura_lp_via_npreg <- function(index,
                                                  ydat,
                                                  h,
                                                  bws,
                                                  spec,
                                                  invalid.penalty,
                                                  localize = TRUE) {
  leaf <- .npindexbw_build_lp_regression_leaf(
    index = index,
    ydat = ydat,
    h = h,
    bws = bws,
    spec = spec
  )

  out <- tryCatch(
    .npindexbw_with_inner_bandwidth_progress_suppressed(
      .npregbw_eval_only(
        xdat = leaf$xdat,
        ydat = ydat,
        bws = leaf$bws,
        invalid.penalty = "baseline",
        penalty.multiplier = 10,
        localize = localize
      )
    ),
    error = function(e) NULL
  )

  if (is.null(out) || !is.finite(out$objective[1L]))
    return(list(objective = as.numeric(invalid.penalty), num.feval.fast = 0L))

  list(
    objective = as.numeric(out$objective[1L]),
    num.feval.fast = as.numeric(out$num.feval.fast[1L])
  )
}

.npindexbw_eval_kleinspady_lp_via_npreg <- function(index,
                                                     ydat,
                                                     h,
                                                     bws,
                                                     spec,
                                                     invalid.penalty,
                                                     localize = TRUE) {
  leaf <- .npindexbw_build_lp_regression_leaf(
    index = index,
    ydat = ydat,
    h = h,
    bws = bws,
    spec = spec
  )

  out <- tryCatch(
    .npindexbw_with_inner_bandwidth_progress_suppressed(
      .npregbw_eval_only(
        xdat = leaf$xdat,
        ydat = ydat,
        bws = leaf$bws,
        invalid.penalty = "dbmax",
        penalty.multiplier = 10,
        localize = localize,
        objective = "ks"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(out) || !is.finite(out$objective[1L]))
    return(list(objective = as.numeric(invalid.penalty), num.feval.fast = 0L))

  list(
    objective = as.numeric(out$objective[1L]),
    num.feval.fast = as.numeric(out$num.feval.fast[1L])
  )
}

.npindexbw_ichimura_lp_service_context <- function(bws,
                                                   spec,
                                                   bandwidth.compute = TRUE,
                                                   comm = 1L,
                                                   service_id = "npindex_ichimura_lp",
                                                   child_service_ids = character(0L)) {
  rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) 0L)
  size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) 1L)
  if (!is.finite(rank)) rank <- 0L
  if (!is.finite(size) || size < 1L) size <- 1L

  in.context <- isTRUE(.npRmpi_autodispatch_in_context()) ||
    isTRUE(.npRmpi_manual_bcast_in_context()) ||
    isTRUE(.npRmpi_autodispatch_called_from_bcast())
  pool.ok <- if (isTRUE(in.context)) {
    size > 1L
  } else {
    isTRUE(.npRmpi_has_active_slave_pool(comm = comm))
  }

  active <- isTRUE(bandwidth.compute) &&
    identical(as.character(bws$method[1L]), "ichimura") &&
    isTRUE(pool.ok) &&
    !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
    isTRUE(in.context)

  list(
    active = isTRUE(active) && size > 1L,
    root = identical(rank, 0L),
    rank = rank,
    size = size,
    comm = as.integer(comm),
    service_id = as.character(service_id)[1L],
    child_service_ids = unique(as.character(child_service_ids))
  )
}

.npindexbw_service_eval_preflight <- function(ok,
                                              comm = 1L) {
  local.invalid <- if (isTRUE(ok)) 0.0 else 1.0
  total.invalid <- mpi.allreduce(as.double(local.invalid),
                                 type = 2,
                                 op = "sum",
                                 comm = comm)
  as.numeric(total.invalid[1L]) > 0.0
}

.npindexbw_service_param_ok <- function(param, expected.length) {
  param <- as.numeric(param)
  expected.length <- as.integer(expected.length)[1L]
  is.finite(expected.length) &&
    length(param) == expected.length &&
    all(is.finite(param))
}

.npindexbw_service_task_error <- function(role,
                                          message,
                                          task,
                                          ctx) {
  task.kind <- if (is.list(task) && !is.null(task$kind))
    as.character(task$kind)[1L]
  else
    NA_character_
  task.service <- if (is.list(task) && !is.null(task$service_id))
    as.character(task$service_id)[1L]
  else
    NA_character_
  err <- list(
    kind = "service_error",
    service_id = ctx$service_id,
    rank = ctx$rank,
    task_kind = task.kind,
    task_service_id = task.service,
    message = as.character(message)[1L]
  )
  .npRmpi_transport_trace(
    role = role,
    event = "task.error",
    fields = err
  )
  invisible(err)
}

.npindexbw_eval_objective_service_traced <- function(param,
                                                     xmat,
                                                     ydat,
                                                     bws,
                                                     spec,
                                                     ctx,
                                                     eval_id = NA_integer_) {
  assignments <- .splitIndices(nrow(xmat), ctx$size)
  local.idx <- if (length(assignments) >= (ctx$rank + 1L)) {
    assignments[[ctx$rank + 1L]]
  } else {
    integer(0L)
  }
  started <- proc.time()[3L]
  .npRmpi_transport_trace(
    role = "npindex.npreg_loo.service",
    event = "eval.enter",
    fields = list(
      service_id = ctx$service_id,
      executor = if (is.null(ctx$executor)) "npreg_loo" else ctx$executor,
      eval_id = eval_id,
      rank = ctx$rank,
      size = ctx$size,
      nominal_partition_rows = length(local.idx),
      objective_rows = nrow(xmat),
      localize = FALSE,
      n = nrow(xmat)
    )
  )
  out <- tryCatch(
    .npindexbw_eval_objective(
      param = param,
      xmat = xmat,
      ydat = ydat,
      bws = bws,
      spec = spec,
      localize = FALSE
    ),
    error = function(e) {
      list(
        objective = NA_real_,
        num.feval.fast = 0L,
        error = conditionMessage(e)
      )
    }
  )
  .npRmpi_transport_trace(
    role = "npindex.npreg_loo.service",
    event = "eval.exit",
    fields = list(
      service_id = ctx$service_id,
      executor = if (is.null(ctx$executor)) "npreg_loo" else ctx$executor,
      eval_id = eval_id,
      rank = ctx$rank,
      size = ctx$size,
      nominal_partition_rows = length(local.idx),
      objective_rows = nrow(xmat),
      localize = FALSE,
      elapsed = proc.time()[3L] - started,
      objective = if (is.null(out$objective)) NA_real_ else as.numeric(out$objective[1L]),
      ok = is.null(out$error),
      error = if (is.null(out$error)) "" else out$error
    )
  )
  if (!is.null(out$error))
    stop(out$error, call. = FALSE)
  out
}

.npindexbw_ichimura_lp_service_task_error <- function(message,
                                                      task,
                                                      ctx) {
  .npindexbw_service_task_error(
    role = "npindex.ichimura.lp.service",
    message = message,
    task = task,
    ctx = ctx
  )
}

.npindexbw_ichimura_lp_service_worker_loop <- function(xmat,
                                                       ydat,
                                                       bws,
                                                       spec,
                                                       ctx) {
  repeat {
    task <- mpi.bcast.Robj(rank = 0L, comm = ctx$comm)
    if (!is.list(task) || is.null(task$kind)) {
      .npindexbw_ichimura_lp_service_task_error(
        "malformed npindex Ichimura LP service task",
        task = task,
        ctx = ctx
      )
      next
    }

    task.service <- if (is.null(task$service_id)) ctx$service_id else as.character(task$service_id)[1L]
    child.service <- !identical(task.service, ctx$service_id) &&
      task.service %in% ctx$child_service_ids
    if (!identical(task.service, ctx$service_id) && !isTRUE(child.service)) {
      .npindexbw_ichimura_lp_service_task_error(
        "unexpected npindex Ichimura LP service identifier",
        task = task,
        ctx = ctx
      )
      next
    }
    task.ctx <- ctx
    ## The NOMAD owner may run a declared fixed-degree child service for
    ## Powell refinement; only the parent service returns to the public caller.
    if (isTRUE(child.service))
      task.ctx$service_id <- task.service

    if (identical(task$kind, "eval")) {
      param <- if (is.null(task$param)) numeric(0L) else as.numeric(task$param)
      invalid.eval <- .npindexbw_service_eval_preflight(
        ok = .npindexbw_service_param_ok(param, ncol(xmat)),
        comm = task.ctx$comm
      )
      if (isTRUE(invalid.eval)) {
        .npindexbw_ichimura_lp_service_task_error(
          "malformed npindex Ichimura LP eval task",
          task = task,
          ctx = task.ctx
        )
        next
      }
      task.spec <- if (is.null(task$spec)) spec else task$spec
      .npindexbw_eval_objective_service_traced(
        param = param,
        xmat = xmat,
        ydat = ydat,
        bws = bws,
        spec = task.spec,
        ctx = task.ctx,
        eval_id = if (is.null(task$eval_id)) NA_integer_ else as.integer(task$eval_id[1L])
      )
      next
    }

    if (identical(task$kind, "result")) {
      if (isTRUE(child.service))
        next
      return(task$value)
    }
    if (identical(task$kind, "error")) {
      if (isTRUE(child.service))
        next
      stop(as.character(task$message)[1L], call. = FALSE)
    }
    if (identical(task$kind, "stop")) {
      if (isTRUE(child.service))
        next
      return(invisible(NULL))
    }

    .npindexbw_ichimura_lp_service_task_error(
      "unknown npindex Ichimura LP service task",
      task = task,
      ctx = ctx
    )
  }
}

.npindexbw_ichimura_lp_service_eval <- function(param,
                                                xmat,
                                                ydat,
                                                bws,
                                                spec,
                                                ctx,
                                                eval_id = NA_integer_) {
  mpi.bcast.Robj(
    list(
      kind = "eval",
      service_id = ctx$service_id,
      eval_id = eval_id,
      param = as.numeric(param),
      spec = spec
    ),
    rank = 0L,
    comm = ctx$comm
  )

  invalid.eval <- .npindexbw_service_eval_preflight(
    ok = .npindexbw_service_param_ok(param, ncol(xmat)),
    comm = ctx$comm
  )
  if (isTRUE(invalid.eval))
    return(list(
      objective = 10 * mean(ydat^2),
      num.feval.fast = 0L
    ))

  .npindexbw_eval_objective_service_traced(
    param = as.numeric(param),
    xmat = xmat,
    ydat = ydat,
    bws = bws,
    spec = spec,
    ctx = ctx,
    eval_id = eval_id
  )
}

.npindexbw_ichimura_lp_service_result <- function(value, ctx) {
  mpi.bcast.Robj(
    list(kind = "result", service_id = ctx$service_id, value = value),
    rank = 0L,
    comm = ctx$comm
  )
  value
}

.npindexbw_ichimura_lp_service_error <- function(message, ctx) {
  try(
    mpi.bcast.Robj(
      list(kind = "error", service_id = ctx$service_id, message = as.character(message)[1L]),
      rank = 0L,
      comm = ctx$comm
    ),
    silent = TRUE
  )
  invisible(NULL)
}

.npindexbw_kleinspady_lp_service_context <- function(bws,
                                                     spec,
                                                     policy,
                                                     bandwidth.compute = TRUE,
                                                     comm = 1L,
                                                     service_id = "npindex_kleinspady_lp") {
  rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) 0L)
  size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) 1L)
  if (!is.finite(rank)) rank <- 0L
  if (!is.finite(size) || size < 1L) size <- 1L

  in.context <- isTRUE(.npRmpi_autodispatch_in_context()) ||
    isTRUE(.npRmpi_manual_bcast_in_context()) ||
    isTRUE(.npRmpi_autodispatch_called_from_bcast())
  pool.ok <- if (isTRUE(in.context)) {
    size > 1L
  } else {
    isTRUE(.npRmpi_has_active_slave_pool(comm = comm))
  }

  executor <- as.character(policy$executor[1L])

  active <- isTRUE(bandwidth.compute) &&
    identical(as.character(bws$method[1L]), "kleinspady") &&
    identical(executor, "npreg_loo") &&
    isTRUE(pool.ok) &&
    !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
    isTRUE(in.context)

  list(
    active = isTRUE(active) && size > 1L,
    root = identical(rank, 0L),
    rank = rank,
    size = size,
    comm = as.integer(comm),
    service_id = as.character(service_id)[1L],
    executor = executor
  )
}

.npindexbw_kleinspady_lp_service_task_error <- function(message,
                                                        task,
                                                        ctx) {
  .npindexbw_service_task_error(
    role = "npindex.kleinspady.lp.service",
    message = message,
    task = task,
    ctx = ctx
  )
}

.npindexbw_kleinspady_lp_service_worker_loop <- function(xmat,
                                                         ydat,
                                                         bws,
                                                         spec,
                                                         ctx) {
  repeat {
    task <- mpi.bcast.Robj(rank = 0L, comm = ctx$comm)
    if (!is.list(task) || is.null(task$kind)) {
      .npindexbw_kleinspady_lp_service_task_error(
        "malformed npindex Klein-Spady LP service task",
        task = task,
        ctx = ctx
      )
      next
    }

    task.service <- if (is.null(task$service_id)) ctx$service_id else as.character(task$service_id)[1L]
    if (!identical(task.service, ctx$service_id)) {
      .npindexbw_kleinspady_lp_service_task_error(
        "unexpected npindex Klein-Spady LP service identifier",
        task = task,
        ctx = ctx
      )
      next
    }

    if (identical(task$kind, "eval")) {
      beta <- if (is.null(task$beta)) numeric(0L) else as.double(task$beta)
      h <- if (is.null(task$h)) NA_real_ else as.double(task$h[1L])
      invalid.eval <- .npindexbw_service_eval_preflight(
        ok = .npindexbw_service_param_ok(c(beta, h), ncol(xmat)),
        comm = ctx$comm
      )
      if (isTRUE(invalid.eval)) {
        .npindexbw_kleinspady_lp_service_task_error(
          "malformed npindex Klein-Spady LP eval task",
          task = task,
          ctx = ctx
        )
        next
      }
      task.spec <- if (is.null(task$spec)) spec else task$spec
      .npindexbw_eval_objective_service_traced(
        param = c(beta, h),
        xmat = xmat,
        ydat = ydat,
        bws = bws,
        spec = task.spec,
        ctx = ctx,
        eval_id = if (is.null(task$eval_id)) NA_integer_ else as.integer(task$eval_id[1L])
      )
      next
    }

    if (identical(task$kind, "result"))
      return(task$value)
    if (identical(task$kind, "error"))
      stop(as.character(task$message)[1L], call. = FALSE)
    if (identical(task$kind, "stop"))
      return(invisible(NULL))

    .npindexbw_kleinspady_lp_service_task_error(
      "unknown npindex Klein-Spady LP service task",
      task = task,
      ctx = ctx
    )
  }
}

.npindexbw_kleinspady_lp_service_eval <- function(beta,
                                                  h,
                                                  xmat,
                                                  ydat,
                                                  bws,
                                                  spec,
                                                  ctx,
                                                  eval_id = NA_integer_) {
  mpi.bcast.Robj(
    list(
      kind = "eval",
      service_id = ctx$service_id,
      eval_id = eval_id,
      beta = as.numeric(beta),
      h = as.numeric(h)[1L],
      spec = spec
    ),
    rank = 0L,
    comm = ctx$comm
  )

  invalid.eval <- .npindexbw_service_eval_preflight(
    ok = .npindexbw_service_param_ok(c(beta, h), ncol(xmat)),
    comm = ctx$comm
  )
  if (isTRUE(invalid.eval))
    return(list(
      objective = sqrt(.Machine$double.xmax),
      invalid = TRUE
    ))

  objective <- .npindexbw_eval_objective_service_traced(
    param = c(as.numeric(beta), as.numeric(h)[1L]),
    xmat = xmat,
    ydat = ydat,
    bws = bws,
    spec = spec,
    ctx = ctx,
    eval_id = eval_id
  )
  list(
    objective = as.numeric(objective$objective[1L]),
    invalid = !is.finite(as.numeric(objective$objective[1L]))
  )
}

.npindexbw_kleinspady_lp_service_result <- function(value, ctx) {
  mpi.bcast.Robj(
    list(kind = "result", service_id = ctx$service_id, value = value),
    rank = 0L,
    comm = ctx$comm
  )
  value
}

.npindexbw_kleinspady_lp_service_error <- function(message, ctx) {
  try(
    mpi.bcast.Robj(
      list(kind = "error", service_id = ctx$service_id, message = as.character(message)[1L]),
      rank = 0L,
      comm = ctx$comm
    ),
    silent = TRUE
  )
  invisible(NULL)
}

.npindexbw_build_sibandwidth <- function(xdat,
                                         ydat,
                                         bws,
                                         template,
                                         bandwidth.compute,
                                         reg.args) {
  p <- ncol(xdat)
  out <- do.call(
    sibandwidth,
    c(
      list(
        beta = bws[seq_len(p)],
        h = bws[p + 1L],
        method = template$method,
        ckertype = template$ckertype,
        ckerorder = template$ckerorder,
        ckerbound = template$ckerbound,
        ckerlb = template$ckerlb,
        ckerub = template$ckerub,
        bwtype = template$type,
        nobs = nrow(xdat),
        xdati = untangle(xdat),
        ydati = untangle(data.frame(ydat)),
        xnames = names(xdat),
        ynames = template$ynames,
        bandwidth = bws[p + 1L],
        bandwidth.compute = bandwidth.compute
      ),
      reg.args
    )
  )
  if (!is.null(reg.args$scale.factor.search.lower))
    out$scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      reg.args$scale.factor.search.lower
    )
  out
}

.npindexbw_run_fixed_degree <- function(xdat, ydat, bws, template, reg.args, opt.args) {
  tbw <- .npindexbw_build_sibandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    template = template,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args
  )

  do.call(npindexbw.sibandwidth, c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args))
}

.npindexbw_eval_objective <- function(param,
                                      xmat,
                                      ydat,
                                      bws,
                                      spec,
                                      localize = TRUE) {
  p <- ncol(xmat)
  beta.idx <- if (p > 1L) seq_len(p - 1L) else integer(0)
  beta <- if (length(beta.idx)) as.double(param[beta.idx]) else numeric(0)
  h <- as.double(param[p])
  nobs <- nrow(xmat)
  policy <- .npindex_objective_policy(
    bws = bws,
    spec = spec,
    bandwidth.compute = TRUE,
    where = "npindexbw objective"
  )
  .npindexbw_check_index_bound_contract(
    bws = bws,
    policy = policy,
    where = "npindexbw objective"
  )
  spec <- policy$objective.spec

  if (identical(bws$method, "ichimura")) {
    invalid.penalty <- 10 * mean(ydat^2)
  } else {
    invalid.penalty <- sqrt(.Machine$double.xmax)
  }

  h.candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bws$type, nobs = nobs)
  if (!h.candidate$ok)
    return(list(objective = invalid.penalty, num.feval.fast = 0L))
  h <- h.candidate$value

  index <- xmat %*% c(1, beta)

  if (identical(bws$method, "ichimura")) {
    return(.npindexbw_eval_ichimura_lp_via_npreg(
      index = index,
      ydat = ydat,
      h = h,
      bws = bws,
      spec = spec,
      invalid.penalty = invalid.penalty,
      localize = localize
    ))
  }

  if (identical(bws$method, "kleinspady")) {
    return(.npindexbw_eval_kleinspady_lp_via_npreg(
      index = index,
      ydat = ydat,
      h = h,
      bws = bws,
      spec = spec,
      invalid.penalty = invalid.penalty,
      localize = localize
    ))
  }

  stop("unsupported npindex method", call. = FALSE)
}

.npindexbw_nomad_search <- function(xdat,
                                    ydat,
                                    bws,
                                    template,
                                    reg.args,
                                    opt.args,
                                    degree.search,
                                    nomad.inner.nmulti = 0L,
                                    nomad.opts = list(),
                                    source = "explicit",
                                    reason = NULL,
                                    handoff_before_build = FALSE,
                                    progress_label = NULL) {
  if (is.null(opt.args$nomad.opts) && length(nomad.opts))
    opt.args$nomad.opts <- nomad.opts

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template.reg.args$regtype.engine <- "lp"
  template.reg.args$degree.engine <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

  baseline.bws <- .npindexbw_build_sibandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    template = template,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args
  )

  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE
  x.clean <- toMatrix(xdat[keep.rows, , drop = FALSE])
  y.clean <- ydat[keep.rows]
  if (is.factor(y.clean))
    y.clean <- dlev(y.clean)[as.integer(y.clean)]
  else
    y.clean <- as.double(y.clean)

  p <- ncol(x.clean)
  beta.coord <- .npindex_beta_coordinate_setup(x.clean)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(ncol(xdat)) else npValidateNmulti(opt.args$nmulti[1L])
  scale.factor.search.lower <- npResolveScaleFactorLowerBound(opt.args$scale.factor.search.lower)
  h.start.controls <- .npindexbw_h_start_controls(
    scale.factor.init.lower = if (is.null(opt.args$scale.factor.init.lower)) 0.1 else opt.args$scale.factor.init.lower,
    scale.factor.init.upper = if (is.null(opt.args$scale.factor.init.upper)) 2.0 else opt.args$scale.factor.init.upper,
    scale.factor.init = if (is.null(opt.args$scale.factor.init)) 0.5 else opt.args$scale.factor.init,
    scale.factor.search.lower = scale.factor.search.lower,
    where = "npindexbw"
  )
  fixed.nomad <- identical(baseline.bws$type, "fixed")
  fixed.setup <- if (fixed.nomad) {
    .npindexbw_nomad_fixed_start_setup(
      xmat = x.clean,
      ydat = y.clean,
      baseline.bws = baseline.bws,
      degree.search = degree.search,
      nmulti = nomad.nmulti,
      random.seed = if (!is.null(opt.args$random.seed)) opt.args$random.seed else 42L,
      h.start.controls = h.start.controls
    )
  } else {
    NULL
  }
  beta.free <- if (fixed.nomad) fixed.setup$beta.free else if (p > 1L) seq_len(p - 1L) else integer(0)
  beta.start <- if (fixed.nomad) {
    if (length(beta.free)) as.double(fixed.setup$start_matrix.raw[1L, beta.free]) else numeric(0)
  } else {
    if (length(beta.free)) as.double(baseline.bws$beta[beta.free + 1L]) else numeric(0)
  }
  beta.search.start <- beta.coord$to_search(beta.start)
  h.start.raw <- if (fixed.nomad) {
    as.double(fixed.setup$start_matrix.raw[1L, fixed.setup$h.col])
  } else {
    as.double(baseline.bws$bw[1L])
  }
  beta.lower <- if (length(beta.start)) {
    if (fixed.nomad) {
      beta.coord$to_search(pmin(0.5 * fixed.setup$ols.beta, 1.5 * fixed.setup$ols.beta))
    } else {
      beta.coord$to_search(-pmax(10, 10 * abs(beta.start)))
    }
  } else {
    numeric(0)
  }
  beta.upper <- if (length(beta.start)) {
    if (fixed.nomad) {
      beta.coord$to_search(pmax(0.5 * fixed.setup$ols.beta, 1.5 * fixed.setup$ols.beta))
    } else {
      beta.coord$to_search(pmax(10, 10 * abs(beta.start)))
    }
  } else {
    numeric(0)
  }
  h.lower.raw <- if (fixed.nomad) 1e-3 else 2
  h.upper.raw <- if (fixed.nomad) {
    max(1e6, abs(h.start.raw) * 1e3, 1)
  } else {
    max(2L, as.integer(nrow(x.clean)) - 1L)
  }
  h.lower <- if (fixed.nomad) {
    h.start.controls$scale.factor.search.lower
  } else {
    h.lower.raw
  }
  h.upper <- if (fixed.nomad) {
    h.upper.raw / fixed.setup$h.scale
  } else {
    h.upper.raw
  }
  h.integer <- !identical(baseline.bws$type, "fixed")
  coordinate.roles <- c(
    rep.int("continuous_real", length(beta.lower)),
    if (isTRUE(h.integer)) "continuous_nn_index" else "continuous_fixed_scale",
    rep.int("degree", length(degree.search$lower))
  )

  x0 <- if (fixed.nomad) {
    start <- as.numeric(fixed.setup$start_matrix.point[1L, ])
    if (length(beta.free))
      start[seq_along(beta.free)] <- beta.coord$to_search(start[seq_along(beta.free)])
    start
  } else {
    c(beta.search.start, h.start.raw, as.integer(degree.search$start.degree))
  }
  lb <- c(beta.lower, h.lower, degree.search$lower)
  ub <- c(beta.upper, h.upper, degree.search$upper)
  start.lb <- lb
  start.ub <- ub
  if (fixed.nomad) {
    h.idx <- length(beta.start) + 1L
    start.lb[h.idx] <- max(start.lb[h.idx], h.start.controls$scale.factor.init.lower)
    start.ub[h.idx] <- min(start.ub[h.idx], h.start.controls$scale.factor.init.upper)
    if (start.ub[h.idx] < start.lb[h.idx]) {
      stop("npindexbw: effective NOMAD fixed-bandwidth random-start interval is empty after applying search bounds",
           call. = FALSE)
    }
  }
  bbin <- c(rep.int(0L, length(beta.start)), if (isTRUE(h.integer)) 1L else 0L, 1L)
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0
  service.spec <- template.reg.args
  service.ctx <- .npindexbw_ichimura_lp_service_context(
    bws = baseline.bws,
    spec = service.spec,
    bandwidth.compute = TRUE,
    comm = 1L,
    service_id = "npindex_ichimura_lp_nomad",
    child_service_ids = if (identical(degree.search$engine, "nomad+powell") ||
                            isTRUE(degree.search$verify)) {
      "npindex_ichimura_lp_fixed"
    } else {
      character(0L)
    }
  )
  if (isTRUE(service.ctx$active) && !isTRUE(service.ctx$root))
    return(.npindexbw_ichimura_lp_service_worker_loop(
      xmat = x.clean,
      ydat = y.clean,
      bws = baseline.bws,
      spec = service.spec,
      ctx = service.ctx
    ))
  service.eval.counter <- 0L
  service.done <- FALSE
  if (isTRUE(service.ctx$active) && isTRUE(service.ctx$root)) {
    on.exit({
      if (!isTRUE(service.done))
        .npindexbw_ichimura_lp_service_error("npindex Ichimura LP NOMAD service stopped before returning a result", service.ctx)
    }, add = TRUE)
  }

  .np_nomad_baseline_note(degree.search$start.degree)

  point_h_to_raw <- function(h.point) {
    h.raw <- if (fixed.nomad) as.double(h.point) * fixed.setup$h.scale else as.double(h.point)
    as.double(h.raw)
  }

  point_to_public <- function(point) {
    point <- as.numeric(point)
    if (length(beta.free))
      point[seq_along(beta.free)] <- beta.coord$to_public(point[seq_along(beta.free)])
    if (length(point) >= (length(beta.free) + 1L))
      point[length(beta.free) + 1L] <- point_h_to_raw(point[length(beta.free) + 1L])
    point
  }

  bandwidth_point_to_public <- function(point) {
    point <- as.numeric(point)
    if (length(beta.free))
      point[seq_along(beta.free)] <- beta.coord$to_public(point[seq_along(beta.free)])
    if (length(point))
      point[length(point)] <- point_h_to_raw(point[length(point)])
    point
  }

  point_to_param <- function(point) {
    beta.tail <- if (length(beta.free)) beta.coord$to_public(point[seq_along(beta.free)]) else numeric(0)
    h <- point_h_to_raw(point[length(beta.free) + 1L])
    h <- .npindex_finalize_bandwidth(
      h = h,
      bwtype = baseline.bws$type,
      nobs = nrow(x.clean),
      lower = if (fixed.nomad) h.start.controls$scale.factor.search.lower * fixed.setup$h.scale else NULL,
      where = "npindexbw"
    )
    c(beta.tail, h)
  }

  point_to_bws <- function(point) {
    param <- point_to_param(point)
    beta.tail <- if (length(beta.free)) param[seq_along(beta.free)] else numeric(0)
    h <- param[length(beta.free) + 1L]
    c(1.0, beta.tail, h)
  }

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- .np_degree_clip_to_grid(
      as.integer(round(point[length(point)])),
      degree.search$candidates
    )
    eval.spec <- reg.args
    eval.spec$regtype.engine <- "lp"
    eval.spec$degree.engine <- degree
    eval.spec$bernstein.basis.engine <- degree.search$bernstein.basis
    eval.spec$basis.engine <- reg.args$basis.engine
    service.eval.counter <<- service.eval.counter + 1L
    objective <- if (isTRUE(service.ctx$active) && isTRUE(service.ctx$root)) {
      .npindexbw_ichimura_lp_service_eval(
        param = point_to_param(point),
        xmat = x.clean,
        ydat = y.clean,
        bws = baseline.bws,
        spec = eval.spec,
        ctx = service.ctx,
        eval_id = service.eval.counter
      )
    } else {
      .npindexbw_eval_objective(
        param = point_to_param(point),
        xmat = x.clean,
        ydat = y.clean,
        bws = baseline.bws,
        spec = eval.spec
      )
    }
    nomad.num.feval.total <<- nomad.num.feval.total + 1L
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(objective$num.feval.fast[1L])
    list(
      objective = as.numeric(objective$objective[1L]),
      degree = as.integer(degree),
      num.feval = 1,
      num.feval.fast = as.numeric(objective$num.feval.fast[1L])
    )
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    degree <- as.integer(best_record$degree)
    bw.vec <- point_to_bws(point)
    powell.elapsed <- NA_real_

    build_direct_payload <- function() {
      final.reg.args <- reg.args
      final.reg.args$regtype <- "lp"
      final.reg.args$degree <- degree
      final.reg.args$bernstein.basis <- degree.search$bernstein.basis
      final.reg.args$regtype.engine <- "lp"
      final.reg.args$degree.engine <- degree
      final.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

      tbw <- .npindexbw_build_sibandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bw.vec,
        template = baseline.bws,
        bandwidth.compute = FALSE,
        reg.args = final.reg.args
      )
      tbw$fval <- as.numeric(best_record$objective)
      tbw$num.feval <- as.numeric(nomad.num.feval.total)
      tbw$num.feval.fast <- as.numeric(nomad.num.feval.fast.total)
      tbw$total.time <- NA_real_
      payload <- npindexbw.sibandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = tbw,
        bandwidth.compute = FALSE
      )
      if (!is.null(payload$method) && length(payload$method))
        payload$pmethod <- bwmToPrint(as.character(payload$method[1L]))
      payload
    }

    direct.payload <- build_direct_payload()
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.reg.args$regtype.engine <- "lp"
      hot.reg.args$degree.engine <- degree
      hot.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
      hot.opt.args <- .np_nomad_powell_hotstart_opt_args(
        opt.args,
        strategy = "single_iteration",
        remin = isTRUE(opt.args$powell.remin)
      )
      hot.opt.args$bwsolver <- NULL
      powell.start <- proc.time()[3L]
      hot.payload <- .np_nomad_with_powell_progress(
        degree = degree,
        best_record = best_record,
        expr = .npindexbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          bws = bw.vec,
          template = template,
          reg.args = hot.reg.args,
          opt.args = hot.opt.args
        )
      )
      powell.elapsed <- proc.time()[3L] - powell.start
      direct.payload$num.feval <- as.numeric(nomad.num.feval.total) + as.numeric(hot.payload$num.feval[1L])
      direct.payload$num.feval.fast <- as.numeric(nomad.num.feval.fast.total) + as.numeric(hot.payload$num.feval.fast[1L])
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

  search.result <- .np_nomad_search(
    engine = degree.search$engine,
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
    random.seed = if (!is.null(opt.args$random.seed)) opt.args$random.seed else 42L,
    remin = isTRUE(opt.args$nomad.remin),
    source = source,
    reason = reason,
    progress_label = progress_label,
    handoff_before_build = isTRUE(handoff_before_build),
    nomad.opts = if (is.null(opt.args$nomad.opts)) list() else opt.args$nomad.opts,
    native.r.bridge = TRUE,
    start.lower = start.lb,
    start.upper = start.ub,
    coordinate.roles = coordinate.roles,
    degree_spec = list(
      initial = degree.search$start.degree,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      user_supplied = degree.search$start.user
    )
  )

  if (isTRUE(degree.search$verify)) {
    verify.started <- proc.time()[3L]
    verify.records <- list()
    verify.best <- search.result$best
    verify.best.payload <- search.result$best_payload
    verify.eval.id <- 0L

    visit_degree <- function(degree.vec) {
      verify.eval.id <<- verify.eval.id + 1L
      verify.reg.args <- reg.args
      verify.reg.args$regtype <- "lp"
      verify.reg.args$degree <- as.integer(degree.vec)
      verify.reg.args$bernstein.basis <- degree.search$bernstein.basis
      verify.reg.args$regtype.engine <- "lp"
      verify.reg.args$degree.engine <- as.integer(degree.vec)
      verify.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

      started <- proc.time()[3L]
      payload <- tryCatch(
        .npindexbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          template = template,
          reg.args = verify.reg.args,
          opt.args = opt.args
        ),
        error = function(e) e
      )
      status <- if (inherits(payload, "error")) "error" else "ok"
      objective <- if (identical(status, "ok")) as.numeric(payload$fval[1L]) else NA_real_
      rec <- list(
        eval_id = verify.eval.id,
        degree = as.integer(degree.vec),
        objective = objective,
        status = status,
        cached = FALSE,
        message = if (inherits(payload, "error")) conditionMessage(payload) else NULL,
        elapsed = proc.time()[3L] - started,
        num.feval = if (identical(status, "ok") && !is.null(payload$num.feval)) {
          as.numeric(payload$num.feval[1L])
        } else {
          NA_real_
        },
        phase = "verify"
      )
      verify.records[[length(verify.records) + 1L]] <<- rec
      if (identical(status, "ok") &&
          (is.null(verify.best) ||
           .np_degree_better(objective, verify.best$objective, direction = "min"))) {
        verify.best <<- rec
        verify.best.payload <<- payload
      }
      TRUE
    }

    .np_degree_check_grid_budget(
      candidates = degree.search$candidates,
      method = "coordinate",
      verify = TRUE
    )
    .np_degree_iterate_grid(degree.search$candidates, visit_degree)
    verify.elapsed <- proc.time()[3L] - verify.started

    search.result$verify <- TRUE
    search.result$certified <- TRUE
    search.result$verify.time <- verify.elapsed
    search.result$verify.results <- verify.records
    search.result$best <- verify.best
    search.result$best_payload <- verify.best.payload
    search.result$n.unique <- as.numeric(search.result$n.unique) + verify.eval.id
    search.result$n.visits <- as.numeric(search.result$n.visits) + verify.eval.id
    search.result$grid.size <- .np_degree_grid_size(degree.search$candidates)
    search.result$optim.time <- sum(c(search.result$optim.time, verify.elapsed), na.rm = TRUE)
  }

  if (fixed.nomad) {
    if (!is.null(search.result$restart.starts)) {
      search.result$restart.starts <- lapply(search.result$restart.starts, point_to_public)
    }
    if (!is.null(search.result$restart.bandwidth.starts)) {
      search.result$restart.bandwidth.starts <- lapply(search.result$restart.bandwidth.starts, bandwidth_point_to_public)
    }
    if (!is.null(search.result$best_point) && length(search.result$best_point)) {
      search.result$best_point <- point_to_public(search.result$best_point)
    }
    if (!is.null(search.result$restart.results) && length(search.result$restart.results)) {
      search.result$restart.results <- lapply(search.result$restart.results, function(rec) {
        if (!is.null(rec$start) && length(rec$start))
          rec$start <- point_to_public(rec$start)
        if (!is.null(rec$solution) && length(rec$solution))
          rec$solution <- point_to_public(rec$solution)
        rec
      })
    }
  }

  if (isTRUE(service.ctx$active) && isTRUE(service.ctx$root)) {
    service.done <- TRUE
    return(.npindexbw_ichimura_lp_service_result(search.result, service.ctx))
  }

  search.result
}

.npindexbw_degree_search_controls <- function(regtype,
                                              regtype.named,
                                              bandwidth.compute,
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
    ncon = 1L,
    source = nomad.source,
    auto.filled = nomad.auto.filled
  )
  search.engine <- .np_degree_search_engine_controls(resolved$search.engine)
  degree.select <- resolved$degree.select

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (!isTRUE(bandwidth.compute))
    stop("automatic degree search requires bandwidth.compute=TRUE")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = 1L,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  baseline.degree <- 0L
  default.start.degree <- if (identical(search.engine, "cell")) baseline.degree else 1L
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, default.start.degree))
  } else {
    start.raw <- npValidateGlpDegree(
      regtype = "lp",
      degree = degree.start,
      ncon = 1L,
      argname = "degree.start"
    )
    if (!(start.raw[1L] %in% bounds$candidates[[1L]]))
      stop("degree.start must lie within the searched degree candidates")
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
    baseline.degree = as.integer(baseline.degree),
    start.degree = as.integer(start.degree),
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

.npindexbw_attach_degree_search <- function(bws, search_result) {
  metadata <- list(
    mode = search_result$method,
    source = if (!is.null(search_result$source)) search_result$source else "explicit",
    reason = if (!is.null(search_result$reason)) search_result$reason else NULL,
    engine = if (!is.null(search_result$engine)) search_result$engine else search_result$method,
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
    singleton = isTRUE(search_result$singleton),
    fixed.degree = search_result$fixed.degree,
    best.restart = search_result$best.restart,
    restart.starts = search_result$restart.starts,
    restart.degree.starts = search_result$restart.degree.starts,
    restart.bandwidth.starts = search_result$restart.bandwidth.starts,
    restart.start.info = search_result$restart.start.info,
    restart.results = search_result$restart.results,
    nn.cache = search_result$nn.cache,
    trace = search_result$trace
  )

  if (isTRUE(search_result$native) &&
      isTRUE(getOption("np.developer.native.nomad.diagnostics", FALSE)) &&
      !is.null(search_result$native.diagnostics)) {
    attr(bws, "native.nomad.diagnostics") <- search_result$native.diagnostics
  }

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

npindexbw.default <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws,
           bandwidth.compute = TRUE,
           basis = c("glp", "additive", "tensor"),
           bernstein.basis = FALSE,
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
           nmulti,
           nomad.remin = FALSE,
           powell.remin = TRUE,
           only.optimize.beta,
           optim.abstol,
           optim.maxattempts,
           optim.maxit,
           optim.method,
           optim.reltol,
           random.seed,
           regtype = c("lc", "ll", "lp"),
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           scale.factor.search.lower = NULL,
           ...){
    .npRmpi_require_active_slave_pool(where = "npindexbw()")
    mc <- match.call(expand.dots = FALSE)
    search.mc.names <- names(mc)
    dots <- list(...)
    npRejectUnsupportedBwsolver(dots, "npindexbw")
    dot.names <- names(dots)
    nomad.mode <- npValidateNomadControl(nomad, "nomad")
    degree.select.value <- if (nomad.mode %in% c("true", "auto")) {
      "coordinate"
    } else if ("degree.select" %in% search.mc.names) {
      degree.select
    } else {
      "manual"
    }
    automatic.degree.search <- !identical(match.arg(degree.select.value, c("manual", "coordinate", "exhaustive")), "manual")
    method.value <- if ("method" %in% dot.names) {
      match.arg(as.character(dots$method[[1L]]), c("ichimura", "kleinspady"))
    } else {
      "ichimura"
    }
    collective.degree.search <- isTRUE(automatic.degree.search) &&
      identical(method.value, "kleinspady")
    nomad.requested <- nomad.mode %in% c("true", "auto")
    regtype.value <- if ("regtype" %in% search.mc.names) {
      match.arg(regtype, c("lc", "ll", "lp"))
    } else if (isTRUE(nomad.requested)) {
      "lp"
    } else {
      "lc"
    }
    search.engine.value <- if ("search.engine" %in% search.mc.names) {
      match.arg(search.engine, c("nomad+powell", "cell", "nomad"))
    } else {
      "nomad+powell"
    }
    ichimura.lp.nomad.degree.search <- isTRUE(automatic.degree.search) &&
      identical(method.value, "ichimura") &&
      identical(regtype.value, "lp") &&
      search.engine.value %in% c("nomad", "nomad+powell")
    if (.npRmpi_autodispatch_active() &&
        (!isTRUE(automatic.degree.search) ||
           isTRUE(collective.degree.search) ||
           isTRUE(ichimura.lp.nomad.degree.search) ||
           .npRmpi_safe_int(mpi.comm.size(1L)) > 2L))
      return(.npRmpi_autodispatch_call(mc, parent.frame()))

    xdat <- toFrame(xdat)

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")

      .np_warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    if (coarseclass(bws) != "numeric" || length(bws) != ncol(xdat)+1)
      stop(paste("manually specified 'bws' must be a numeric vector of length ncol(xdat)+1.",
                 "See documentation for details."))

    p <- ncol(xdat)
    mc.names <- names(mc)
    nomad.shortcut <- .np_prepare_nomad_shortcut(
      nomad = nomad,
      call_names = unique(c(mc.names, dot.names)),
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
        bwtype = if ("bwtype" %in% dot.names) dots$bwtype else NULL,
        degree = if ("degree" %in% mc.names) degree else NULL
      ),
      where = "npindexbw"
    )

    if (isTRUE(nomad.shortcut$enabled)) {
      if ("degree" %in% mc.names)
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% dot.names &&
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
    }

    random.seed.value <- if ("random.seed" %in% mc.names) {
      npValidateNonNegativeInteger(random.seed, "random.seed")
    } else {
      42L
    }
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else degree.select.value
    degree.search <- .npindexbw_degree_search_controls(
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      regtype.named = isTRUE(nomad.shortcut$enabled) || ("regtype" %in% mc.names),
      bandwidth.compute = bandwidth.compute,
      nobs = NROW(xdat),
      basis = basis,
      degree.select = degree.select.value,
      search.engine = if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell",
      degree.min = nomad.shortcut$values$degree.min,
      degree.max = nomad.shortcut$values$degree.max,
      degree.start = if ("degree.start" %in% mc.names) degree.start else NULL,
      degree.restarts = if ("degree.restarts" %in% mc.names) degree.restarts else 0L,
      degree.max.cycles = if ("degree.max.cycles" %in% mc.names) degree.max.cycles else 20L,
      degree.verify = if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE,
      bernstein.basis = if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else bernstein.basis,
      bernstein.named = isTRUE(nomad.shortcut$enabled) || ("bernstein.basis" %in% mc.names),
      nomad.source = nomad.shortcut$metadata$source,
      nomad.auto.filled = nomad.shortcut$metadata$auto.filled
    )
    nomad.inner <- .np_nomad_validate_inner_multistart(
      call_names = mc.names,
      dot.args = dots,
      nomad.nmulti = nomad.nmulti,
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      automatic.degree.search = !is.null(degree.search),
      search.engine = if (is.null(degree.search)) "" else degree.search$engine
    )
    nomad.inner.named <- nomad.inner$named
    nomad.inner.nmulti <- nomad.inner$nmulti
    degree.setup <- npSetupGlpDegree(
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      degree = degree,
      ncon = 1L,
      degree.select = degree.select.value
    )
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)
    spec.mc.names <- mc.names
    if (isTRUE(nomad.shortcut$enabled))
      spec.mc.names <- unique(c(spec.mc.names, "regtype", "bernstein.basis"))

    spec <- npResolveCanonicalConditionalRegSpec(
      mc.names = spec.mc.names,
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      basis = basis,
      degree = degree.setup,
      bernstein.basis = if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else bernstein.basis,
      ncon = 1L,
      where = "npindexbw"
    )
    if (!is.null(degree.search)) {
      spec$bernstein.basis <- degree.search$bernstein.basis
      spec$bernstein.basis.engine <- degree.search$bernstein.basis
    }
    initial.bwtype <- if (!is.null(nomad.shortcut$values$bwtype)) {
      as.character(nomad.shortcut$values$bwtype)[1L]
    } else if ("bwtype" %in% dot.names) {
      as.character(dots$bwtype)[1L]
    } else {
      "fixed"
    }
    if (isTRUE(bandwidth.compute) &&
        !is.null(degree.search) &&
        initial.bwtype %in% c("generalized_nn", "adaptive_nn")) {
      h.candidate <- .npindex_nn_candidate_bandwidth(
        h = bws[p + 1L],
        bwtype = initial.bwtype,
        nobs = nrow(xdat)
      )
      bws[p + 1L] <- if (isTRUE(h.candidate$ok)) {
        h.candidate$value
      } else {
        .npindex_default_start_bandwidth(
          fit = NULL,
          bwtype = initial.bwtype,
          nobs = nrow(xdat)
        )
      }
    }
    tbw <- sibandwidth(beta = bws[seq_len(p)],
                       h = bws[p+1L], ...,
                       regtype = spec$regtype,
                       basis = spec$basis,
                       degree = spec$degree,
                       bernstein.basis = spec$bernstein.basis,
                       nobs = dim(xdat)[1],
                       xdati = untangle(xdat),
                       ydati = untangle(data.frame(ydat)),
                       xnames = names(xdat),
                       ynames = deparse(substitute(ydat)),
                       bandwidth = bws[p+1L],
                       bandwidth.compute = bandwidth.compute)
    tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)

    if (tbw$method == "kleinspady" && !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    margs <- c("nmulti", "nomad.remin", "powell.remin", "random.seed", "optim.method", "optim.maxattempts",
               "nomad.opts",
               "optim.reltol", "optim.abstol", "optim.maxit", "only.optimize.beta",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "scale.factor.search.lower")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)
    bwsel.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
    if (any.m) {
      nms <- mc.names[m]
      bwsel.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }

    if (bandwidth.compute && !is.null(degree.search)) {
      reg.args <- list(
        regtype = spec$regtype,
        basis = spec$basis,
        degree = spec$degree,
        bernstein.basis = spec$bernstein.basis,
        regtype.engine = spec$regtype.engine,
        basis.engine = spec$basis.engine,
        degree.engine = spec$degree.engine,
        bernstein.basis.engine = spec$bernstein.basis.engine,
        scale.factor.search.lower = scale.factor.search.lower
      )
      opt.args <- c(
        list(bandwidth.compute = bandwidth.compute),
        bwsel.args[setdiff(names(bwsel.args), c("xdat", "ydat", "bws"))]
      )
      if ("nomad.opts" %in% names(dots))
        opt.args$nomad.opts <- dots$nomad.opts
      opt.args$scale.factor.search.lower <- scale.factor.search.lower

      eval_fun <- function(degree.vec) {
        cell.reg.args <- reg.args
        cell.reg.args$regtype <- "lp"
        cell.reg.args$degree <- as.integer(degree.vec)
        cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
        cell.reg.args$regtype.engine <- "lp"
        cell.reg.args$degree.engine <- as.integer(degree.vec)
        cell.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
        cell.bws <- .npindexbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          template = tbw,
          reg.args = cell.reg.args,
          opt.args = opt.args
        )
        list(
          objective = as.numeric(cell.bws$fval[1L]),
          payload = cell.bws,
          num.feval = if (!is.null(cell.bws$num.feval)) as.numeric(cell.bws$num.feval[1L]) else NA_real_,
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
          objective_name = "fval",
          source = degree.search$source,
          reason = degree.search$reason
        )
      } else {
        search.result <- .npindexbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          template = tbw,
          reg.args = reg.args,
          opt.args = utils::modifyList(opt.args, list(random.seed = random.seed.value)),
          degree.search = degree.search,
          nomad.inner.nmulti = nomad.inner.nmulti,
          nomad.opts = if (is.null(opt.args$nomad.opts)) list() else opt.args$nomad.opts,
          source = degree.search$source,
          reason = degree.search$reason,
          handoff_before_build = identical(degree.search$engine, "nomad+powell"),
          progress_label = .np_degree_search_label(degree.search$engine, degree.search$source)
        )
      }
      tbw <- .npindexbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
    } else if (bandwidth.compute) {
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting single-index bandwidth",
        do.call(npindexbw.sibandwidth, bwsel.args)
      )
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)
  }

npindexbw.sibandwidth <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws,
           bandwidth.compute = TRUE,
           nmulti,
           only.optimize.beta = FALSE,
           optim.abstol = .Machine$double.eps,
           optim.maxattempts = 10,
           optim.maxit = 500,
           optim.method = c("Nelder-Mead", "BFGS", "CG"),
           optim.reltol = sqrt(.Machine$double.eps),
           random.seed = 42,
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           scale.factor.search.lower = NULL,
           ...){

    dots <- list(...)
    npRejectUnsupportedBwsolver(dots, "npindexbw")

    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)
    on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)


    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(ncol(xdat))
    }
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    only.optimize.beta <- npValidateScalarLogical(only.optimize.beta, "only.optimize.beta")
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    h.start.controls <- .npindexbw_h_start_controls(
      scale.factor.init.lower = scale.factor.init.lower,
      scale.factor.init.upper = scale.factor.init.upper,
      scale.factor.init = scale.factor.init,
      scale.factor.search.lower = scale.factor.search.lower,
      where = "npindexbw"
    )
    .npRmpi_require_active_slave_pool(where = "npindexbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (bws$method == "kleinspady" && !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")

      .np_warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    optim.method <- match.arg(optim.method)

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- attr(na.omit(data.frame(xdat,ydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]

    ## convert to numeric
    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat = toMatrix(xdat)
    p <- ncol(xdat)
    beta.idx <- if (p > 1L) seq_len(p - 1L) else integer(0)
    nobs <- nrow(xdat)
    spec <- .npindex_resolve_spec(bws, where = "npindexbw")
    objective.policy <- .npindex_objective_policy(
      bws = bws,
      spec = spec,
      bandwidth.compute = bandwidth.compute,
      where = "npindexbw"
    )
    .npindexbw_check_index_bound_contract(
      bws = bws,
      policy = objective.policy,
      where = "npindexbw"
    )
    objective.spec <- objective.policy$objective.spec
    service.ctx <- .npindexbw_ichimura_lp_service_context(
      bws = bws,
      spec = objective.spec,
      bandwidth.compute = bandwidth.compute,
      comm = 1L,
      service_id = "npindex_ichimura_lp_fixed"
    )
    if (isTRUE(service.ctx$active) && !isTRUE(service.ctx$root))
      return(.npindexbw_ichimura_lp_service_worker_loop(
        xmat = xdat,
        ydat = ydat,
        bws = bws,
        spec = objective.spec,
        ctx = service.ctx
      ))
    ks.service.ctx <- .npindexbw_kleinspady_lp_service_context(
      bws = bws,
      spec = objective.spec,
      policy = objective.policy,
      bandwidth.compute = bandwidth.compute,
      comm = 1L,
      service_id = "npindex_kleinspady_lp_fixed"
    )
    if (isTRUE(ks.service.ctx$active) && !isTRUE(ks.service.ctx$root))
      return(.npindexbw_kleinspady_lp_service_worker_loop(
        xmat = xdat,
        ydat = ydat,
        bws = bws,
        spec = objective.spec,
        ctx = ks.service.ctx
      ))
    service.eval.counter <- 0L
    service.done <- FALSE
    if (isTRUE(service.ctx$active) && isTRUE(service.ctx$root)) {
      on.exit({
        if (!isTRUE(service.done))
          .npindexbw_ichimura_lp_service_error("npindex Ichimura LP fixed-degree service stopped before returning a result", service.ctx)
      }, add = TRUE)
    }
    ks.service.done <- FALSE
    if (isTRUE(ks.service.ctx$active) && isTRUE(ks.service.ctx$root)) {
      on.exit({
        if (!isTRUE(ks.service.done))
          .npindexbw_kleinspady_lp_service_error("npindex Klein-Spady LP fixed-degree service stopped before returning a result", ks.service.ctx)
      }, add = TRUE)
    }

    total.time <-
      system.time({

        if(bandwidth.compute){

          ## Invariant objects used by objective evaluations.
          xmat <- xdat
          beta.coord <- .npindex_beta_coordinate_setup(xmat)
          bandwidth_eval_count <- 0L
          objective.cache.enabled <- npObjectiveCacheEnabled()
          r.nn.cache.surface <- !isTRUE(only.optimize.beta) &&
            identical(bws$type %in% c("generalized_nn", "adaptive_nn"), TRUE)
          r.nn.cache.eligible <- r.nn.cache.surface &&
            objective.cache.enabled
          r.nn.cache <- if (r.nn.cache.surface) {
            .np_r_nn_cache_new(r.nn.cache.eligible, key.length = length(beta.idx) + 1L)
          } else {
            NULL
          }

          bandwidth_progress_step <- function() {
            bandwidth_eval_count <<- bandwidth_eval_count + 1L
            .np_progress_bandwidth_activity_step(done = bandwidth_eval_count)
            invisible(NULL)
          }
          r_nn_cache_lookup <- function(beta, h) {
            if (!is.environment(r.nn.cache) || !isTRUE(r.nn.cache$enabled))
              return(list(hit = FALSE, token = NULL, value = NULL))
            token <- .np_r_nn_cache_param_key(doubles = beta, integers = as.integer(h))
            .np_r_nn_cache_get_token(r.nn.cache, token)
          }
          r_nn_cache_store <- function(token, value, penalty) {
            if (is.finite(value) && value < penalty)
              .np_r_nn_cache_put(r.nn.cache, token, value)
            invisible(NULL)
          }
          fixed.h.lower <- NULL

          ## Note - there are two methods currently implemented, Ichimura's
          ## least squares approach and Klein and Spady's likelihood approach.

          ##We define ichimura's objective function. Since we normalize beta_1
          ##to equal 1, we only worry about beta_2...beta_k in the index
          ##function for these internals, and when computing the leave-one-out
          ##objective function use c(1,beta). However, we do indeed return
          ##c(1,beta) which can be used in the index.model function above.

          ichimuraMaxPenalty <- 10*mean(ydat^2)

          ichimura <- function(param) {
            bandwidth_progress_step()

            ##Define the leave-one-out objective function, sum (y - \hat
            ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
            ## parameters in `param') and then let h denote the kth column.
            beta <- param[beta.idx]
            h <- param[p]
            h.candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bws$type, nobs = nobs)
            if (!h.candidate$ok)
              return(ichimuraMaxPenalty)
            h <- h.candidate$value
            if (!is.null(fixed.h.lower) && h < fixed.h.lower)
              return(ichimuraMaxPenalty)
            cache.hit <- r_nn_cache_lookup(beta, h)
            if (isTRUE(cache.hit$hit)) {
              num.feval.fast.overall <<- num.feval.fast.overall + 1L
              return(cache.hit$value)
            }

            ## Next we define the sum of squared leave-one-out residuals
            sum.squares.leave.one.out <- function(beta, h) {
              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)
              index <- xmat %*% c(1, beta)
              service.eval.counter <<- service.eval.counter + 1L
              objective <- if (isTRUE(service.ctx$active) && isTRUE(service.ctx$root)) {
                .npindexbw_ichimura_lp_service_eval(
                  param = c(beta, h),
                  xmat = xmat,
                  ydat = ydat,
                  bws = bws,
                  spec = objective.spec,
                  ctx = service.ctx,
                  eval_id = service.eval.counter
                )
              } else {
                .npindexbw_eval_ichimura_lp_via_npreg(
                  index = index,
                  ydat = ydat,
                  h = h,
                  bws = bws,
                  spec = objective.spec,
                  invalid.penalty = ichimuraMaxPenalty
                )
              }
              num.feval.fast.overall <<- num.feval.fast.overall +
                as.numeric(objective$num.feval.fast[1L])
              as.numeric(objective$objective[1L])
            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              fv <- sum.squares.leave.one.out(beta,h)
              r_nn_cache_store(cache.hit$token, fv, ichimuraMaxPenalty)
              return(fv)
            } else {
              return(ichimuraMaxPenalty)
            }

          }

          ichimura.nobw <- function(param,h){ return(ichimura(c(param,h))) }

          ## We define ichimura's objective function. Since we normalize beta_1
          ## to equal 1, we only worry about beta_2...beta_k in the index
          ## function for these internals, and when computing the leave-one-out
          ## objective function use c(1,beta). However, we do indeed return
          ## c(1,beta) which can be used in the index.model function above.

          kleinspady <- function(param) {
            bandwidth_progress_step()

            ## Define the leave-one-out objective function, sum (y - \hat
            ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
            ## parameters in `param') and then let h denote the kth column.
            beta <- param[beta.idx]
            h <- param[p]
            h.candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bws$type, nobs = nobs)
            if (!h.candidate$ok)
              return(sqrt(.Machine$double.xmax))
            h <- h.candidate$value
            if (!is.null(fixed.h.lower) && h < fixed.h.lower)
              return(sqrt(.Machine$double.xmax))
            cache.hit <- r_nn_cache_lookup(beta, h)
            if (isTRUE(cache.hit$hit)) {
              num.feval.fast.overall <<- num.feval.fast.overall + 1L
              return(cache.hit$value)
            }

            ## Next we define the sum of logs
            sum.log.leave.one.out <- function(beta, h) {
              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)
              index <- xmat %*% c(1, beta)
              service.eval.counter <<- service.eval.counter + 1L
              objective <- if (isTRUE(ks.service.ctx$active) && isTRUE(ks.service.ctx$root)) {
                collective <- .npindexbw_kleinspady_lp_service_eval(
                  beta = beta,
                  h = h,
                  xmat = xmat,
                  ydat = ydat,
                  bws = bws,
                  spec = objective.spec,
                  ctx = ks.service.ctx,
                  eval_id = service.eval.counter
                )
                list(
                  objective = if (isTRUE(collective$invalid))
                    sqrt(.Machine$double.xmax)
                  else
                    as.numeric(collective$objective[1L]),
                  num.feval.fast = if (.npindexbw_fast_eligible(h = h, bws = bws, eval.index = index)) 1L else 0L
                )
              } else {
                .npindexbw_eval_kleinspady_lp_via_npreg(
                  index = index,
                  ydat = ydat,
                  h = h,
                  bws = bws,
                  spec = objective.spec,
                  invalid.penalty = sqrt(.Machine$double.xmax)
                )
              }
              num.feval.fast.overall <<- num.feval.fast.overall +
                as.numeric(objective$num.feval.fast[1L])
              as.numeric(objective$objective[1L])
            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              fv <- sum.log.leave.one.out(beta,h)
              r_nn_cache_store(cache.hit$token, fv, sqrt(.Machine$double.xmax))
              return(fv)
            } else {
              ## No natural counterpart to var of y here, unlike Ichimura above...
              return(sqrt(.Machine$double.xmax))
            }

          }

          kleinspady.nobw <- function(param,h){ return(kleinspady(c(param,h))) }
          ## Now we implement multistarting

          fval.min <- .Machine$double.xmax
          numimp <- 0
          fval.value <- numeric(nmulti)
          num.feval.fast.overall <- 0L

          if(bws$method == "ichimura"){
            optim.fn <- if(only.optimize.beta) ichimura.nobw else ichimura
            optim.control <- list(abstol=optim.abstol,
                                  reltol=optim.reltol,
                                  maxit=optim.maxit)
          } else if(bws$method == "kleinspady"){
            optim.fn <- if(only.optimize.beta) kleinspady.nobw else  kleinspady
            optim.control <- list(reltol=optim.reltol,maxit=optim.maxit)
          }
          optim.fn.search <- if (only.optimize.beta) {
            function(param, h) optim.fn(beta.coord$to_public(param), h)
          } else {
            function(param) {
              param <- as.double(param)
              if (length(beta.idx))
                param[beta.idx] <- beta.coord$to_public(param[beta.idx])
              optim.fn(param)
            }
          }

          for (i in seq_len(nmulti)) {
            ## We use the nlm command to minimize the objective function using
            ## starting values. Note that since we normalize beta_1=1 here beta
            ## is the k-1 vector containing beta_2...beta_k

            if(i == 1) {

              ## Initial values taken from OLS fit with a constant used for
              ## multistart 1
              ols.fit <- lm(ydat~xdat,x=TRUE)

              if (p != 1L){
                if (setequal(bws$beta[2:p], c(0)))
                  beta <- .npindex_ols_beta_tail(ols.fit)
                else
                  beta = bws$beta[2:p]
              } else { beta = numeric(0) }
              fit <- .npindex_index_from_beta_tail(xmat, beta)
              fixed.h.lower <- if (identical(bws$type, "fixed")) {
                h.start.controls$scale.factor.search.lower * .npindex_start_bandwidth_scale(fit = fit, nobs = nobs)
              } else {
                NULL
              }

              if (bws$bw == 0)
                h <- .npindex_default_start_bandwidth(
                  fit = fit,
                  bwtype = bws$type,
                  nobs = nobs,
                  start.controls = h.start.controls
                )
              else
                h <- tryCatch(
                  .npindex_finalize_bandwidth(
                    h = bws$bw,
                    bwtype = bws$type,
                    nobs = nobs,
                    lower = fixed.h.lower,
                    where = "npindexbw"
                  ),
                  error = function(e) .npindex_default_start_bandwidth(
                    fit = fit,
                    bwtype = bws$type,
                    nobs = nobs,
                    start.controls = h.start.controls
                  )
                )
            } else {
              ## Random initialization used for remaining multistarts

              ols.beta <- .npindex_ols_beta_tail(ols.fit)
              beta.length <- length(ols.beta)
              beta <- runif(beta.length,min=0.5,max=1.5)*ols.beta
              if (!only.optimize.beta)
                h <- .npindex_random_start_bandwidth(
                  fit = fit,
                  bwtype = bws$type,
                  nobs = nobs,
                  start.controls = h.start.controls
                )
            }

            beta.search <- beta.coord$to_search(beta)
            optim.parm <- if(only.optimize.beta) beta.search else c(beta.search,h)

            optim.base.args <- list(
              par = optim.parm,
              fn = optim.fn.search,
              gr = NULL,
              method = optim.method,
              control = optim.control
            )
            if (only.optimize.beta) {
              optim.base.args$h <- h
            }
            suppressWarnings(optim.return <- do.call(optim, optim.base.args))
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              ols.beta <- .npindex_ols_beta_tail(ols.fit)
              beta.length <- length(ols.beta)
              beta <- runif(beta.length,min=0.5,max=1.5)*ols.beta
              if(!only.optimize.beta)
                h <- .npindex_random_start_bandwidth(
                  fit = fit,
                  bwtype = bws$type,
                  nobs = nobs,
                  start.controls = h.start.controls
                )

              if(optim.return$convergence == 1){
                if(optim.control$maxit < (2^32/10))
                  optim.control$maxit <- 10*optim.control$maxit
                else
                  stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))
              }

              if(optim.return$convergence == 10){
                optim.control$reltol <-  10.0*optim.control$reltol
                if(!is.null(optim.control$abstol))
                  optim.control$abstol <-  10.0*optim.control$abstol
              }

              beta.search <- beta.coord$to_search(beta)
              optim.parm <- if(only.optimize.beta) beta.search else c(beta.search,h)
              optim.base.args$par <- optim.parm
              optim.base.args$control <- optim.control
              if (!only.optimize.beta && ("h" %in% names(optim.base.args))) {
                optim.base.args$h <- NULL
              }
              if (only.optimize.beta) {
                optim.base.args$h <- h
              }
              suppressWarnings(optim.return <- do.call(optim, optim.base.args))
            }

            if(optim.return$convergence != 0)
              stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))

            fval.value[i] <- optim.return$value
            if(optim.return$value < fval.min) {
              param <- if(only.optimize.beta) {
                c(beta.coord$to_public(optim.return$par), h)
              } else {
                out <- as.double(optim.return$par)
                if (length(beta.idx))
                  out[beta.idx] <- beta.coord$to_public(out[beta.idx])
                out
              }
              fval.min <- optim.return$value
              numimp <- numimp + 1
              best <- i
            }

            .np_progress_bandwidth_multistart_step(done = i, total = nmulti)

          }

          bws$beta <- c(1.0, param[beta.idx])
          bws$bw <- .npindex_finalize_bandwidth(
            h = param[p],
            bwtype = bws$type,
            nobs = nobs,
            lower = fixed.h.lower,
            where = "npindexbw"
          )
          bws$fval <- fval.min
          bws$ifval <- best
          bws$num.feval <- bandwidth_eval_count
          bws$num.feval.fast <- num.feval.fast.overall
          bws$nn.cache <- .np_r_nn_cache_stats(r.nn.cache)
          bws$numimp <- numimp
          bws$fval.vector <- fval.value
        }
      })[["elapsed"]]
    ## Return a list with beta (we append the restricted value of
    ## beta_1=1), the bandwidth h, the value of the objective function at
    ## its minimum, the number of restarts that resulted in an improved
    ## value of the objective function, the restart that resulted in the
    ## smallest value, and the vector of objective function values.

    ## Restore seed

    .np_seed_exit(seed.state, remove_if_absent = TRUE)
    nn.cache <- bws$nn.cache

    bws <- sibandwidth(beta = bws$beta,
                       h = bws$bw,
                       method = bws$method,
                       regtype = bws$regtype,
                       basis = bws$basis,
                       degree = bws$degree,
                       bernstein.basis = bws$bernstein.basis,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       ckerbound = bws$ckerbound,
                       ckerlb = bws$ckerlb,
                       ckerub = bws$ckerub,
                       bwtype = bws$type,
                       fval = bws$fval,
                       ifval = bws$ifval,
                       num.feval = bws$num.feval,
                       num.feval.fast = bws$num.feval.fast,
                       numimp = bws$numimp,
                       fval.vector = bws$fval.vector,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       bandwidth = bws$bw,
                       rows.omit = rows.omit,
                       bandwidth.compute = bandwidth.compute,
                       optim.method = optim.method,
                       only.optimize.beta = only.optimize.beta,
                       total.time = total.time)
    bws$nn.cache <- nn.cache
    bws <- npSetScaleFactorSearchLower(bws, scale.factor.search.lower)

    if (isTRUE(service.ctx$active) && isTRUE(service.ctx$root)) {
      service.done <- TRUE
      return(.npindexbw_ichimura_lp_service_result(bws, service.ctx))
    }
    if (isTRUE(ks.service.ctx$active) && isTRUE(ks.service.ctx$root)) {
      ks.service.done <- TRUE
      return(.npindexbw_kleinspady_lp_service_result(bws, ks.service.ctx))
    }

    bws

  }
