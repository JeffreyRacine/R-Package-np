nplsqreg <-
  function(bws, ...) {
    args <- list(...)

    if (!missing(bws)) {
      if (isa(bws, "lsqregressionbandwidth"))
        UseMethod("nplsqreg", bws)
      if (inherits(bws, "formula") && is.null(args$txdat))
        UseMethod("nplsqreg", bws)
      if (is.recursive(bws)) {
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("nplsqreg", bws$formula)
        else if (!is.call(bws))
          UseMethod("nplsqreg", bws)
        else
          UseMethod("nplsqreg", NULL)
      } else {
        UseMethod("nplsqreg", NULL)
      }
    } else {
      UseMethod("nplsqreg", NULL)
    }
  }

nplsqregbw <-
  function(bws, ...) {
    args <- list(...)

    if (!missing(bws)) {
      if (isa(bws, "lsqregressionbandwidth"))
        UseMethod("nplsqregbw", bws)
      if (inherits(bws, "formula") && is.null(args$xdat))
        UseMethod("nplsqregbw", bws)
      if (is.recursive(bws) && !is.call(bws))
        UseMethod("nplsqregbw", bws)
      else
        UseMethod("nplsqregbw", NULL)
    } else {
      UseMethod("nplsqregbw", NULL)
    }
  }

.npRmpi_nplsqreg_reset_worker_comm_state <- function(comm = 1L) {
  reset <- get0(".npRmpi_npqreg_reset_worker_comm_state",
                envir = asNamespace("npRmpi"),
                inherits = FALSE)
  if (is.function(reset))
    reset(comm = comm)
}

.npRmpi_nplsqreg_should_autodispatch <- function(mc) {
  call.name <- try(.npRmpi_autodispatch_call_name(mc), silent = TRUE)
  .npRmpi_autodispatch_active() &&
    !isTRUE(.npRmpi_autodispatch_in_context()) &&
    !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
    !inherits(call.name, "try-error") &&
    !identical(call.name, "autodispatch-call")
}

.nplsqreg_validate_tau_values <- function(tau, allow.duplicates = FALSE) {
  if (!is.numeric(tau) || !length(tau) || anyNA(tau) ||
      any(!is.finite(tau)) || any(tau <= 0) || any(tau >= 1))
    stop("'tau' must contain finite numeric values in (0, 1)",
         call. = FALSE)
  tau <- as.numeric(tau)
  if (!isTRUE(allow.duplicates) && anyDuplicated(tau))
    stop("duplicate 'tau' values are not supported in this nplsqreg tranche",
         call. = FALSE)
  tau
}

.nplsqreg_validate_tau <- function(tau) {
  tau <- .nplsqreg_validate_tau_values(tau, allow.duplicates = TRUE)
  if (length(tau) != 1L)
    stop("'tau' must be a single finite numeric value in (0, 1)",
         call. = FALSE)
  tau
}

.nplsqreg_tau_labels <- function(tau) {
  paste0("tau=", format(tau, trim = TRUE, scientific = FALSE))
}

.nplsqreg_omit_length <- function(omit) {
  if (is.null(omit)) 0L else length(omit)
}

.nplsqreg_omit_rows <- function(omit) {
  if (!.nplsqreg_omit_length(omit)) NA_integer_ else as.vector(omit)
}

.nplsqreg_napredict_eval <- function(omit, x) {
  if (!.nplsqreg_omit_length(omit))
    return(x)
  omit <- as.integer(omit)
  keep <- seq_len(NROW(x) + length(omit))[-omit]
  if (is.null(dim(x))) {
    out <- x[rep(NA_integer_, length(x) + length(omit))]
    out[keep] <- x
    attr(out, "na.action") <- NULL
    return(out)
  }
  if (is.data.frame(x)) {
    out <- x[rep(NA_integer_, nrow(x) + length(omit)), , drop = FALSE]
    out[keep, ] <- x
    attr(out, "na.action") <- NULL
    return(out)
  }
  if (length(dim(x)) <= 2L) {
    out <- x[rep(NA_integer_, nrow(x) + length(omit)), , drop = FALSE]
    out[keep, ] <- x
    attr(out, "na.action") <- NULL
    return(out)
  }
  d <- dim(x)
  dn <- dimnames(x)
  new.dim <- c(d[1L] + length(omit), d[-1L])
  if (!is.null(dn)) {
    dn[[1L]] <- NULL
    if (length(dn) != length(new.dim) ||
        any(vapply(seq_along(dn), function(i) {
          !is.null(dn[[i]]) && length(dn[[i]]) != new.dim[[i]]
        }, logical(1L))))
      dn <- NULL
  }
  out <- array(NA_real_, dim = new.dim, dimnames = dn)
  out[keep, , ] <- x
  out
}

.nplsqreg_prepare_train_data <- function(xdat, ydat, scale = NULL) {
  xdat <- toFrame(xdat)
  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector")
  if (nrow(xdat) != length(ydat))
    stop("number of explanatory data and response data do not match")
  if (is.factor(ydat))
    stop("nplsqreg requires a numeric dependent variable", call. = FALSE)
  ydat <- as.numeric(ydat)

  rows.omit <- attr(na.omit(data.frame(xdat, .nplsqreg_y = ydat)),
                    "na.action")
  if (.nplsqreg_omit_length(rows.omit)) {
    keep.rows <- rep(TRUE, nrow(xdat))
    keep.rows[as.integer(rows.omit)] <- FALSE
    xdat <- xdat[keep.rows, , drop = FALSE]
    ydat <- ydat[keep.rows]
    if (!is.null(scale) && length(scale) == length(keep.rows))
      scale <- scale[keep.rows]
  }
  if (any(!is.finite(ydat)))
    stop("'ydat' must be finite")
  list(xdat = xdat, ydat = ydat, scale = scale, omit = rows.omit)
}

.nplsqreg_prepare_eval_data <- function(exdat) {
  exdat <- toFrame(exdat)
  rows.omit <- attr(na.omit(exdat), "na.action")
  if (.nplsqreg_omit_length(rows.omit)) {
    keep.rows <- rep(TRUE, nrow(exdat))
    keep.rows[as.integer(rows.omit)] <- FALSE
    exdat <- exdat[keep.rows, , drop = FALSE]
  }
  list(exdat = exdat, omit = rows.omit)
}

.nplsqreg_same_training_data <- function(xdat, ydat, bws) {
  xdat <- toFrame(xdat)
  bx <- toFrame(bws$xdat)
  row.names(xdat) <- NULL
  row.names(bx) <- NULL
  identical(xdat, bx) && identical(as.numeric(ydat), as.numeric(bws$ydat))
}

.nplsqreg_assert_reuse_training_data <- function(bws, txdat, tydat) {
  if (!.nplsqreg_same_training_data(txdat, tydat, bws))
    stop("explicit txdat/tydat do not match the stored nplsqreg bandwidth training data",
         call. = FALSE)
  invisible(TRUE)
}

.nplsqreg_record_omit <- function(obj, omit) {
  if (.nplsqreg_omit_length(omit)) {
    obj$omit <- omit
    obj$rows.omit <- as.vector(omit)
    obj$nobs.omit <- length(omit)
  } else {
    obj$rows.omit <- NA_integer_
    obj$nobs.omit <- 0L
  }
  obj
}

.nplsqreg_pad_fit_outputs <- function(out, omit) {
  if (!.nplsqreg_omit_length(omit))
    return(out)
  out$quantile <- .nplsqreg_napredict_eval(omit, out$quantile)
  out$quanterr <- .nplsqreg_napredict_eval(omit, out$quanterr)
  if (isTRUE(out$gradients)) {
    out$quantgrad <- .nplsqreg_napredict_eval(omit, out$quantgrad)
    out$quantgerr <- .nplsqreg_napredict_eval(omit, out$quantgerr)
  }
  out$xeval <- .nplsqreg_napredict_eval(omit, out$xeval)
  out$nobs <- nrow(out$xeval)
  out
}

.nplsqreg_attach_native_diagnostics <- function(obj, diagnostics) {
  if (!is.null(diagnostics))
    attr(obj, "native.nomad.diagnostics") <- diagnostics
  obj
}

.nplsqreg_native_diagnostics_from_search <- function(search.result) {
  if (!isTRUE(getOption("np.developer.native.nomad.diagnostics", FALSE)) ||
      is.null(search.result) ||
      is.null(search.result$native.diagnostics))
    return(NULL)
  search.result$native.diagnostics
}

.nplsqreg_collect_native_diagnostics <- function(objects, labels) {
  child <- lapply(objects, function(z) attr(z, "native.nomad.diagnostics", exact = TRUE))
  names(child) <- labels
  keep <- vapply(child, is.list, logical(1L))
  if (!any(keep))
    return(NULL)
  child <- child[keep]
  list(
    route_native = all(vapply(child, function(z) isTRUE(z$route_native), logical(1L))),
    source = "nplsqreg-tau-composition",
    children = child
  )
}

.nplsqreg_validate_tau_search <- function(tau.search) {
  match.arg(tau.search, c("full", "refined"))
}

.nplsqreg_check_loss <- function(u, tau) {
  mean(u * (tau - (u < 0)), na.rm = TRUE)
}

.nplsqreg_validate_scale <- function(scale, n) {
  if (!is.numeric(scale) || length(scale) != n ||
      any(!is.finite(scale)) || any(scale <= 0))
    stop("'scale' must be a finite positive numeric vector matching the training data length",
         call. = FALSE)
  as.numeric(scale)
}

.nplsqreg_scale_pilot_fit <- function(xdat, ydat, pilot.dots,
                                      regtype.pilot, nomad.pilot) {
  mean.bw <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat), pilot.dots))
  mean.fit <- npreg(bws = mean.bw, txdat = xdat, tydat = ydat)
  res <- as.numeric(ydat) - as.numeric(fitted(mean.fit))
  scale.fit <- npreg(bws = mean.fit$bws, txdat = xdat, tydat = res^2)
  variance.hat <- as.numeric(fitted(scale.fit))
  floor <- .Machine$double.eps
  variance.hat <- pmax(floor, variance.hat)
  variance.hat[!is.finite(variance.hat)] <- floor
  scale <- sqrt(variance.hat)
  list(scale = scale, mean.fit = mean.fit, scale.fit = scale.fit,
       regtype.pilot = regtype.pilot, nomad.pilot = nomad.pilot)
}

.nplsqreg_scale_pilot <- function(xdat, ydat, dots,
                                  regtype.pilot = c("auto", "ll", "lc", "lp"),
                                  nomad.pilot = FALSE,
                                  pilot.args = list()) {
  xdat <- toFrame(xdat)
  regtype.pilot <- match.arg(regtype.pilot)
  nomad.pilot <- npValidateScalarLogical(nomad.pilot, "nomad.pilot")
  if (!is.list(pilot.args))
    stop("'pilot.args' must be a list", call. = FALSE)

  ncon <- sum(vapply(xdat, function(z) inherits(z, c("integer", "numeric")),
                     logical(1)))
  if (isTRUE(nomad.pilot)) {
    if (!(regtype.pilot %in% c("auto", "lp")))
      stop("nomad.pilot=TRUE requires regtype.pilot='auto' or 'lp'",
           call. = FALSE)
    regtype.pilot <- "lp"
  } else if (identical(regtype.pilot, "auto")) {
    regtype.pilot <- if (ncon > 0L) "ll" else "lc"
  }

  pilot.dots <- utils::modifyList(dots, pilot.args)
  pilot.dots$regtype <- regtype.pilot
  if (isTRUE(nomad.pilot))
    pilot.dots$nomad <- TRUE

  pilot.bwtype <- if (is.null(pilot.dots$bwtype)) {
    "fixed"
  } else {
    as.character(pilot.dots$bwtype[1L])
  }
  mpi.size <- tryCatch(as.integer(mpi.comm.size(1L)), error = function(e) 1L)

  if (identical(pilot.bwtype, "generalized_nn") &&
      !identical(regtype.pilot, "lc") &&
      isTRUE(getOption("npRmpi.autodispatch.context", FALSE)) &&
      is.finite(mpi.size) &&
      mpi.size > 1L) {
    rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
    if (isTRUE(rank == 0L)) {
      out <- .npRmpi_with_local_regression(
        .nplsqreg_scale_pilot_fit(
          xdat = xdat,
          ydat = ydat,
          pilot.dots = pilot.dots,
          regtype.pilot = regtype.pilot,
          nomad.pilot = nomad.pilot
        )
      )
      mpi.bcast.Robj(out, rank = 0L, comm = 1L)
      return(out)
    }
    return(mpi.bcast.Robj(rank = 0L, comm = 1L))
  }

  .nplsqreg_scale_pilot_fit(
    xdat = xdat,
    ydat = ydat,
    pilot.dots = pilot.dots,
    regtype.pilot = regtype.pilot,
    nomad.pilot = nomad.pilot
  )
}

.nplsqreg_rebuild_rbandwidth <- function(template, bw, ydat, xdat = NULL,
                                         fval = NA, num.feval = NA,
                                         bandwidth.compute = FALSE) {
  if (!is.null(xdat)) {
    tbw <- template
    tbw$bw <- as.numeric(bw)
    out <- npregbw.rbandwidth(xdat = xdat, ydat = ydat, bws = tbw,
                              bandwidth.compute = FALSE)
    out$fval <- fval
    out$num.feval <- num.feval
    out
  } else {
  out <- rbandwidth(
    bw = as.numeric(bw),
    regtype = template$regtype,
    basis = template$basis,
    degree = template$degree,
    bernstein.basis = template$bernstein.basis,
    bwmethod = if (is.null(template$method)) "cv.ls" else template$method,
    bwscaling = template$scaling,
    bwtype = template$type,
    ckertype = template$ckertype,
    ckerorder = template$ckerorder,
    ckerbound = template$ckerbound,
    ckerlb = template$ckerlb,
    ckerub = template$ckerub,
    ukertype = template$ukertype,
    okertype = template$okertype,
    fval = fval,
    num.feval = num.feval,
    nobs = template$nobs,
    xdati = template$xdati,
    ydati = untangle(data.frame(ydat)),
    xnames = template$xnames,
    ynames = template$ynames,
    sfactor = as.numeric(bw),
    bandwidth = as.numeric(bw),
    rows.omit = template$rows.omit,
    nconfac = template$nconfac,
    ncatfac = template$ncatfac,
    sdev = template$sdev,
    bandwidth.compute = bandwidth.compute,
    timing = template$timing,
    total.time = template$total.time)
  out
  }
}

.nplsqreg_optimizer_controls <- function(dots, optim.control) {
  itmax <- if (!is.null(optim.control$maxit)) {
    npValidatePositiveInteger(optim.control$maxit, "optim.control$maxit")
  } else {
    10000L
  }
  if (!is.null(dots$itmax))
    itmax <- npValidatePositiveInteger(dots$itmax, "itmax")

  out <- list(
    nmulti = if (!is.null(dots$nmulti)) npValidateNmulti(dots$nmulti) else 1L,
    itmax = itmax,
    ftol = if (!is.null(dots$ftol)) npValidatePositiveFiniteNumeric(dots$ftol, "ftol") else 1.490116e-07,
    tol = if (!is.null(dots$tol)) npValidatePositiveFiniteNumeric(dots$tol, "tol") else 1.490116e-04,
    small = if (!is.null(dots$small)) npValidatePositiveFiniteNumeric(dots$small, "small") else 1.490116e-05,
    powell.remin = if (!is.null(dots$powell.remin)) npValidateScalarLogical(dots$powell.remin, "powell.remin") else TRUE
  )
  out
}

.nplsqreg_strip_optimizer_dots <- function(dots) {
  optimizer.names <- c(
    "nmulti", "itmax", "ftol", "tol", "small", "powell.remin",
    "nomad", "nomad.nmulti", "nomad.remin", "search.engine",
    "degree.select", "degree.min", "degree.max", "degree.start",
    "degree.restarts", "degree.max.cycles", "degree.verify",
    "optim.control", "delta.bounds", "tau.search",
    "random.seed", "nomad.opts"
  )
  dots[setdiff(names(dots), optimizer.names)]
}

.nplsqreg_normalize_dots <- function(dots, where = "nplsqreg") {
  if (!length(dots))
    return(dots)
  npRejectUnsupportedBwsolver(dots, where)
  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep("", length(dots))
  if (any(!nzchar(dot.names)))
    stop(sprintf("%s requires named optional arguments", where), call. = FALSE)

  x.aliases <- c(cxkertype = "ckertype",
                 cxkerorder = "ckerorder",
                 cxkerbound = "ckerbound",
                 cxkerlb = "ckerlb",
                 cxkerub = "ckerub",
                 uxkertype = "ukertype",
                 oxkertype = "okertype")
  y.aliases <- c("cykertype", "cykerorder", "cykerbound", "cykerlb",
                 "cykerub", "uykertype", "oykertype")

  y.hits <- intersect(dot.names, y.aliases)
  if (length(y.hits))
    stop(sprintf(
      "%s does not use response-side conditional-distribution kernel controls: %s",
      where, paste(y.hits, collapse = ", ")
    ), call. = FALSE)

  for (alias in intersect(names(x.aliases), dot.names)) {
    target <- x.aliases[[alias]]
    if (target %in% dot.names && !identical(dots[[alias]], dots[[target]]))
      stop(sprintf("conflicting '%s' and '%s' arguments", alias, target),
           call. = FALSE)
    dots[[target]] <- dots[[alias]]
    dots[[alias]] <- NULL
    dot.names <- names(dots)
  }

  optimizer.names <- c(
    "nmulti", "itmax", "ftol", "tol", "small", "powell.remin",
    "nomad", "nomad.nmulti", "nomad.remin", "search.engine",
    "degree.select", "degree.min", "degree.max", "degree.start",
    "degree.restarts", "degree.max.cycles", "degree.verify",
    "optim.control", "delta.bounds", "tau.search",
    "random.seed", "nomad.opts"
  )
  reg.names <- setdiff(names(formals(utils::getS3method("npregbw", "default"))),
                       c("xdat", "ydat", "bws", "bandwidth.compute", "..."))
  eval.names <- c("newdata", "exdat", "data")
  allowed <- unique(c(reg.names, optimizer.names))
  bad <- setdiff(names(dots), allowed)
  if (length(bad)) {
    if (any(bad %in% eval.names))
      stop(sprintf(
        "%s evaluation argument '%s' reached bandwidth selection unexpectedly",
        where, bad[bad %in% eval.names][1L]
      ), call. = FALSE)
    stop(sprintf("unused %s argument '%s'", where, bad[1L]), call. = FALSE)
  }
  dots
}

.nplsqreg_call_fixed_degree_core <- function(xdat, ydat, scale, tau, bws,
                                             delta, delta.bounds, opt.args,
                                             bandwidth.compute) {
  opt.value <- function(name, default)
    if (is.null(opt.args[[name]])) default else opt.args[[name]]

  xdat <- toFrame(xdat)
  xmat <- toMatrix(xdat)
  ydat <- as.double(ydat)
  scale <- as.double(scale)
  if (length(scale) != length(ydat))
    stop("'scale' must have the same length as 'ydat'", call. = FALSE)

  nmulti <- if (isTRUE(bandwidth.compute)) {
    opt.value("nmulti", npDefaultNmulti(dim(xdat)[2L]))
  } else {
    1L
  }
  nmulti <- npValidateNmulti(nmulti)
  itmax <- npValidatePositiveInteger(opt.value("itmax", 10000L), "itmax")
  remin <- npValidateScalarLogical(opt.value("powell.remin", TRUE), "powell.remin")
  scale.init.categorical.sample <-
    npValidateScalarLogical(opt.value("scale.init.categorical.sample", FALSE),
                            "scale.init.categorical.sample")
  transform.bounds <-
    npValidateScalarLogical(opt.value("transform.bounds", FALSE),
                            "transform.bounds")
  ftol <- npValidatePositiveFiniteNumeric(opt.value("ftol", 1.490116e-07), "ftol")
  tol <- npValidatePositiveFiniteNumeric(opt.value("tol", 1.490116e-04), "tol")
  small <- npValidatePositiveFiniteNumeric(opt.value("small", 1.490116e-05), "small")
  penalty.multiplier <-
    npValidatePositiveFiniteNumeric(opt.value("penalty.multiplier", 10),
                                    "penalty.multiplier")
  invalid.penalty <- match.arg(opt.value("invalid.penalty", "baseline"),
                               c("baseline", "dbmax"))
  penalty.mode <- if (invalid.penalty == "baseline") 1L else 0L
  scale.factor.search.lower <- npResolveScaleFactorLowerBound(
    if (is.null(opt.args$scale.factor.search.lower)) {
      npGetScaleFactorSearchLower(bws)
    } else {
      opt.args$scale.factor.search.lower
    }
  )
  .np_progress_bandwidth_set_total(nmulti)

  runo <- xmat[, bws$iuno, drop = FALSE]
  rcon <- xmat[, bws$icon, drop = FALSE]
  rord <- xmat[, bws$iord, drop = FALSE]
  mysd <- EssDee(rcon)
  nrow <- dim(xmat)[1L]
  nconfac <- nrow^(-1.0 / (2.0 * bws$ckerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$ckerorder + bws$ncon))

  reg.spec <- npCanonicalConditionalRegSpec(
    regtype = bws$regtype,
    basis = bws$basis,
    degree = bws$degree,
    bernstein.basis = bws$bernstein.basis,
    ncon = bws$ncon,
    where = "nplsqregbw"
  )
  reg.c <- npRegtypeToC(regtype = reg.spec$regtype.engine,
                        degree = reg.spec$degree.engine,
                        ncon = bws$ncon,
                        context = "nplsqregbw")
  npCheckRegressionDesignCondition(reg.code = reg.c$code,
                                   xcon = rcon,
                                   basis = reg.spec$basis.engine,
                                   degree = reg.spec$degree.engine,
                                   bernstein.basis = reg.spec$bernstein.basis.engine,
                                   where = "nplsqregbw")
  degree.c <- if (bws$ncon > 0L) {
    as.integer(if (is.null(reg.c$degree)) rep.int(0L, bws$ncon) else reg.c$degree)
  } else {
    integer(1L)
  }
  cont.start <- npContinuousSearchStartControls(
    opt.value("scale.factor.init.lower", 0.1),
    opt.value("scale.factor.init.upper", 2.0),
    opt.value("scale.factor.init", 0.5),
    scale.factor.search.lower,
    where = "nplsqregbw"
  )
  myopti <- list(
    num_obs_train = nrow,
    iMultistart = if (isTRUE(bandwidth.compute)) IMULTI_TRUE else IMULTI_FALSE,
    iNum_Multistart = nmulti,
    int_use_starting_values = if (all(bws$bw == 0)) USE_START_NO else USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_reg_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = itmax,
    int_RESTART_FROM_MIN = if (remin) RE_MIN_TRUE else RE_MIN_FALSE,
    int_MINIMIZE_IO = if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
    bwmethod = BWM_CVLS,
    kerneval = switch(bws$ckertype,
      gaussian = CKER_GAUSS + bws$ckerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$ckerorder / 2 - 1,
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
    int_do_tree = .npregbw_tree_code(bws, ncon = bws$ncon,
                                     ncat = bws$nuno + bws$nord),
    scale.init.categorical.sample = scale.init.categorical.sample,
    dfc.dir = opt.value("dfc.dir", 3),
    transform.bounds = transform.bounds
  )
  myoptd <- list(
    ftol = ftol,
    tol = tol,
    small = small,
    lbc.dir = opt.value("lbc.dir", 0.5),
    cfac.dir = opt.value("cfac.dir", 2.5 * (3.0 - sqrt(5))),
    initc.dir = opt.value("initc.dir", 1.0),
    lbd.dir = opt.value("lbd.dir", 0.1),
    hbd.dir = opt.value("hbd.dir", 1),
    dfac.dir = opt.value("dfac.dir", 0.25 * (3.0 - sqrt(5))),
    initd.dir = opt.value("initd.dir", 1.0),
    lbc.init = cont.start$scale.factor.init.lower,
    hbc.init = cont.start$scale.factor.init.upper,
    cfac.init = cont.start$scale.factor.init,
    lbd.init = opt.value("lbd.init", 0.1),
    hbd.init = opt.value("hbd.init", 0.9),
    dfac.init = opt.value("dfac.init", 0.375),
    nconfac = nconfac,
    ncatfac = ncatfac,
    scale.factor.lower.bound = scale.factor.search.lower
  )
  cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon],
                                         bws$ckerub[bws$icon])
  out <- if (isTRUE(bandwidth.compute)) {
    .Call(
      "C_np_lsqregression_bw",
      as.double(runo), as.double(rord), as.double(rcon),
      as.double(ydat), as.double(scale),
      as.double(tau), as.double(delta), as.double(delta.bounds),
      as.double(mysd), as.integer(myopti), as.double(myoptd),
      as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
      as.integer(nmulti),
      as.integer(penalty.mode), as.double(penalty.multiplier),
      as.integer(degree.c),
      as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
      as.integer(npLpBasisCode(reg.spec$basis.engine)),
      as.double(cker.bounds.c$lb), as.double(cker.bounds.c$ub),
      PACKAGE = "npRmpi"
    )
  } else {
    .Call(
      "C_np_lsqregression_bw_eval",
      as.double(runo), as.double(rord), as.double(rcon),
      as.double(ydat), as.double(scale),
      as.double(tau), as.double(delta), as.double(delta.bounds),
      as.double(mysd), as.integer(myopti), as.double(myoptd),
      as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
      as.integer(nmulti),
      as.integer(penalty.mode), as.double(penalty.multiplier),
      as.integer(degree.c),
      as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
      as.integer(npLpBasisCode(reg.spec$basis.engine)),
      as.double(cker.bounds.c$lb), as.double(cker.bounds.c$ub),
      PACKAGE = "npRmpi"
    )
  }
  rorder <- numeric(length(bws$bw))
  ord.idx <- seq_along(bws$bw)
  rorder[c(ord.idx[bws$icon], ord.idx[bws$iuno], ord.idx[bws$iord])] <- ord.idx
  list(
    bw = as.numeric(out$bw[rorder]),
    delta = as.numeric(out$delta[1L]),
    objective = as.numeric(out$fval[1L]),
    ifval = as.numeric(out$fval[2L]),
    num.feval = sum(out$eval.history[is.finite(out$eval.history)]),
    num.feval.fast = as.numeric(out$fast.history[1L]),
    nn.cache = .np_nn_cache_stats(out$nn.cache),
    fval.history = out$fval.history,
    eval.history = out$eval.history,
    invalid.history = out$invalid.history,
    timing = out$timing
  )
}

.nplsqreg_nomad_search <- function(xdat, ydat, scale, tau, template,
                                   delta, delta.bounds, opt.args,
                                   degree.search,
                                   nomad.inner.nmulti = 0L,
                                   random.seed = 42L,
                                   nomad.opts = list()) {
  if (isTRUE(degree.search$verify))
    stop("nplsqregbw nomad=TRUE does not support degree.verify=TRUE",
         call. = FALSE)

  template$regtype <- "lp"
  template$degree <- as.integer(degree.search$start.degree)
  template$bernstein.basis <- degree.search$bernstein.basis
  if (!(template$type %in% c("fixed", "generalized_nn", "adaptive_nn")))
    stop("nplsqregbw nomad=TRUE requires bwtype='fixed', 'generalized_nn', or 'adaptive_nn'",
         call. = FALSE)

  setup <- .npregbw_nomad_bw_setup(xdat = xdat, template = template, allow.extended.nn = TRUE)
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) {
    npDefaultNmulti(dim(toFrame(xdat))[2L])
  } else {
    npValidateNmulti(opt.args$nmulti[1L])
  }
  opt.value <- function(name, default) {
    if (is.null(opt.args[[name]])) default else opt.args[[name]]
  }
  bw.bounds <- .npregbw_nomad_bw_bounds(template = template, setup = setup)
  bw.start.bounds <- .np_nomad_bw_restart_start_bounds(
    bounds = bw.bounds,
    setup = setup,
    opt.value = opt.value,
    where = "nplsqregbw"
  )
  point.start <- if (all(template$bw == 0)) {
    NULL
  } else {
    .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup)
  }
  bw.start <- .npregbw_nomad_complete_bw_start_point(
    point = point.start,
    bounds = bw.bounds,
    setup = setup
  )
  x0 <- c(bw.start, delta, as.integer(degree.search$start.degree))
  lb <- c(bw.bounds$lower, delta.bounds[1L], degree.search$lower)
  ub <- c(bw.bounds$upper, delta.bounds[2L], degree.search$upper)
  bbin <- c(bw.bounds$bbin, 0L, rep.int(1L, ncon))
  coordinate.roles <- c(
    .np_nomad_coordinate_roles(bw.bounds),
    "continuous_real",
    rep.int("degree", length(degree.search$lower))
  )
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  eval_fun <- function(point) {
    point <- as.numeric(point)
    bw.idx <- seq_len(ncon + ncat)
    delta.idx <- ncon + ncat + 1L
    deg.idx <- delta.idx + seq_len(ncon)
    degree <- as.integer(round(point[deg.idx]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw.vec <- .npregbw_nomad_point_to_bw(point[bw.idx],
                                         template = template,
                                         setup = setup)
    tbw <- template
    tbw$bw <- bw.vec
    tbw$regtype <- "lp"
    tbw$degree <- degree
    tbw$bernstein.basis <- degree.search$bernstein.basis
    out <- .nplsqreg_call_fixed_degree_core(
      xdat = xdat,
      ydat = ydat,
      scale = scale,
      tau = tau,
      bws = tbw,
      delta = point[delta.idx],
      delta.bounds = delta.bounds,
      opt.args = utils::modifyList(opt.args, list(nmulti = 1L)),
      bandwidth.compute = FALSE
    )
    nomad.num.feval.total <<- nomad.num.feval.total + as.numeric(out$num.feval[1L])
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])
    list(objective = out$objective, degree = degree, num.feval = out$num.feval)
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    bw.idx <- seq_len(ncon + ncat)
    delta.idx <- ncon + ncat + 1L
    degree <- as.integer(best_record$degree)
    bw.vec <- .npregbw_nomad_point_to_bw(point[bw.idx],
                                         template = template,
                                         setup = setup)
    tbw <- template
    tbw$bw <- bw.vec
    tbw$regtype <- "lp"
    tbw$degree <- degree
    tbw$bernstein.basis <- degree.search$bernstein.basis
    direct <- .nplsqreg_call_fixed_degree_core(
      xdat = xdat,
      ydat = ydat,
      scale = scale,
      tau = tau,
      bws = tbw,
      delta = point[delta.idx],
      delta.bounds = delta.bounds,
      opt.args = utils::modifyList(opt.args, list(nmulti = 1L)),
      bandwidth.compute = FALSE
    )
    direct$num.feval <- as.numeric(nomad.num.feval.total)
    direct$num.feval.fast <- as.numeric(nomad.num.feval.fast.total)
    direct$bws <- tbw
    direct$bws$bw <- direct$bw
    direct$bws$degree <- degree
    direct$bws$regtype <- "lp"
    direct$bws$bernstein.basis <- degree.search$bernstein.basis
    direct.objective <- as.numeric(direct$objective[1L])
    powell.elapsed <- NA_real_

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.opt.args <- .np_nomad_powell_hotstart_opt_args(
        opt.args,
        strategy = "disable_multistart",
        remin = isTRUE(opt.args$powell.remin)
      )
      powell.start <- proc.time()[3L]
      hot <- .npregbw_with_powell_refinement_progress(degree, local({
        .nplsqreg_call_fixed_degree_core(
          xdat = xdat,
          ydat = ydat,
          scale = scale,
          tau = tau,
          bws = direct$bws,
          delta = direct$delta,
          delta.bounds = delta.bounds,
          opt.args = hot.opt.args,
          bandwidth.compute = TRUE
        )
      }))
      powell.elapsed <- proc.time()[3L] - powell.start
      hot$bws <- direct$bws
      hot$bws$bw <- hot$bw
      hot$bws$degree <- degree
      hot$bws$regtype <- "lp"
      hot$bws$bernstein.basis <- degree.search$bernstein.basis
      hot$num.feval <- as.numeric(direct$num.feval[1L]) + as.numeric(hot$num.feval[1L])
      hot$num.feval.fast <- as.numeric(direct$num.feval.fast[1L]) + as.numeric(hot$num.feval.fast[1L])
      if (is.finite(hot$objective) &&
          .np_degree_better(hot$objective, direct.objective, direction = "min")) {
        return(list(payload = hot, objective = hot$objective,
                    powell.time = powell.elapsed))
      }
    }
    list(payload = direct, objective = direct.objective,
         powell.time = powell.elapsed)
  }

  search.engine.used <- if (identical(degree.search$engine, "nomad+powell")) {
    "nomad"
  } else {
    degree.search$engine
  }
  result <- .np_nomad_search(
    engine = search.engine.used,
    baseline_record = NULL,
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
    nomad.inner.nmulti = npValidateNonNegativeInteger(nomad.inner.nmulti,
                                                      "nomad.nmulti"),
    random.seed = random.seed,
    handoff_before_build = identical(degree.search$engine, "nomad+powell"),
    remin = isTRUE(opt.args$nomad.remin),
    nomad.opts = nomad.opts,
    native.r.bridge = TRUE,
    start.lower = c(bw.start.bounds$lower, delta.bounds[1L], degree.search$lower),
    start.upper = c(bw.start.bounds$upper, delta.bounds[2L], degree.search$upper),
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
  if (!identical(search.engine.used, degree.search$engine))
    result$method <- degree.search$engine
  result
}

.nplsqreg_combine_bandwidths <- function(bw.list, tau, tau.search = "full",
                                         fit.order = seq_along(tau),
                                         warm.start.from = rep(NA_integer_, length(tau)),
                                         warm.start.degree = vector("list", length(tau)),
                                         tau.search.controls = NULL,
                                         call = NULL) {
  if (!length(bw.list))
    stop("internal error: no scalar nplsqreg bandwidths to combine",
         call. = FALSE)
  tau <- .nplsqreg_validate_tau_values(tau)
  if (length(bw.list) != length(tau))
    stop("internal error: tau and bandwidth-list lengths differ",
         call. = FALSE)
  labels <- .nplsqreg_tau_labels(tau)
  first <- bw.list[[1L]]
  out <- first
  out$bw <- do.call(cbind, lapply(bw.list, `[[`, "bw"))
  colnames(out$bw) <- labels
  out$tau <- tau
  out$delta <- stats::setNames(vapply(bw.list, `[[`, numeric(1L), "delta"),
                              labels)
  out$objective <- stats::setNames(
    vapply(bw.list, function(z) as.numeric(z$objective[1L]), numeric(1L)),
    labels)
  out$qdat <- do.call(cbind, lapply(bw.list, `[[`, "qdat"))
  colnames(out$qdat) <- labels
  out$tau.bws <- stats::setNames(bw.list, labels)
  out$tau.search <- tau.search
  out$fit.order <- as.integer(fit.order)
  out$warm.start.from <- as.integer(warm.start.from)
  out$warm.start.degree <- warm.start.degree
  out$tau.search.controls <- tau.search.controls
  out$pilot.shared <- TRUE
  out$child.degree.common <- all(vapply(
    bw.list,
    function(z) identical(z$reg.bws$degree, first$reg.bws$degree),
    logical(1L)))
  out$total.time <- .nplsqreg_sum_times(lapply(bw.list, function(z) z$total.time))
  out$optim.time <- out$total.time
  out$fit.time <- NA_real_
  out$call <- call
  class(out) <- "lsqregressionbandwidth"
  out <- .nplsqreg_attach_native_diagnostics(
    out,
    .nplsqreg_collect_native_diagnostics(bw.list, labels)
  )
  out
}

.nplsqreg_combine_fits <- function(fit.list, tau, bws, tau.search = "full",
                                   call = NULL) {
  if (!length(fit.list))
    stop("internal error: no scalar nplsqreg fits to combine", call. = FALSE)
  tau <- .nplsqreg_validate_tau_values(tau)
  if (length(fit.list) != length(tau))
    stop("internal error: tau and fit-list lengths differ", call. = FALSE)
  labels <- .nplsqreg_tau_labels(tau)
  first <- fit.list[[1L]]
  out <- first
  out$bw <- bws$bw
  out$bws <- bws
  out$reg.bws <- stats::setNames(lapply(fit.list, `[[`, "reg.bws"), labels)
  out$fit <- stats::setNames(lapply(fit.list, `[[`, "fit"), labels)
  out$tau <- tau
  out$delta <- bws$delta
  out$quantile <- do.call(cbind, lapply(fit.list, `[[`, "quantile"))
  colnames(out$quantile) <- labels
  out$quanterr <- do.call(cbind, lapply(fit.list, `[[`, "quanterr"))
  colnames(out$quanterr) <- labels
  if (isTRUE(first$gradients)) {
    p <- ncol(fit.list[[1L]]$quantgrad)
    grad.names <- colnames(fit.list[[1L]]$quantgrad)
    out$quantgrad <- array(NA_real_,
                           dim = c(nrow(fit.list[[1L]]$quantgrad), p, length(tau)),
                           dimnames = list(NULL, grad.names, labels))
    out$quantgerr <- array(NA_real_,
                           dim = c(nrow(fit.list[[1L]]$quantgerr), p, length(tau)),
                           dimnames = list(NULL, grad.names, labels))
    for (j in seq_along(tau)) {
      out$quantgrad[, , j] <- fit.list[[j]]$quantgrad
      out$quantgerr[, , j] <- fit.list[[j]]$quantgerr
    }
  } else {
    out$quantgrad <- NA
    out$quantgerr <- NA
  }
  out$objective <- bws$objective
  out$optim <- stats::setNames(lapply(fit.list, `[[`, "optim"), labels)
  out$tau.fits <- stats::setNames(fit.list, labels)
  out$tau.search <- tau.search
  out$fit.order <- bws$fit.order
  out$warm.start.from <- bws$warm.start.from
  out$warm.start.degree <- bws$warm.start.degree
  out$tau.search.controls <- bws$tau.search.controls
  out$pilot.shared <- TRUE
  out$child.degree.common <- isTRUE(bws$child.degree.common)
  optim.time <- .nplsqreg_sum_times(bws$total.time)
  fit.time <- .nplsqreg_sum_times(lapply(fit.list, function(z) z$fit.time))
  out$optim.time <- optim.time
  out$fit.time <- fit.time
  out$total.time <- .nplsqreg_sum_times(optim.time, fit.time)
  out$call <- call
  class(out) <- "lsqregression"
  out <- .nplsqreg_attach_native_diagnostics(
    out,
    attr(bws, "native.nomad.diagnostics", exact = TRUE)
  )
  out
}

nplsqregbw.formula <-
  function(bws, data = NULL, tau = 0.5, subset, na.action, ...) {
    .npRmpi_require_active_slave_pool(where = "nplsqregbw()")

    tt <- terms(bws)
    mc <- match.call(expand.dots = FALSE)
    m <- match(c("bws", "data", "subset", "na.action"),
               names(mc), nomatch = 0)
    tmf <- mc[c(1, m)]
    if ("bws" %in% names(tmf))
      names(tmf)[names(tmf) == "bws"] <- "formula"
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    mf <- do.call(stats::model.frame, mf.args, envir = environment(tt))
    train.omit <- attr(mf, "na.action")
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"), "term.labels"), drop = FALSE]
    out <- nplsqregbw(xdat = xdat, ydat = ydat, tau = tau, ...)
    out <- .nplsqreg_record_omit(out, attr(mf, "na.action"))
    out$formula <- bws
    out <- .nplsqreg_set_response_name(
      out, .nplsqreg_formula_response_name(bws))
    out$call <- match.call(expand.dots = FALSE)
    environment(out$call) <- parent.frame()
    out
  }

nplsqregbw.lsqregressionbandwidth <- function(bws, tau = bws$tau, ...) {
  tau <- .nplsqreg_validate_tau_values(tau)
  if (length(tau) != length(bws$tau) || !isTRUE(all.equal(tau, bws$tau)))
    stop("cross-tau nplsqreg bandwidth-object reuse is not supported",
         call. = FALSE)
  bws
}

nplsqregbw.NULL <- function(...) {
  nplsqregbw.default(...)
}

nplsqregbw.default <-
  function(bws,
           xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           tau = 0.5,
           tau.search = c("full", "refined"),
           delta = NULL,
           scale = NULL,
           regtype = c("lc", "ll", "lp"),
           regtype.pilot = c("auto", "ll", "lc", "lp"),
           nomad = FALSE,
           nomad.pilot = FALSE,
           pilot.args = list(),
           bandwidth.compute = TRUE,
           delta.bounds = c(1e-4, 1 - 1e-4),
           optim.control = list(maxit = 50L),
           ...,
           nomad.opts = list()) {

    nomad.opts <- .np_nomad_normalize_user_opts(nomad.opts, "nplsqregbw")
    .npRmpi_require_active_slave_pool(where = "nplsqregbw()")
    mc.dispatch <- match.call()
    tau.dispatch <- try(.nplsqreg_validate_tau_values(tau),
                        silent = TRUE)
    if (.npRmpi_nplsqreg_should_autodispatch(mc.dispatch) &&
        !inherits(tau.dispatch, "try-error") &&
        length(tau.dispatch) == 1L)
      return(.npRmpi_autodispatch_call(mc.dispatch, environment()))

    elapsed.start <- proc.time()[3]
    progress.wrapped <- isTRUE(.np_progress_runtime$nplsqreg_bw_wrapped)
    if (isTRUE(bandwidth.compute) &&
        !progress.wrapped &&
        !.np_progress_bandwidth_active()) {
      mc <- match.call(expand.dots = TRUE)
      old.progress.wrapped <- .np_progress_runtime$nplsqreg_bw_wrapped
      .np_progress_runtime$nplsqreg_bw_wrapped <- TRUE
      on.exit({
        .np_progress_runtime$nplsqreg_bw_wrapped <- old.progress.wrapped
      }, add = TRUE)
      return(.np_progress_select_bandwidth_enhanced(
        "Selecting least-squares quantile regression bandwidth",
        eval.parent(mc)
      ))
    }

    tau.raw <- .nplsqreg_validate_tau_values(tau)
    tau.search <- .nplsqreg_validate_tau_search(tau.search)
    regtype.supplied <- !missing(regtype)
    if (regtype.supplied)
      regtype <- match.arg(regtype)
    if (length(tau.raw) > 1L) {
      bws.missing <- missing(bws)
      bw.list <- vector("list", length(tau.raw))
      shared.scale <- scale
      if (identical(tau.search, "refined")) {
        central <- which.min(abs(tau.raw - 0.5))
        fit.order <- c(
          central,
          setdiff(order(abs(tau.raw - tau.raw[[central]])), central)
        )
      } else {
        fit.order <- seq_along(tau.raw)
      }
      warm.start.from <- rep(NA_integer_, length(tau.raw))
      warm.start.degree <- vector("list", length(tau.raw))
      previous.idx <- NA_integer_
      full.extra.args <- .nplsqreg_normalize_dots(
        c(list(...), if (length(nomad.opts)) list(nomad.opts = nomad.opts) else list()),
        where = "nplsqregbw")
      refined.extra.args <- full.extra.args
      tau.search.controls <- NULL
      if (identical(tau.search, "refined")) {
        refined.extra.args$nmulti <- 1L
        refined.extra.args$powell.remin <- FALSE
        refined.extra.args$nomad.remin <- FALSE
        refined.extra.args$nomad.nmulti <- 0L
        tau.search.controls <- list(
          nmulti = 1L,
          powell.remin = FALSE,
          nomad.remin = FALSE,
          nomad.nmulti = 0L,
          anchor.controls = list(
            nmulti = if (is.null(full.extra.args$nmulti)) NA_integer_ else full.extra.args$nmulti,
            powell.remin = if (is.null(full.extra.args$powell.remin)) NA else full.extra.args$powell.remin,
            nomad.remin = if (is.null(full.extra.args$nomad.remin)) NA else full.extra.args$nomad.remin,
            nomad.nmulti = if (is.null(full.extra.args$nomad.nmulti)) NA_integer_ else full.extra.args$nomad.nmulti
          )
        )
      }
      for (j in fit.order) {
        extra.args <- if (identical(tau.search, "refined") && identical(j, central)) {
          full.extra.args
        } else {
          refined.extra.args
        }
        one.args <- list(
          xdat = xdat,
          ydat = ydat,
          tau = tau.raw[[j]],
          tau.search = "full",
          delta = delta,
          scale = shared.scale,
          regtype.pilot = regtype.pilot,
          nomad = nomad,
          nomad.pilot = nomad.pilot,
          pilot.args = pilot.args,
          bandwidth.compute = bandwidth.compute,
          delta.bounds = delta.bounds,
          optim.control = optim.control
        )
        if (regtype.supplied)
          one.args$regtype <- regtype
        if (!bws.missing) {
          one.args$bws <- bws
        } else if (identical(tau.search, "refined") &&
                   !is.na(previous.idx) &&
                   !is.null(bw.list[[previous.idx]]$reg.bws)) {
          warm.bws <- bw.list[[previous.idx]]$reg.bws
          one.args$bws <- as.numeric(warm.bws$bw)
          if (!is.null(warm.bws$degree) &&
              isTRUE(warm.bws$ncon > 0L) &&
              length(warm.bws$degree) == warm.bws$ncon) {
            warm.start.degree[[j]] <- as.integer(warm.bws$degree)
            extra.args$degree.start <- warm.bws$degree
          }
          warm.start.from[[j]] <- previous.idx
        }
        bw.list[[j]] <- do.call("nplsqregbw.default", c(one.args, extra.args))
        if (is.null(shared.scale))
          shared.scale <- bw.list[[j]]$scale
        previous.idx <- j
      }
      out <- .nplsqreg_combine_bandwidths(
        bw.list = bw.list,
        tau = tau.raw,
        tau.search = tau.search,
        fit.order = fit.order,
        warm.start.from = warm.start.from,
        warm.start.degree = warm.start.degree,
        tau.search.controls = tau.search.controls,
        call = match.call(expand.dots = FALSE))
      environment(out$call) <- parent.frame()
      return(out)
    }
    tau <- .nplsqreg_validate_tau(tau.raw)
    regtype.pilot <- match.arg(regtype.pilot)
    nomad <- npValidateScalarLogical(nomad, "nomad")
    nomad.pilot <- npValidateScalarLogical(nomad.pilot, "nomad.pilot")
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute,
                                                 "bandwidth.compute")
    dots <- .nplsqreg_normalize_dots(
      c(list(...), if (length(nomad.opts)) list(nomad.opts = nomad.opts) else list()),
      where = "nplsqregbw")
    controls <- .nplsqreg_optimizer_controls(dots, optim.control)
    prepared <- .nplsqreg_prepare_train_data(xdat, ydat, scale = scale)
    xdat <- prepared$xdat
    ydat <- prepared$ydat
    scale <- prepared$scale
    rows.omit <- prepared$omit
    if (!is.numeric(delta.bounds) || length(delta.bounds) != 2L ||
        any(!is.finite(delta.bounds)) || delta.bounds[1] <= 0 ||
        delta.bounds[2] >= 1 || delta.bounds[1] >= delta.bounds[2])
      stop("'delta.bounds' must be a two-value numeric interval inside (0, 1)",
           call. = FALSE)

    if (isTRUE(nomad) && !isTRUE(bandwidth.compute))
      stop("nplsqregbw nomad=TRUE requires bandwidth.compute=TRUE",
           call. = FALSE)
    if (isTRUE(nomad) && regtype.supplied && !identical(regtype, "lp"))
      stop("nplsqregbw nomad=TRUE requires regtype='lp' when regtype is supplied",
           call. = FALSE)
    reg.dots <- .nplsqreg_strip_optimizer_dots(dots)
    if (regtype.supplied)
      reg.dots$regtype <- regtype
    if (isTRUE(nomad)) {
      reg.dots$regtype <- NULL
      reg.dots$degree <- NULL
      reg.dots$basis <- NULL
      reg.dots$bernstein.basis <- NULL
    }
    opt.args <- utils::modifyList(dots, controls)

    if (is.null(scale)) {
      pilot <- .nplsqreg_scale_pilot(
        xdat = xdat,
        ydat = ydat,
        dots = reg.dots,
        regtype.pilot = regtype.pilot,
        nomad.pilot = nomad.pilot,
        pilot.args = pilot.args)
      scale <- pilot$scale
      scale.type <- if (isTRUE(pilot$nomad.pilot)) {
        "pilot-nomad"
      } else {
        paste0("pilot-", pilot$regtype.pilot)
      }
      mean.fit <- pilot$mean.fit
      scale.fit <- pilot$scale.fit
    } else {
      scale <- .nplsqreg_validate_scale(scale, length(ydat))
      scale.type <- "supplied"
      mean.fit <- NULL
      scale.fit <- NULL
    }

    if (missing(bws)) {
      start.bws <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat), reg.dots))
    } else if (isa(bws, "lsqregressionbandwidth")) {
      return(nplsqregbw.lsqregressionbandwidth(bws, tau = tau))
    } else if (isa(bws, "rbandwidth")) {
      if (identical(bws$method, "cv.check"))
        stop("fit$reg.bws is internal nplsqreg regression state; use fit$bws for exact nplsqreg bandwidth reuse or pass an ordinary npregbw rbandwidth object as a starting value",
             call. = FALSE)
      start.bws <- bws
    } else {
      start.bws <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat,
                                           bws = bws,
                                           bandwidth.compute = FALSE), reg.dots))
    }

    if (!(start.bws$type %in% c("fixed", "adaptive_nn", "generalized_nn")))
      stop("nplsqregbw requires bwtype='fixed', 'adaptive_nn', or 'generalized_nn'",
           call. = FALSE)

    if (is.null(delta))
      delta <- 0.5
    if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
        delta <= delta.bounds[1] || delta >= delta.bounds[2])
      stop("'delta' must be a single finite numeric value inside delta.bounds",
           call. = FALSE)

    search.result <- NULL
    if (isTRUE(nomad)) {
      degree.search <- .npregbw_degree_search_controls(
        regtype = "lp",
        regtype.named = TRUE,
        ncon = start.bws$ncon,
        nobs = nrow(xdat),
        basis = if (is.null(dots$basis)) start.bws$basis else dots$basis,
        degree.select = if (is.null(dots$degree.select)) "coordinate" else dots$degree.select,
        search.engine = if (is.null(dots$search.engine)) "nomad+powell" else dots$search.engine,
        degree.min = dots$degree.min,
        degree.max = dots$degree.max,
        degree.start = dots$degree.start,
        degree.restarts = if (is.null(dots$degree.restarts)) 0L else dots$degree.restarts,
        degree.max.cycles = if (is.null(dots$degree.max.cycles)) 25L else dots$degree.max.cycles,
        degree.verify = if (is.null(dots$degree.verify)) FALSE else dots$degree.verify,
        bernstein.basis = if (is.null(dots$bernstein.basis)) TRUE else dots$bernstein.basis,
        bernstein.named = !is.null(dots$bernstein.basis)
      )
      if (isTRUE(degree.search$singleton)) {
        fixed.reg.dots <- reg.dots
        fixed.reg.dots$regtype <- "lp"
        fixed.reg.dots$degree <- as.integer(degree.search$fixed.degree)
        fixed.reg.dots$bernstein.basis <- degree.search$bernstein.basis
        fixed.bws <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat),
                                        fixed.reg.dots))
        fixed.bws$regtype <- "lp"
        fixed.bws$degree <- as.integer(degree.search$fixed.degree)
        fixed.bws$bernstein.basis <- degree.search$bernstein.basis
        fixed.core <- .nplsqreg_call_fixed_degree_core(
          xdat = xdat,
          ydat = ydat,
          scale = scale,
          tau = tau,
          bws = fixed.bws,
          delta = delta,
          delta.bounds = delta.bounds,
          opt.args = opt.args,
          bandwidth.compute = bandwidth.compute
        )
        fixed.core$bws <- fixed.bws
        fixed.core$bws$bw <- fixed.core$bw
        fixed.core$bws$degree <- as.integer(degree.search$fixed.degree)
        fixed.core$bws$regtype <- "lp"
        fixed.core$bws$bernstein.basis <- degree.search$bernstein.basis
        search.result <- .np_degree_singleton_search_result(
          degree.search = degree.search,
          eval_result = list(
            objective = as.numeric(fixed.core$objective[1L]),
            payload = fixed.core,
            num.feval = if (!is.null(fixed.core$num.feval)) as.numeric(fixed.core$num.feval[1L]) else NA_real_
          ),
          direction = "min",
          objective_name = "fval"
        )
      } else {
        search.result <- .nplsqreg_nomad_search(
          xdat = xdat,
          ydat = ydat,
          scale = scale,
          tau = tau,
          template = start.bws,
          delta = delta,
          delta.bounds = delta.bounds,
          opt.args = opt.args,
          degree.search = degree.search,
          nomad.inner.nmulti = if (is.null(dots$nomad.nmulti)) 0L else dots$nomad.nmulti,
          random.seed = if (is.null(dots$random.seed)) 42L else dots$random.seed,
          nomad.opts = if (is.null(dots$nomad.opts)) list() else dots$nomad.opts
        )
      }
      core <- search.result$best_payload
      start.bws <- core$bws
    } else {
      core <- .nplsqreg_call_fixed_degree_core(
        xdat = xdat,
        ydat = ydat,
        scale = scale,
        tau = tau,
        bws = start.bws,
        delta = delta,
        delta.bounds = delta.bounds,
        opt.args = opt.args,
        bandwidth.compute = bandwidth.compute)
    }
    best.bw <- core$bw
    best.delta <- core$delta
    best.value <- core$objective
    qdat <- ydat + scale * stats::qnorm(best.delta)
    reg.bws <- .nplsqreg_rebuild_rbandwidth(
      template = start.bws,
      bw = best.bw,
      ydat = qdat,
      xdat = xdat,
      fval = best.value,
      num.feval = core$num.feval,
      bandwidth.compute = bandwidth.compute)
    reg.bws$method <- "cv.check"
    reg.bws$pmethod <- "Check-Loss Cross-Validation"
    reg.bws$ifval <- core$ifval
    reg.bws$num.feval.fast <- core$num.feval.fast
    reg.bws$nn.cache <- core$nn.cache
    reg.bws$fval.history <- core$fval.history
    reg.bws$eval.history <- core$eval.history
    reg.bws$invalid.history <- core$invalid.history
    reg.bws$timing <- core$timing
    reg.bws$total.time <- proc.time()[3] - elapsed.start
    reg.bws$rows.omit <- .nplsqreg_omit_rows(rows.omit)
    reg.bws$nobs.omit <- .nplsqreg_omit_length(rows.omit)
    if (!is.null(search.result)) {
      reg.bws <- .npregbw_attach_degree_search(reg.bws, search.result)
      reg.bws$nomad.shortcut <- list(enabled = TRUE, preset = "lp_nomad")
    }

    out <- lsqregressionbandwidth(
      reg.bws = reg.bws,
      tau = tau,
      delta = best.delta,
      scale = scale,
      scale.type = scale.type,
      xdat = xdat,
      ydat = ydat,
      qdat = qdat,
      objective = best.value,
      optim = core,
      mean.fit = mean.fit,
      scale.fit = scale.fit,
      formula = NULL,
      call = match.call(expand.dots = FALSE))
    out <- .nplsqreg_record_omit(out, rows.omit)
    out <- .nplsqreg_attach_native_diagnostics(
      out,
      .nplsqreg_native_diagnostics_from_search(search.result)
    )
    environment(out$call) <- parent.frame()
    out
  }

nplsqreg.formula <-
  function(bws, data = NULL, newdata = NULL, tau = 0.5,
           gradients = FALSE, residuals = FALSE, subset, na.action, ...) {

    .npRmpi_require_active_slave_pool(where = "nplsqreg()")

    tt <- terms(bws)
    dots <- list(...)
    response.name <- .nplsqreg_formula_response_name(bws)
    native.exdat <- dots$exdat
    dots$exdat <- NULL
    mc <- match.call(expand.dots = FALSE)
    m <- match(c("bws", "data", "subset", "na.action"),
               names(mc), nomatch = 0)
    tmf <- mc[c(1, m)]
    if ("bws" %in% names(tmf))
      names(tmf)[names(tmf) == "bws"] <- "formula"
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    mf <- do.call(stats::model.frame, mf.args, envir = environment(tt))
    train.omit <- attr(mf, "na.action")
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"), "term.labels"), drop = FALSE]

    has.eval <- !is.null(native.exdat) || !is.null(newdata)
    if (!is.null(native.exdat)) {
      npValidateNewdataFormula(native.exdat, delete.response(tt),
                               include.response = FALSE)
      emf <- do.call(stats::model.frame,
                     list(formula = delete.response(tt), data = native.exdat),
                     envir = parent.frame())
      eval.omit <- attr(emf, "na.action")
      exdat <- emf[, attr(attr(emf, "terms"), "term.labels"), drop = FALSE]
    } else if (has.eval) {
      npValidateNewdataFormula(newdata, delete.response(tt),
                               include.response = FALSE)
      emf <- do.call(stats::model.frame,
                     list(formula = delete.response(tt), data = newdata),
                     envir = parent.frame())
      eval.omit <- attr(emf, "na.action")
      exdat <- emf[, attr(attr(emf, "terms"), "term.labels"), drop = FALSE]
    } else {
      eval.omit <- NULL
    }

    out <- if (has.eval) {
      do.call(nplsqreg.default,
              c(list(txdat = xdat, tydat = ydat, tau = tau,
                     exdat = exdat, gradients = gradients,
                     residuals = residuals),
                dots))
    } else {
      do.call(nplsqreg.default,
              c(list(txdat = xdat, tydat = ydat, tau = tau,
                     gradients = gradients, residuals = residuals),
                dots))
    }
    out$call <- match.call(expand.dots = FALSE)
    environment(out$call) <- parent.frame()
    out$bws$formula <- bws
    out <- .nplsqreg_set_response_name(out, response.name)
    out <- .nplsqreg_record_omit(out, train.omit)
    out$bws <- .nplsqreg_record_omit(out$bws, train.omit)
    if (!has.eval)
      out <- .nplsqreg_pad_fit_outputs(out, train.omit)
    if (isTRUE(residuals) && .nplsqreg_omit_length(train.omit))
      out$resid <- .nplsqreg_napredict_eval(train.omit, out$resid)
    if (has.eval && .nplsqreg_omit_length(eval.omit)) {
      out$eval.omit <- eval.omit
      out$eval.rows.omit <- as.vector(eval.omit)
      out$eval.nobs.omit <- length(eval.omit)
      out <- .nplsqreg_pad_fit_outputs(out, eval.omit)
    }
    if (!is.null(out$bws$tau.bws))
      for (j in seq_along(out$bws$tau.bws))
        out$bws$tau.bws[[j]]$formula <- bws
    if (!is.null(out$tau.fits))
      for (j in seq_along(out$tau.fits))
        out$tau.fits[[j]]$bws$formula <- bws
    out
  }

nplsqreg.lsqregressionbandwidth <-
  function(bws, txdat = NULL, tydat = NULL, tau = bws$tau, ...) {
    tau <- .nplsqreg_validate_tau_values(tau)
    if (length(tau) != length(bws$tau) || !isTRUE(all.equal(tau, bws$tau)))
      stop("cross-tau nplsqreg bandwidth-object reuse is not supported",
           call. = FALSE)
    if (is.null(txdat))
      txdat <- bws$xdat
    if (is.null(tydat))
      tydat <- bws$ydat
    if (length(tau) > 1L) {
      if (is.null(bws$tau.bws) || length(bws$tau.bws) != length(tau))
        stop("vector nplsqreg bandwidth object lacks per-tau bandwidth state",
             call. = FALSE)
      fit.list <- lapply(seq_along(tau), function(j) {
        one.bws <- bws$tau.bws[[j]]
        one.bws$formula <- bws$formula
        nplsqreg.default(bws = one.bws,
                         txdat = txdat,
                         tydat = tydat,
                         tau = tau[[j]],
                         ...)
      })
      out <- .nplsqreg_combine_fits(
        fit.list = fit.list,
        tau = tau,
        bws = bws,
        tau.search = if (is.null(bws$tau.search)) "full" else bws$tau.search,
        call = match.call(expand.dots = FALSE))
      environment(out$call) <- parent.frame()
      return(out)
    }
    nplsqreg.default(bws = bws, txdat = txdat, tydat = tydat, tau = tau, ...)
  }

nplsqreg.NULL <- function(...) {
  nplsqreg.default(...)
}

nplsqreg.default <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           tau = 0.5,
           exdat,
           gradients = FALSE,
           residuals = FALSE,
           ...) {

    .npRmpi_require_active_slave_pool(where = "nplsqreg()")
    dots.dispatch <- list(...)
    if (length(dots.dispatch)) {
      dot.names <- names(dots.dispatch)
      if (!is.null(dot.names))
        # Make S3 `..n` placeholders resolvable by the MPI call materializer.
        for (nm in dot.names[!is.na(dot.names) & nzchar(dot.names)])
          if (!exists(nm, envir = environment(), inherits = FALSE))
            assign(nm, dots.dispatch[[nm]], envir = environment())
    }
    mc.dispatch <- match.call()
    tau.dispatch <- try(.nplsqreg_validate_tau_values(tau),
                        silent = TRUE)
    if (.npRmpi_nplsqreg_should_autodispatch(mc.dispatch) &&
        !inherits(tau.dispatch, "try-error") &&
        length(tau.dispatch) == 1L)
      return(.npRmpi_autodispatch_call(mc.dispatch, environment()))

    tau.raw <- .nplsqreg_validate_tau_values(tau)
    gradients <- npValidateScalarLogical(gradients, "gradients")
    residuals <- npValidateScalarLogical(residuals, "residuals")
    dots <- dots.dispatch
    native.newdata <- dots$newdata
    dots$newdata <- NULL
    dots$exdat <- NULL

    if (missing(bws) || !isa(bws, "lsqregressionbandwidth")) {
      bw.args <- list(xdat = txdat, ydat = tydat, tau = tau.raw)
      if (!missing(bws))
        bw.args$bws <- bws
      bw <- do.call("nplsqregbw", c(bw.args, dots))
      fit.args <- list(bws = bw, txdat = txdat, tydat = tydat,
                       gradients = gradients, residuals = residuals)
      if (!missing(exdat))
        fit.args$exdat <- exdat
      else if (!is.null(native.newdata))
        fit.args$exdat <- native.newdata
      return(do.call("nplsqreg.lsqregressionbandwidth", fit.args))
    }

    if (length(tau.raw) > 1L) {
      reuse.args <- list(bws = bws,
                         txdat = txdat,
                         tydat = tydat,
                         tau = tau.raw,
                         gradients = gradients,
                         residuals = residuals)
      if (!missing(exdat))
        reuse.args$exdat <- exdat
      else if (!is.null(native.newdata))
        reuse.args$newdata <- native.newdata
      return(do.call(nplsqreg.lsqregressionbandwidth,
                     c(reuse.args, dots)))
    }
    tau <- .nplsqreg_validate_tau(tau.raw)

    if (!isTRUE(all.equal(tau, bws$tau)))
      stop("cross-tau nplsqreg bandwidth-object reuse is not supported",
           call. = FALSE)

    prepared <- .nplsqreg_prepare_train_data(txdat, tydat)
    txdat <- prepared$xdat
    tydat <- prepared$ydat
    train.omit <- prepared$omit
    .nplsqreg_assert_reuse_training_data(bws, txdat, tydat)

    fit.args <- list(bws = bws$reg.bws, txdat = txdat, tydat = bws$qdat,
                     gradients = gradients)
    eval.present <- !missing(exdat) || !is.null(native.newdata)
    eval.omit <- NULL
    if (!missing(exdat)) {
      eval.prepared <- .nplsqreg_prepare_eval_data(exdat)
      fit.args$exdat <- eval.prepared$exdat
      eval.omit <- eval.prepared$omit
    } else if (!is.null(native.newdata)) {
      eval.prepared <- .nplsqreg_prepare_eval_data(native.newdata)
      fit.args$exdat <- eval.prepared$exdat
      eval.omit <- eval.prepared$omit
    }
    fit <- do.call(npreg, c(fit.args, dots))

    quant <- fitted(fit)
    qerr <- se(fit)
    qgrad <- if (gradients) gradients(fit) else NA
    qgerr <- if (gradients) gradients(fit, errors = TRUE) else NA
    resid.out <- if (residuals) {
      if (eval.present) {
        tydat - fitted(npreg(bws = bws$reg.bws,
                             txdat = txdat,
                             tydat = bws$qdat))
      } else {
        tydat - fitted(fit)
      }
    } else {
      NA
    }

    out <- lsqregression(
      bws = bws,
      fit = fit,
      xeval = fit$eval,
      tau = tau,
      delta = bws$delta,
      quantile = quant,
      quanterr = qerr,
      quantgrad = qgrad,
      quantgerr = qgerr,
      ntrain = nrow(txdat),
      trainiseval = !eval.present,
      gradients = gradients,
      residuals = residuals,
      resid = resid.out,
      call = match.call(expand.dots = FALSE))
    out <- .nplsqreg_attach_native_diagnostics(
      out,
      attr(bws, "native.nomad.diagnostics", exact = TRUE)
    )
    out <- .nplsqreg_record_omit(out, train.omit)
    if (!eval.present)
      out <- .nplsqreg_pad_fit_outputs(out, train.omit)
    if (isTRUE(residuals) && .nplsqreg_omit_length(train.omit))
      out$resid <- .nplsqreg_napredict_eval(train.omit, out$resid)
    if (eval.present && .nplsqreg_omit_length(eval.omit)) {
      out$eval.omit <- eval.omit
      out$eval.rows.omit <- as.vector(eval.omit)
      out$eval.nobs.omit <- length(eval.omit)
      out <- .nplsqreg_pad_fit_outputs(out, eval.omit)
    }
    environment(out$call) <- parent.frame()
    out
  }
