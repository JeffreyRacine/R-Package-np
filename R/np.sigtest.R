# This function implements an individual test of significance for both
# discrete (Racine, hart, Li, 2006, ER) and continuous variables
# (Racine, 1997, JBES). It accepts a data frame for explanatory data
# (mixed datatypes allowed), a vector for y for a regression model, an
# npregbw object, and a set of indices for the columns of X for which
# the test is to be run (default = all).

if (getRversion() >= "2.15.1")
  utils::globalVariables(".npsig_worker")

.npRmpi_npsig_extract_xy_from_bws <- function(obj) {
  if (!is.null(obj$call)) {
    call.names <- names(obj$call)
    if (!is.null(call.names) &&
        any(call.names == "xdat") &&
        any(call.names == "ydat")) {
      xdat <- .np_eval_bws_call_arg(obj, "xdat")
      ydat <- .np_eval_bws_call_arg(obj, "ydat")
      return(list(xdat = xdat, ydat = ydat))
    }
  }

  if (!is.null(obj$formula) && !is.null(obj$call)) {
    tt <- terms(obj$formula)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(obj$call), nomatch = 0)
    tmf <- obj$call[c(1, m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"), "term.labels"), drop = FALSE]
    return(list(xdat = xdat, ydat = ydat))
  }

  stop("unable to extract xdat/ydat from bandwidth object")
}

.npRmpi_npsig_validate_index <- function(index, xdat) {
  if(anyNA(index)) stop("index must not contain missing values")
  if(any(index < 1 | index > NCOL(xdat), na.rm = TRUE)) stop(paste("invalid index provided: index entries must lie between 1 and ",NCOL(xdat),sep=""))
  if(length(unique(index)) < length(index)) stop("index contains repeated values (must be unique)")
  invisible(TRUE)
}

.npRmpi_with_local_regression <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  old.local <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.autodispatch.context = TRUE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old.local), add = TRUE)
  old.mode <- .Call("C_np_set_local_regression_mode", TRUE, PACKAGE = "npRmpi")
  on.exit(.Call("C_np_set_local_regression_mode", old.mode, PACKAGE = "npRmpi"), add = TRUE)
  force(expr)
}

.npRmpi_npsig_npreg_local <- function(...) {
  .npRmpi_with_local_regression(npreg(...))
}

.npRmpi_npsig_do_local <- function(extra.args = NULL, ...) {
  args <- c(list(...), if (length(extra.args)) extra.args else NULL)
  do.call(.npRmpi_npsig_npreg_local, args)
}

.npRmpi_npsig_do_leaf <- function(fun, extra.args = NULL, ...) {
  args <- c(list(...), if (length(extra.args)) extra.args else NULL)
  .npRmpi_autodispatch_untag(do.call(fun, args))
}

.npRmpi_npsig_npreg_leaf <- function(extra.args = NULL, ...) {
  .npRmpi_npsig_do_leaf(npreg, extra.args = extra.args, ...)
}

.npRmpi_npsig_collective_context <- function() {
  isTRUE(.npRmpi_autodispatch_called_from_bcast())
}

.npRmpi_npsig_bootstrap_seed_plan <- function(num.obs,
                                              boot.num,
                                              boot.method,
                                              draw.wild.mult,
                                              a,
                                              b,
                                              p.a) {
  seeds <- vector("list", boot.num)
  for (i.star in seq_len(boot.num)) {
    seeds[[i.star]] <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (boot.method == "iid" || boot.method == "pairwise") {
      sample.int(num.obs, replace = TRUE)
    } else if (boot.method == "wild") {
      draw.wild.mult(num.obs, a, b, p.a)
    } else if (boot.method == "wild-rademacher") {
      draw.wild.mult(num.obs, -1, 1, p.a)
    } else {
      stop(sprintf("unsupported bootstrap method '%s'", boot.method), call. = FALSE)
    }
  }
  seeds
}

.npRmpi_npsig_gather_rank_chunks <- function(gathered, size) {
  size <- as.integer(size)[1L]
  if (is.na(size) || size < 1L)
    stop("invalid MPI gather size")

  if (is.matrix(gathered)) {
    if (!identical(ncol(gathered), size))
      stop("npsigtest MPI gather returned malformed matrix output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.array(gathered)) {
    dims <- dim(gathered)
    if (length(dims) < 2L || !identical(dims[[length(dims)]], size))
      stop("npsigtest MPI gather returned malformed array output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.list(gathered)) {
    if (!identical(length(gathered), size))
      stop("npsigtest MPI gather returned malformed list output", call. = FALSE)
    return(gathered)
  }

  chunks <- as.list(gathered)
  if (!identical(length(chunks), size))
    stop("npsigtest MPI gather returned malformed atomic output", call. = FALSE)
  chunks
}

.npRmpi_npsig_parallel_boot_values_collective <- function(boot.seeds,
                                                          worker,
                                                          comm = 1L) {
  n.boot <- length(boot.seeds)
  if (n.boot < 1L)
    return(numeric(0))

  size <- mpi.comm.size(comm)
  if (size < 2L) {
    local.idx <- seq_len(n.boot)
    return(as.numeric(worker(local.idx, boot.seeds)))
  }

  rank <- mpi.comm.rank(comm)
  local.idx <- seq.int(rank + 1L, n.boot, by = size)
  local.vals <- if (length(local.idx)) {
    as.numeric(worker(local.idx, boot.seeds))
  } else {
    numeric(0)
  }

  gathered <- mpi.gather.Robj(local.vals, root = 0L, comm = comm)
  if (rank == 0L) {
    out <- numeric(n.boot)
    gathered <- .npRmpi_npsig_gather_rank_chunks(gathered = gathered, size = size)
    for (r in seq_len(size)) {
      idx.r <- seq.int(r, n.boot, by = size)
      vals.r <- as.numeric(gathered[[r]])
      if (length(idx.r) != length(vals.r))
        stop("npsigtest MPI gather returned mismatched bootstrap chunk lengths", call. = FALSE)
      if (length(idx.r))
        out[idx.r] <- vals.r
    }
    mpi.bcast.Robj(out, rank = 0L, comm = comm)
    out
  } else {
    mpi.bcast.Robj(rank = 0L, comm = comm)
  }
}

.npRmpi_npsig_bootstrap_tasks <- function(n.boot, chunk.size) {
  starts <- seq.int(1L, n.boot, by = chunk.size)
  lapply(starts, function(start) {
    list(
      start = as.integer(start),
      bsz = as.integer(min(chunk.size, n.boot - start + 1L))
    )
  })
}

.npRmpi_npsig_parallel_boot_values <- function(boot.seeds,
                                               worker,
                                               required.bindings = NULL,
                                               what = "npsigtest",
                                               profile.where = NA_character_,
                                               progress.context = NULL,
                                               comm = 1L) {
  n.boot <- length(boot.seeds)
  if (n.boot < 1L)
    return(numeric(0))

  if (.npRmpi_npsig_collective_context()) {
    return(.npRmpi_npsig_parallel_boot_values_collective(
      boot.seeds = boot.seeds,
      worker = worker,
      comm = comm
    ))
  }

  if (!isTRUE(.npRmpi_has_active_slave_pool(comm = comm))) {
    return(as.numeric(worker(seq_len(n.boot), boot.seeds)))
  }

  workers <- max(1L, .npRmpi_bootstrap_worker_count(comm = comm))
  chunk.size <- max(1L, as.integer(floor(n.boot / workers)))
  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = n.boot,
    chunk.size = chunk.size,
    comm = comm,
    include.master = TRUE
  )
  tasks <- .npRmpi_npsig_bootstrap_tasks(n.boot = n.boot, chunk.size = chunk.size)

  bindings <- c(list(.npsig_worker = worker), required.bindings)

  worker.chunk <- function(task, boot.seeds) {
    idx <- seq.int(task$start, length.out = task$bsz)
    vals <- as.numeric(.npsig_worker(idx, boot.seeds))
    matrix(vals, nrow = task$bsz, ncol = 1L)
  }

  out <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker.chunk,
    ncol.out = 1L,
    what = what,
    profile.where = profile.where,
    comm = comm,
    master_local_chunk = TRUE,
    required.bindings = bindings,
    progress.context = progress.context,
    boot.seeds = boot.seeds
  )

  as.numeric(out[, 1L])
}

.npRmpi_npsig_bootstrap_bw_reselect <- function(xdat,
                                                ydat,
                                                bws.seed,
                                                extra.args = list(),
                                                bootstrap.iter,
                                                bw.fun = npregbw,
                                                localize = TRUE) {
  bw.args <- if (length(extra.args)) extra.args else list()
  bw.args[c("xdat", "ydat", "bws")] <- NULL

  user.nmulti <- !is.null(names(bw.args)) &&
    "nmulti" %in% names(bw.args) &&
    !is.null(bw.args$nmulti)

  if (!user.nmulti && bootstrap.iter > 1L)
    bw.args$nmulti <- 1L

  call.args <- c(list(xdat = xdat, ydat = ydat, bws = bws.seed), bw.args)

  result <- .np_progress_with_legacy_suppressed(
    if (localize) {
      .npRmpi_with_local_regression(do.call(bw.fun, call.args))
    } else {
      do.call(bw.fun, call.args)
    }
  )

  .npRmpi_autodispatch_untag(result)
}

npsigtest <-
  function(bws, ..., B = 399){
    args <- list(...)
    npRejectLegacyBootstrapCount(names(args), "npsigtest")

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$xdat))
          UseMethod("npsigtest",bws$formula)
#        else if (!is.null(bws$call) && is.null(args$xdat) && (class(bws) != "npregression"))
        else if (!is.null(bws$call) && is.null(args$xdat) && (!isa(bws,"npregression")))
          UseMethod("npsigtest",bws$call)
        else if (!is.call(bws))
          UseMethod("npsigtest",bws)
        else
          UseMethod("npsigtest",NULL)
      } else {
        UseMethod("npsigtest", NULL)
      }
    } else {
      UseMethod("npsigtest", NULL)
    }
  }

npsigtest.formula <-
  function(bws, data = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    if (!is.null(data))
      tmf[["data"]] <- substitute(data)
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]

    ev <- npsigtest(xdat = xdat, ydat = ydat, bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    return(ev)
  }

npsigtest.call <-
  function(bws, ...) {
    ev <- npsigtest(xdat = .np_eval_bws_call_arg(bws, "xdat"),
                    ydat = .np_eval_bws_call_arg(bws, "ydat"),
                    bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npsigtest.npregression <-
  function(bws, ...){
    ev <- npsigtest(bws$bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

.np_npsig_pivot_plan <- function(pivot, xdat, index, joint) {
  categorical <- vapply(
    xdat[index],
    function(variable) is.factor(variable) || is.ordered(variable),
    logical(1L)
  )

  if (!is.null(pivot)) {
    pivot <- npValidateScalarLogical(pivot, "pivot")
    if (pivot && any(categorical)) {
      stop(sprintf(
        paste0(
          "pivot = TRUE is supported only when every tested predictor is ",
          "continuous; the published categorical test is unstandardized ",
          "(categorical predictor%s: %s). Use pivot = NULL (the ",
          "method-aware default) or pivot = FALSE."
        ),
        if (sum(categorical) == 1L) "" else "s",
        paste(names(xdat)[index][categorical], collapse = ", ")
      ), call. = FALSE)
    }
    effective <- rep.int(pivot, length(index))
  } else if (joint == TRUE) {
    effective <- rep.int(!any(categorical), length(index))
  } else {
    effective <- !categorical
  }

  names(effective) <- names(xdat)[index]
  list(requested = pivot, effective = effective)
}

.np_npsig_statistic <- function(fit, index, pivot) {
  gradient <- fit$grad[, index, drop = FALSE]
  if (any(!is.finite(gradient)))
    stop("npsigtest cannot construct a statistic from non-finite gradient estimates", call. = FALSE)

  if (!pivot)
    return(mean(gradient^2))

  if (is.null(fit$gerr))
    stop("npsigtest cannot construct a pivotal statistic because gradient standard errors are unavailable", call. = FALSE)

  gradient.stderr <- fit$gerr[, index, drop = FALSE]
  if (any(!is.finite(gradient.stderr)) || any(gradient.stderr <= 0.0)) {
    stop(paste0(
      "npsigtest cannot construct a pivotal statistic because a required ",
      "gradient standard error is non-positive or non-finite"
    ), call. = FALSE)
  }

  statistic <- mean((gradient / gradient.stderr)^2)
  if (!is.finite(statistic))
    stop("npsigtest pivotal statistic is non-finite", call. = FALSE)
  statistic
}

.np_npsig_progress_promote <- function(state, total, done) {
  if (!isTRUE(state$known_total)) {
    state$known_total <- TRUE
    state$total <- total
    state$throttle_sec <- .np_progress_interval_sec(
      known_total = TRUE,
      domain = state$domain
    )
  }
  .np_progress_step_at(
    state,
    now = .np_progress_now(),
    done = done,
    force = TRUE
  )
}

.np_npsig_upper_tail_p <- function(bootstrap, observed) {
  mean(bootstrap >= observed)
}

.np_npsig_streamed_iid_eligible <- function(bws,
                                             xdat,
                                             index,
                                             joint,
                                             boot.type,
                                             boot.method,
                                             pivot,
                                             extra.args) {
  tested <- xdat[index]
  categorical <- vapply(
    tested,
    function(variable) is.factor(variable) || is.ordered(variable),
    logical(1L)
  )
  equivalent.pivot <- if (joint) {
    is.null(pivot) || identical(pivot, FALSE) ||
      (identical(pivot, TRUE) && !any(categorical))
  } else {
    is.null(pivot) ||
      (identical(pivot, FALSE) && all(categorical)) ||
      (identical(pivot, TRUE) && !any(categorical))
  }

  regression.engine <- bws[["regtype.engine", exact = TRUE]]
  engine.supported <- regression.engine %in% c("lc", "lp")

  if (!identical(boot.type, "I") ||
      !boot.method %in% c("iid", "wild", "wild-rademacher") ||
      !equivalent.pivot ||
      length(extra.args) ||
      !bws[["type", exact = TRUE]] %in%
        c("fixed", "generalized_nn", "adaptive_nn") ||
      !engine.supported ||
      identical(bws[["ckertype", exact = TRUE]], "beta"))
    return(FALSE)

  degree <- bws[["degree.engine", exact = TRUE]]
  if (bws[["ncon", exact = TRUE]] > 0L) {
    if (!is.numeric(degree) ||
        length(degree) != bws[["ncon", exact = TRUE]] ||
        anyNA(degree) || any(!is.finite(degree)) ||
        any(degree < 0) || any(degree != floor(degree)) ||
        (identical(regression.engine, "lc") && any(degree != 0L)))
      return(FALSE)
  }

  all(vapply(
    tested,
    function(variable)
      is.factor(variable) || is.ordered(variable) ||
        inherits(variable, c("integer", "numeric")),
    logical(1L)
  ))
}

.np_npsig_streamed_iid_tile <- function(bws,
                                         xdat,
                                         tested.index,
                                         donor.index = NULL,
                                         response.matrix = NULL,
                                         null.mean,
                                         residual.pool,
                                         pivotal = NULL) {
  continuous <- which(bws[["icon", exact = TRUE]])
  unordered <- which(bws[["iuno", exact = TRUE]])
  ordered <- which(bws[["iord", exact = TRUE]])
  gradient.order <- integer(bws[["ncon", exact = TRUE]])

  if (tested.index %in% continuous) {
    mode <- 1L
    coordinate <- match(tested.index, continuous) - 1L
    gradient.order[coordinate + 1L] <- 1L
  } else if (tested.index %in% unordered) {
    mode <- 2L
    coordinate <- match(tested.index, unordered) - 1L
  } else if (tested.index %in% ordered) {
    mode <- 3L
    coordinate <- match(tested.index, ordered) - 1L
  } else {
    stop("private npsigtest tile received an unsupported predictor", call. = FALSE)
  }

  if (is.null(pivotal))
    pivotal <- identical(mode, 1L)
  pivotal <- npValidateScalarLogical(pivotal, "pivotal")
  if (pivotal && !identical(mode, 1L))
    stop("private npsigtest categorical tiles cannot be pivotal", call. = FALSE)

  response.ready <- !is.null(response.matrix)
  if (response.ready == !is.null(donor.index))
    stop("private npsigtest tile requires exactly one response payload", call. = FALSE)

  as.numeric(.npreghat_exact_lp_apply_from_regression_core(
    bws = bws,
    txdat = xdat,
    y = if (response.ready) response.matrix else donor.index,
    basis = bws[["basis.engine", exact = TRUE]],
    degree = as.integer(bws[["degree.engine", exact = TRUE]]),
    bernstein.basis = bws[["bernstein.basis.engine", exact = TRUE]],
    s = gradient.order,
    sigtest = list(
      mode = mode,
      coordinate = coordinate,
      response.ready = response.ready,
      pivotal = pivotal,
      null.mean = null.mean,
      residual.pool = residual.pool
    )
  ))
}

.np_npsig_streamed_response_statistic <- function(bws,
                                                    xdat,
                                                    index,
                                                    response.matrix,
                                                    pivotal) {
  statistic <- numeric(ncol(response.matrix))
  placeholder <- response.matrix[, 1L]
  for (tested.index in index) {
    statistic <- statistic + .np_npsig_streamed_iid_tile(
      bws = bws,
      xdat = xdat,
      tested.index = tested.index,
      response.matrix = response.matrix,
      null.mean = placeholder,
      residual.pool = placeholder,
      pivotal = pivotal
    ) / length(index)
  }
  statistic
}

.np_npsig_validate_lp_degree <- function(bws, xdat, index) {
  if (!identical(bws[["regtype", exact = TRUE]], "lp"))
    return(invisible(TRUE))

  icon <- bws[["icon", exact = TRUE]]
  degree <- bws[["degree.engine", exact = TRUE]]
  valid.degree <- is.numeric(degree) &&
    all(is.finite(degree)) &&
    all(degree >= 0) &&
    all(degree == floor(degree))
  if (!is.logical(icon) || length(icon) != ncol(xdat) || anyNA(icon) ||
      length(degree) != sum(icon) || !valid.degree) {
    stop("npsigtest() encountered inconsistent local-polynomial degree metadata", call. = FALSE)
  }

  continuous <- which(icon)
  tested.continuous <- index[index %in% continuous]
  position <- match(tested.continuous, continuous)
  zero.degree <- tested.continuous[degree[position] == 0]
  if (length(zero.degree)) {
    stop(sprintf(
      paste0(
        "npsigtest() requires polynomial degree at least 1 for every tested ",
        "continuous coordinate when regtype = 'lp'; degree 0 was selected ",
        "for: %s. If this model was selected via NOMAD, refit with ",
        "degree.min = 1."
      ),
      paste(names(xdat)[zero.degree], collapse = ", ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}

npsigtest.rbandwidth <- function(bws,
                                 xdat = stop("data xdat missing"),
                                 ydat = stop("data ydat missing"),
                                 B = 399,
                                 boot.method = c("iid","wild","wild-rademacher","pairwise"),
                                 boot.type = c("I","II"),
                                 pivot = NULL,
                                 joint = FALSE,
                                 index = seq_len(ncol(xdat)),
                                 random.seed = 42,
                                 ...) {
  .npRmpi_require_active_slave_pool(where = "npsigtest()")

  xdat <- toFrame(xdat)

  if(B < 9) stop("number of bootstrap replications must be >= 9")

  ## catch and destroy NA's
  goodrows <- seq_len(nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat,ydat)), "na.action")
  goodrows[rows.omit] <- 0

  if (all(goodrows==0))
    stop("Data has no rows without NAs")

  xdat <- xdat[goodrows,,drop = FALSE]
  ydat <- ydat[goodrows]

  ## Fast-fail contract for invalid index values must run before any
  ## expensive/distributed bandwidth/regression work.
  .npRmpi_npsig_validate_index(index = index, xdat = xdat)

  pivot.plan <- .np_npsig_pivot_plan(
    pivot = pivot,
    xdat = xdat,
    index = index,
    joint = joint
  )
  pivot <- pivot.plan$requested

  bws <- local({
    old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
    options(npRmpi.autodispatch.disable = TRUE)
    on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
    npregbw(xdat = xdat, ydat = ydat, bws = bws, bandwidth.compute = FALSE)
  })
  bws <- .npRmpi_autodispatch_untag(bws)

  if (is.factor(ydat))
    stop("dependent variable must be continuous.")

  .np_npsig_validate_lp_degree(bws = bws, xdat = xdat, index = index)

  extra.args <- list(...)
  npRejectLegacyBootstrapCount(names(extra.args), "npsigtest")
  boot.type <- match.arg(boot.type)
  boot.method <- match.arg(boot.method)
  direct.statistic <- .np_npsig_streamed_iid_eligible(
    bws = bws,
    xdat = xdat,
    index = index,
    joint = joint,
    boot.type = "I",
    boot.method = "iid",
    pivot = pivot,
    extra.args = extra.args
  )
  streamed.iid <- direct.statistic && identical(boot.type, "I") &&
    boot.method %in% c("iid", "wild", "wild-rademacher")
  direct.pairwise <- direct.statistic && !joint && identical(boot.type, "I") &&
    identical(boot.method, "pairwise")

  tested.names <- names(xdat)[index]
  missing.names <- is.na(tested.names) | !nzchar(tested.names)
  tested.names[missing.names] <- paste("variable", index[missing.names])
  progress <- .np_progress_begin(
    if (joint) "Testing joint significance" else paste("Testing", tested.names[[1L]]),
    surface = "bootstrap"
  )
  progress <- .np_progress_show_now(progress)
  progress.active <- TRUE
  progress.context <- new.env(parent = emptyenv())
  on.exit({
    if (isTRUE(progress.active))
      .np_progress_abort(progress)
  }, add = TRUE)

  ## Save seed prior to setting.  The computational owner is selected before
  ## this boundary so unsupported calls enter the incumbent without RNG work.

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  collective.mode <- .npRmpi_npsig_collective_context()
  if (boot.type == "II")
    bws.original <- bws

  num.obs <- nrow(xdat)
  npreg.eval.fun <- if (boot.type == "II") .npRmpi_npsig_npreg_leaf else .npRmpi_npsig_do_local

  if(!joint) {

    In <- numeric(length(index))
    P <- numeric(length(index))

  }

  ## Some constants for the wild bootstrap

  a <- -0.6180339887499  # (1-sqrt(5))/2
  b <- 1.6180339887499   # (1+sqrt(5))/2
  P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))
  P.rademacher <- 0.5

  draw.wild.mult <- function(n.obs, a, b, p.a) {
    u <- stats::runif(n.obs)
    mult <- rep.int(b, n.obs)
    mult[u <= p.a] <- a
    mult
  }

  ## A vector for storing the resampled statistics

  In.vec <- numeric(B)

  if(joint==TRUE) {

    pivot.use <- pivot.plan$effective[[1L]]

    ## Joint test

    In.mat = matrix(data = 0, ncol = 1, nrow = B)

    if (boot.type == "II")
      bws <- bws.original

    ## Note - xdat must be a data frame

    ## Construct In, the average value of the squared derivatives of
    ## the jth element, discrete or continuous

    progress <- .np_progress_step(progress)
    npreg.out <- npreg.eval.fun(extra.args,
                                txdat = xdat,
                                tydat = ydat,
                                bws = bws,
                                gradients = TRUE,
                                se = pivot.use)

    In <- .np_npsig_statistic(npreg.out, index = index, pivot = pivot.use)
    progress <- .np_progress_step(progress)

    if(boot.method != "pairwise") {

      ## Compute scale and mean of unrestricted residuals

      progress <- .np_progress_step(progress)
      npreg.unres <- npreg.eval.fun(extra.args,
                                    txdat = xdat,
                                    tydat = ydat,
                                    bws = bws,
                                    residuals = TRUE)
      ei.unres <- scale(npreg.unres$resid)
      ei.unres.scale <- attr(ei.unres,"scaled:scale")
      ei.unres.center <- attr(ei.unres,"scaled:center")      
      progress <- .np_progress_step(progress)

      ## We now construct mhat.xi holding constant the variable whose
      ## significance is being tested at its median. First, make a copy
      ## of the data frame xdat
      
      xdat.eval <- xdat
      
      ## Impose the null by evaluating the conditional mean holding
      ## xdat[,i] constant at its median (numeric) or mode
      ## (factor/ordered) using uocquantile()

      for(i in index) {
        xq <- uocquantile(xdat[,i], 0.5)
        if (is.factor(xdat[,i]) || is.ordered(xdat[,i])) {
          xdat.eval[,i] <- cast(xq, xdat[,i], same.levels = TRUE)
        } else {
          xdat.eval[,i] <- xq
        }
      }
      
      progress <- .np_progress_step(progress)
      mhat.xi <-  npreg.eval.fun(extra.args,
                                 txdat = xdat,
                                 tydat = ydat,
                                 exdat = xdat.eval,
                                 bws = bws)$mean
      progress <- .np_progress_step(progress)

      ## Rescale and recenter the residuals under the null to those
      ## under the alternative
      
      ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
      
      ## Recenter the residuals to have mean zero

      ei <- ei - mean(ei)
      
    }

    if (boot.type == "II") {
      bws.boot.prev <- bws.original

      for (i.star in seq_len(B)) {
        if (boot.method == "iid") {
          ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
        } else if (boot.method == "wild") {
          ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
        } else if (boot.method == "wild-rademacher") {
          ydat.star <- mhat.xi + ei * draw.wild.mult(
            num.obs, -1, 1, P.rademacher
          )
        } else {
          boot.index <- sample.int(num.obs, replace = TRUE)
          ydat.star <- ydat[boot.index]
          xdat.star <- xdat[boot.index,]
          for (jj in index)
            xdat.star[, jj] <- xdat[, jj]
        }

        if (boot.method == "pairwise") {
          bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
            xdat = xdat.star,
            ydat = ydat.star,
            bws.seed = bws.boot.prev,
            extra.args = extra.args,
            bootstrap.iter = i.star,
            localize = FALSE
          )
        } else {
          bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
            xdat = xdat,
            ydat = ydat.star,
            bws.seed = bws.boot.prev,
            extra.args = extra.args,
            bootstrap.iter = i.star,
            localize = FALSE
          )
        }

        bws.boot.prev <- bws.boot
        bws <- bws.original
        bws$bw[index] <- bws.boot$bw[index]

        if (boot.method == "pairwise") {
          npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                 txdat = xdat.star,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE,
                                                 se = pivot.use)
        } else {
          npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE,
                                                 se = pivot.use)
        }

        In.vec[i.star] <- .np_npsig_statistic(
          npreg.boot,
          index = index,
          pivot = pivot.use
        )
        if (!isTRUE(progress$known_total)) {
          progress <- .np_npsig_progress_promote(
            progress, total = B, done = i.star
          )
        } else {
          progress <- .np_progress_step(progress, done = i.star)
        }
      }
    } else {
      boot.seeds <- .npRmpi_npsig_bootstrap_seed_plan(
        num.obs = num.obs,
        boot.num = B,
        boot.method = boot.method,
        draw.wild.mult = draw.wild.mult,
        a = a,
        b = b,
        p.a = if (identical(boot.method, "wild-rademacher")) P.rademacher else P.a
      )
      post.boot.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

      joint.eval <- function(task.idx, seed.plan) {
        out <- numeric(length(task.idx))
        if (streamed.iid) {
          tile.width <- 8L
          for (tile.start in seq.int(1L, length(task.idx), by = tile.width)) {
            tile.count <- min(tile.width, length(task.idx) - tile.start + 1L)
            tile.position <- tile.start:(tile.start + tile.count - 1L)
            response.matrix <- vapply(
              tile.position,
              function(position) {
                assign(".Random.seed", seed.plan[[task.idx[position]]],
                       envir = .GlobalEnv)
                if (identical(boot.method, "iid")) {
                  donor <- sample.int(num.obs, replace = TRUE)
                  mhat.xi + ei[donor]
                } else {
                  wild.values <- if (identical(boot.method, "wild"))
                    c(a, b, P.a) else c(-1, 1, P.rademacher)
                  mhat.xi + ei * draw.wild.mult(
                    num.obs, wild.values[[1L]], wild.values[[2L]],
                    wild.values[[3L]]
                  )
                }
              },
              numeric(num.obs)
            )
            out[tile.position] <- .np_npsig_streamed_response_statistic(
              bws = bws,
              xdat = xdat,
              index = index,
              response.matrix = response.matrix,
              pivotal = pivot.use
            )
          }
          return(out)
        }
        for (kk in seq_along(task.idx)) {
          assign(".Random.seed", seed.plan[[task.idx[kk]]], envir = .GlobalEnv)
          if (boot.method == "iid") {
            ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE,
                                                 se = pivot.use)
          } else if (boot.method == "wild") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE,
                                                 se = pivot.use)
          } else if (boot.method == "wild-rademacher") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(
              num.obs, -1, 1, P.rademacher
            )
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE,
                                                 se = pivot.use)
          } else {
            boot.index <- sample.int(num.obs, replace = TRUE)
            ydat.star <- ydat[boot.index]
            xdat.star <- xdat[boot.index,]
            for (jj in index)
              xdat.star[, jj] <- xdat[, jj]
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat.star,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE,
                                                 se = pivot.use)
          }

          out[kk] <- .np_npsig_statistic(
            npreg.boot,
            index = index,
            pivot = pivot.use
          )
        }
        out
      }

      joint.bindings <- list(
        boot.method = boot.method,
        xdat = xdat,
        ydat = ydat,
        bws = bws,
        index = index,
        pivot.use = pivot.use,
        .np_npsig_statistic = .np_npsig_statistic,
        num.obs = num.obs,
        draw.wild.mult = draw.wild.mult,
        a = a,
        b = b,
        P.a = P.a,
        extra.args = extra.args,
        streamed.iid = streamed.iid,
        .np_npsig_streamed_iid_tile = .np_npsig_streamed_iid_tile,
        .np_npsig_streamed_response_statistic =
          .np_npsig_streamed_response_statistic
      )
      if (boot.method != "pairwise")
        joint.bindings <- c(joint.bindings, list(mhat.xi = mhat.xi, ei = ei))

      progress$known_total <- TRUE
      progress$total <- B
      progress$throttle_sec <- .np_progress_interval_sec(
        known_total = TRUE,
        domain = progress$domain
      )
      progress.context$state <- progress
      progress.context$done <- NULL
      progress.context$use.bootstrap.done <- TRUE
      progress.context$force.next <- TRUE

      In.vec <- .npRmpi_npsig_parallel_boot_values(
        boot.seeds = boot.seeds,
        worker = joint.eval,
        required.bindings = joint.bindings,
        what = "npsigtest",
        profile.where = "npsigtest:joint",
        progress.context = progress.context
      )
      assign(".Random.seed", post.boot.seed, envir = .GlobalEnv)
      progress <- progress.context$state
      progress <- .np_npsig_progress_promote(
        progress, total = B, done = B
      )
    }

    ## Compute the P-value

    P <- .np_npsig_upper_tail_p(In.vec, In)

    In.mat[,1] = In.vec

  } else {

    ## Individual test

    ## ii is the counter for successive elements of In and P...

    In.mat = matrix(data = 0, ncol = length(index), nrow = B)

    ii <- 0

    if (streamed.iid) {
      progress <- .np_progress_step(progress)
      streamed.unrestricted <- npreg.eval.fun(
        extra.args,
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = TRUE,
        se = any(pivot.plan$effective)
      )
      progress <- .np_progress_step(progress)
      progress <- .np_progress_step(progress)
      streamed.unres <- npreg.eval.fun(
        extra.args,
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        residuals = TRUE
      )
      streamed.ei.unres <- scale(streamed.unres$resid)
      streamed.ei.unres.scale <- attr(streamed.ei.unres, "scaled:scale")
      streamed.ei.unres.center <- attr(streamed.ei.unres, "scaled:center")
      streamed.unres <- NULL
      streamed.ei.unres <- NULL
      progress <- .np_progress_step(progress)
    }

    for(i in index) {
      
      ## Increment counter...
      
      ii <- ii + 1
      progress$label <- paste("Testing", tested.names[[ii]])
      pivot.use <- pivot.plan$effective[[ii]]

      if (boot.type == "II")
        bws <- bws.original
      
      ## Note - xdat must be a data frame
      
      ## Construct In, the average value of the squared derivatives of
      ## the jth element, discrete or continuous
      
      if (streamed.iid) {
        npreg.out <- streamed.unrestricted
      } else {
        progress <- .np_progress_step(progress)
        npreg.out <- npreg.eval.fun(extra.args,
                                    txdat = xdat,
                                    tydat = ydat,
                                    bws = bws,
                                    gradients = TRUE,
                                    se = pivot.use)
      }
      
      In[ii] <- .np_npsig_statistic(npreg.out, index = i, pivot = pivot.use)
      progress <- .np_progress_step(progress)
      
      if(boot.method != "pairwise") {

        ## Compute scale and mean of unrestricted residuals

        if (streamed.iid) {
          ei.unres.scale <- streamed.ei.unres.scale
          ei.unres.center <- streamed.ei.unres.center
        } else {
          progress <- .np_progress_step(progress)
          npreg.unres <- npreg.eval.fun(extra.args,
                                        txdat = xdat,
                                        tydat = ydat,
                                        bws = bws,
                                        residuals = TRUE)
          ei.unres <- scale(npreg.unres$resid)
          ei.unres.scale <- attr(ei.unres,"scaled:scale")
          ei.unres.center <- attr(ei.unres,"scaled:center")
          progress <- .np_progress_step(progress)
        }

        ## We now construct mhat.xi holding constant the variable whose
        ## significance is being tested at its median. First, make a copy
        ## of the data frame xdat
        
        xdat.eval <- xdat
        
        ## Impose the null by evaluating the conditional mean holding
        ## xdat[,i] constant at its median (numeric) or mode
        ## (factor/ordered) using uocquantile()
        
        xq <- uocquantile(xdat[,i], 0.5)
        if (is.factor(xdat[,i]) || is.ordered(xdat[,i])) {
          xdat.eval[,i] <- cast(xq, xdat[,i], same.levels = TRUE)
        } else {
          xdat.eval[,i] <- xq
        }
        
        progress <- .np_progress_step(progress)
        mhat.xi <-  npreg.eval.fun(extra.args,
                                   txdat = xdat,
                                   tydat = ydat,
                                   exdat = xdat.eval,
                                   bws = bws)$mean
        progress <- .np_progress_step(progress)
        
        ## Rescale and recenter the residuals under the null to those
        ## under the alternative
        
        ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
        
        ## Recenter the residuals to have mean zero
        
        ei <- ei - mean(ei)
        
      }

      if (boot.type == "II") {
        bws.boot.prev <- bws.original

        for (i.star in seq_len(B)) {
          if (boot.method == "iid") {
            ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
          } else if (boot.method == "wild") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
          } else if (boot.method == "wild-rademacher") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(
              num.obs, -1, 1, P.rademacher
            )
          } else {
            boot.index <- sample.int(num.obs, replace = TRUE)
            ydat.star <- ydat[boot.index]
            xdat.star <- xdat
            xdat.star[, -i] <- xdat[boot.index, -i]
          }

          if (boot.method == "pairwise") {
            bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
              xdat = xdat.star,
              ydat = ydat.star,
              bws.seed = bws.boot.prev,
              extra.args = extra.args,
              bootstrap.iter = i.star,
              localize = FALSE
            )
          } else {
            bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
              xdat = xdat,
              ydat = ydat.star,
              bws.seed = bws.boot.prev,
              extra.args = extra.args,
              bootstrap.iter = i.star,
              localize = FALSE
            )
          }

          bws.boot.prev <- bws.boot
          bws <- bws.original
          bws$bw[i] <- bws.boot$bw[i]

          if (boot.method == "pairwise") {
            npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                   txdat = xdat.star,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE,
                                                   se = pivot.use)
          } else {
            npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE,
                                                   se = pivot.use)
          }

          In.vec[i.star] <- .np_npsig_statistic(
            npreg.boot,
            index = i,
            pivot = pivot.use
          )
          if (length(index) == 1L && !isTRUE(progress$known_total)) {
            progress <- .np_npsig_progress_promote(
              progress, total = B, done = i.star
            )
          } else if (length(index) == 1L) {
            progress <- .np_progress_step(progress, done = i.star)
          } else {
            progress <- .np_progress_step(progress)
          }
        }
      } else {
        boot.seeds <- .npRmpi_npsig_bootstrap_seed_plan(
          num.obs = num.obs,
          boot.num = B,
          boot.method = boot.method,
          draw.wild.mult = draw.wild.mult,
          a = a,
          b = b,
          p.a = if (identical(boot.method, "wild-rademacher")) P.rademacher else P.a
        )
        post.boot.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        indiv.eval <- function(task.idx, seed.plan) {
          out <- numeric(length(task.idx))
          if (streamed.iid) {
            tile.width <- 8L
            for (tile.start in seq.int(1L, length(task.idx), by = tile.width)) {
              tile.count <- min(tile.width, length(task.idx) - tile.start + 1L)
              tile.position <- tile.start:(tile.start + tile.count - 1L)
              donor.index <- NULL
              response.matrix <- NULL
              if (identical(boot.method, "iid")) {
                donor.index <- vapply(
                  tile.position,
                  function(position) {
                    assign(".Random.seed", seed.plan[[task.idx[position]]],
                           envir = .GlobalEnv)
                    sample.int(num.obs, replace = TRUE)
                  },
                  integer(num.obs)
                )
              } else {
                wild.values <- if (identical(boot.method, "wild"))
                  c(a, b, P.a) else c(-1, 1, P.rademacher)
                response.matrix <- vapply(
                  tile.position,
                  function(position) {
                    assign(".Random.seed", seed.plan[[task.idx[position]]],
                           envir = .GlobalEnv)
                    mhat.xi + ei * draw.wild.mult(
                      num.obs, wild.values[[1L]], wild.values[[2L]],
                      wild.values[[3L]]
                    )
                  },
                  numeric(num.obs)
                )
              }
              out[tile.position] <- .np_npsig_streamed_iid_tile(
                bws = bws,
                xdat = xdat,
                tested.index = i,
                donor.index = donor.index,
                response.matrix = response.matrix,
                null.mean = mhat.xi,
                residual.pool = ei,
                pivotal = pivot.use
              )
            }
            return(out)
          }
          for (kk in seq_along(task.idx)) {
            assign(".Random.seed", seed.plan[[task.idx[kk]]], envir = .GlobalEnv)
            if (boot.method == "iid") {
              ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
              npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE,
                                                   se = pivot.use)
            } else if (boot.method == "wild") {
              ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
              npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE,
                                                   se = pivot.use)
            } else if (boot.method == "wild-rademacher") {
              ydat.star <- mhat.xi + ei * draw.wild.mult(
                num.obs, -1, 1, P.rademacher
              )
              npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE,
                                                   se = pivot.use)
            } else {
              boot.index <- sample.int(num.obs, replace = TRUE)
              ydat.star <- ydat[boot.index]
              xdat.star <- xdat
              xdat.star[, -i] <- xdat[boot.index, -i]
              if (direct.pairwise) {
                out[kk] <- .np_npsig_streamed_response_statistic(
                  bws = bws,
                  xdat = xdat.star,
                  index = i,
                  response.matrix = matrix(ydat.star, ncol = 1L),
                  pivotal = pivot.use
                )
              } else {
                npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                     txdat = xdat.star,
                                                     tydat = ydat.star,
                                                     bws = bws,
                                                     gradients = TRUE,
                                                     se = pivot.use)
              }
            }

            if (!direct.pairwise)
              out[kk] <- .np_npsig_statistic(
                npreg.boot,
                index = i,
                pivot = pivot.use
              )
          }
          out
        }

        indiv.bindings <- list(
          boot.method = boot.method,
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          i = i,
          pivot.use = pivot.use,
          .np_npsig_statistic = .np_npsig_statistic,
          num.obs = num.obs,
          draw.wild.mult = draw.wild.mult,
          a = a,
          b = b,
          P.a = P.a,
          extra.args = extra.args,
          streamed.iid = streamed.iid,
          direct.pairwise = direct.pairwise,
          .np_npsig_streamed_iid_tile = .np_npsig_streamed_iid_tile,
          .np_npsig_streamed_response_statistic =
            .np_npsig_streamed_response_statistic
        )
        if (boot.method != "pairwise")
          indiv.bindings <- c(indiv.bindings, list(mhat.xi = mhat.xi, ei = ei))

        if (length(index) == 1L) {
          progress$known_total <- TRUE
          progress$total <- B
          progress$throttle_sec <- .np_progress_interval_sec(
            known_total = TRUE,
            domain = progress$domain
          )
        }
        progress.context$state <- progress
        progress.context$done <- progress$last_done
        progress.context$use.bootstrap.done <- length(index) == 1L
        progress.context$force.next <- length(index) == 1L

        In.vec <- .npRmpi_npsig_parallel_boot_values(
          boot.seeds = boot.seeds,
          worker = indiv.eval,
          required.bindings = indiv.bindings,
          what = "npsigtest",
          profile.where = "npsigtest:indiv",
          progress.context = progress.context
        )
        assign(".Random.seed", post.boot.seed, envir = .GlobalEnv)
        progress <- progress.context$state
        if (length(index) == 1L) {
          progress <- .np_npsig_progress_promote(
            progress, total = B, done = B
          )
        }
      }
      
      ## Compute the P-value
      
      P[ii] <- .np_npsig_upper_tail_p(In.vec, In[ii])
      
      In.mat[,ii] = In.vec

      if (length(index) > 1L) {
        if (ii < length(index))
          progress$label <- paste("Testing", tested.names[[ii + 1L]])
        progress <- .np_npsig_progress_promote(
          progress, total = length(index), done = ii
        )
      }
      
    }
    
  } ## End invididual test

  progress <- .np_progress_end(progress)
  progress.active <- FALSE

  ## Return a list containing the statistic and its P-value
  ## bootstrapped In.vec for each variable...

  ## Restore seed

  .np_seed_exit(seed.state, remove_if_absent = TRUE)

  sigtest(In = In,
          In.bootstrap = In.mat,
          P = P,
          bws = bws,
          ixvar = index,
          boot.method = boot.method,
          pivot = pivot,
          pivot.effective = pivot.plan$effective,
          joint = joint,
          boot.type = boot.type,
          boot.num = B)

}

.np_npsig_default_test_args <- c(
  "B", "boot.method", "boot.type", "pivot", "joint", "index", "random.seed"
)

npsigtest.default <- function(bws, xdat, ydat, ...){
  .npRmpi_require_active_slave_pool(where = "npsigtest()")

  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat = if it
  ## was, we need to catch that and map it to dat = otherwise the call
  ## is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  xdat.named <- any(sc.names == "xdat")
  ydat.named <- any(sc.names == "ydat")

  no.bws <- missing(bws)
  no.xdat <- missing(xdat)
  no.ydat <- missing(ydat)

  ## if bws was passed in explicitly, do not compute bandwidths

  ## autodispatch normalizes calls via match.call(), which can turn an
  ## originally unnamed formula first argument into named bws=... .
  ## Preserve legacy formula behavior by rewriting npregbw() call shape.
  if (bws.named && no.xdat && no.ydat && inherits(bws, "formula")) {
    sc$`bws` <- NULL
    sc$formula <- bws
    sc.bw <- sc
    sc.bw[[1]] <- quote(npregbw)
    bws.named <- FALSE
  } else {
    sc.bw <- sc
    sc.bw[[1]] <- quote(npregbw)
  }
  sc.bw[.np_npsig_default_test_args] <- NULL

  if(xdat.named)
    xdat <- toFrame(xdat)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  tbw <- .np_eval_bw_call(sc.bw, caller_env = parent.frame())
  
  call.args <- list(bws = tbw)
  if(!no.xdat)
    call.args$xdat <- xdat
  if(!no.ydat)
    call.args$ydat <- ydat

  dots <- list(...)
  dots[c("bws", "bandwidth.compute", "formula", "data", "xdat", "ydat")] <- NULL

  ev <- do.call(npsigtest, c(call.args, dots))

  ev$call <- match.call(expand.dots = FALSE)
  environment(ev$call) <- parent.frame()
  return(ev)
}
