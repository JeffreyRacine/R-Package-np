npplregbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npplregbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat", "zdat"),
                                     eval_env = parent.frame())
    UseMethod("npplregbw", target)
  }

npplregbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    formula.call <- .np_bw_formula_from_call(call_obj = call, eval_env = parent.frame())
    if (!is.null(formula.call))
      mf[[2]] <- formula.call

    mf.xf <- mf
    
    mf[[1]] <- as.name("model.frame")
    mf.xf[[1]] <- as.name("model.frame")
    
    ## mangle formula ...
    formula.obj <- .np_bw_resolve_formula(formula_obj = formula,
                                        formula_call = formula.call,
                                        eval_env = parent.frame())
    chromoly <- explodePipe(formula.obj, env = environment(formula))

    if (length(chromoly) != 3) ## stop if malformed formula
      stop("invoked with improper formula, please see npplregbw documentation for proper use")

    ## make formula evaluable, then eval
    bronze <- lapply(chromoly, paste, collapse = " + ")

    mf.xf[["formula"]] <- as.formula(paste(" ~ ", bronze[[2]]),
                                    env = environment(formula))

    mf[["formula"]] <- as.formula(paste(bronze[[1]]," ~ ", bronze[[3]]),
                                  env = environment(formula))

    formula.all <- if(missing(data)) {
        terms(as.formula(paste(" ~ ",bronze[[1]]," + ",bronze[[2]], " + ",bronze[[3]]),
                                  env = environment(formula)))
    } else {
        terms(as.formula(paste(" ~ ",bronze[[1]]," + ",bronze[[2]], " + ",bronze[[3]]),
                                  env = environment(formula)), data = data)
    }

    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = formula.all,
                        data = environment(formula.all),
                        eval_env = environment(formula.all))
    else .np_terms_ts_mask(terms_obj = formula.all,
                           data = data,
                           eval_env = environment(formula.all))

    arguments.mfx <- chromoly[[2]]
    arguments.mf <- c(chromoly[[1]],chromoly[[3]])

    mf[["formula"]] <- terms(mf[["formula"]])
    mf.xf[["formula"]] <- terms(mf.xf[["formula"]])
    
    if(all(orig.ts)){
      arguments <- (as.list(attr(formula.all, "variables"))[-1])
      attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mf,arguments)),drop = FALSE])
      attr(mf.xf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mfx,arguments)),drop = FALSE])
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(formula.all, "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      attr(mf[["formula"]], "predvars") <- bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mf,arguments)),drop = FALSE])
      attr(mf.xf[["formula"]], "predvars") <- bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mfx,arguments)),drop = FALSE])
    }
    
    mf.args <- as.list(mf[-1L])
    mf.xf.args <- as.list(mf.xf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())
    mf.xf <- do.call(stats::model.frame, mf.xf.args, envir = parent.frame())

    ydat <- model.response(mf)
    xdat <- mf.xf
    zdat <- mf[, chromoly[[3]], drop = FALSE]

    tbw <- do.call(npplregbw,
                   c(list(xdat = xdat, ydat = ydat, zdat = zdat), list(...)))

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$xterms <- attr(mf.xf,"terms")
    tbw$chromoly <- chromoly

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)

    tbw
  }


npplregbw.NULL =
  function(xdat = stop("invoked without data `xdat'"),
           ydat = stop("invoked without data `ydat'"),
           zdat = stop("invoked without data `zdat'"),
           bws, ...){

    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain z names and 'toFrame'
    zdat <- toFrame(zdat)

    ## bandwidths

    bws = matrix(data = 0, nrow = 1+ncol(xdat), ncol = ncol(zdat))

    tbw <- npplregbw.default(xdat = xdat, ydat = ydat, zdat = zdat,
                             bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <-
      updateBwNameMetadata(nameList = list(ynames = .npplregbw_sanitize_yname(deparse(substitute(ydat)))),
                           bws = tbw)

    tbw
  }

npplregbw.plbandwidth =
  function(xdat = stop("invoked without data `xdat'"),
           ydat = stop("invoked without data `ydat'"),
           zdat = stop("invoked without data `zdat'"),
           bws, nmulti, ...){

    ## n observations
    ## rows are observations, columns are variables
    ## x is n x k
    ## y is n x 1
    ## z is n x p

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(dim(zdat)[2])
    }
    nmulti <- npValidateNmulti(nmulti)
    
    xdat = toFrame(xdat)
    zdat = toFrame(zdat)

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- na.action(na.omit(data.frame(xdat,ydat,zdat)))
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]
    zdat <- zdat[keep.rows,,drop = FALSE]

    ## y on z
    total.groups <- 1L + ncol(xdat)
    .np_progress_bandwidth_set_coordinator(
      total_groups = total.groups,
      local_total = nmulti
    )
    total.time <-
      system.time({
        .np_progress_bandwidth_set_coordinator_group(1L, "E[y|z]")
        bws$bw$yzbw  <- npregbw(xdat = zdat, ydat = ydat,
                                bws = bws$bw$yzbw, nmulti = nmulti, ...)
        
        ## x on z

        for (i in seq_len(ncol(xdat))) {
          .np_progress_bandwidth_set_coordinator_group(i + 1L, sprintf("E[%s|z]", names(xdat)[[i]]))
          bws$bw[[i+1]] <- npregbw(xdat=zdat, ydat=xdat[,i],
                  bws = bws$bw[[i+1]], nmulti = nmulti, ...)
        }
      })[["elapsed"]]
    num.feval <- sum(sapply(bws$bw, function(bwi) {
      if (is.null(bwi$num.feval) || identical(bwi$num.feval, NA)) 0 else bwi$num.feval
    }))
    num.feval.fast <- sum(sapply(bws$bw, function(bwi) {
      if (is.null(bwi$num.feval.fast) || identical(bwi$num.feval.fast, NA)) 0 else bwi$num.feval.fast
    }))
    fval <- {
      fv <- unlist(lapply(bws$bw, function(bwi) bwi$fval))
      if (length(fv) == 0 || all(!is.finite(fv))) NA else sum(fv[is.finite(fv)])
    }
    bws <- plbandwidth(bws = bws$bw,
                       regtype = bws$regtype,
                       basis = bws$basis,
                       degree = bws$degree,
                       bernstein.basis = bws$bernstein.basis,
                       bwmethod = bws$method,
                       bwscaling = bws$scaling,
                       bwtype = bws$type,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       ckerbound = bws$ckerbound,
                       ckerlb = bws$ckerlb,
                       ckerub = bws$ckerub,
                       ukertype = bws$ukertype,
                       okertype = bws$okertype,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       zdati = bws$zdati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       znames = bws$znames,
                       nobs = bws$nobs,
                       fval = fval,
                       num.feval = num.feval,
                       num.feval.fast = num.feval.fast,
                       rows.omit = rows.omit,
                       total.time = total.time)
    bws$degree.policy <- if (isTRUE(bws$child.degree.common)) "common-child-degree" else "child-specific"
    bws

  }

.npplregbw_sanitize_yname <- function(yname) {
  yname <- paste(yname, collapse = " ")
  if (!nzchar(yname) || nchar(yname) > 80L || grepl("^c\\(", yname))
    return("y")
  yname
}

.npplregbw_build_plbandwidth <- function(xdat,
                                         ydat,
                                         zdat,
                                         bws,
                                         bandwidth.compute,
                                         reg.args,
                                         outer.args,
                                         child.degree = NULL) {
  yname <- .npplregbw_sanitize_yname(deparse(substitute(ydat)))

  plband <- list()
  child.reg.args <- reg.args
  if (!is.null(child.degree))
    child.reg.args$degree <- child.degree[[1L]]
  plband$yzbw <- do.call(
    npregbw,
    c(list(xdat = zdat, ydat = ydat, bws = bws[1, ]), child.reg.args)
  )

  for (i in seq_len(dim(xdat)[2])) {
    child.reg.args <- reg.args
    if (!is.null(child.degree))
      child.reg.args$degree <- child.degree[[i + 1L]]
    plband[[i + 1L]] <- do.call(
      npregbw,
      c(list(xdat = zdat, ydat = xdat[, i], bws = bws[i + 1L, ]), child.reg.args)
    )
  }

  plbw.args <- c(
    list(
      bws = plband,
      nobs = dim(xdat)[1],
      fval = {
        fv <- unlist(lapply(plband, function(bwi) bwi$fval))
        if (length(fv) == 0 || all(!is.finite(fv))) NA else sum(fv[is.finite(fv)])
      },
      num.feval = sum(sapply(plband, function(bwi) {
        if (is.null(bwi$num.feval) || identical(bwi$num.feval, NA)) 0 else bwi$num.feval
      })),
      num.feval.fast = sum(sapply(plband, function(bwi) {
        if (is.null(bwi$num.feval.fast) || identical(bwi$num.feval.fast, NA)) 0 else bwi$num.feval.fast
      })),
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      zdati = untangle(zdat),
      xnames = names(xdat),
      ynames = yname,
      znames = names(zdat),
      bandwidth.compute = bandwidth.compute
    ),
    outer.args
  )

  tbw <- do.call(plbandwidth, plbw.args)
  if (!is.null(child.degree)) {
    degree.keys <- vapply(child.degree, paste, collapse = ",", character(1L))
    tbw$degree.policy <- if (length(unique(degree.keys)) <= 1L) "common-child-degree" else "child-specific"
  }
  tbw
}

.npplregbw_run_fixed_degree <- function(xdat,
                                        ydat,
                                        zdat,
                                        bws,
                                        reg.args,
                                        outer.args,
                                        opt.args,
                                        child.degree = NULL) {
  tbw <- .npplregbw_build_plbandwidth(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    bandwidth.compute = TRUE,
    reg.args = reg.args,
    outer.args = outer.args,
    child.degree = child.degree
  )

  do.call(
    npplregbw.plbandwidth,
    c(list(xdat = xdat, ydat = ydat, zdat = zdat, bws = tbw), opt.args)
  )
}

.npplregbw_child_responses <- function(xdat, ydat, yname = deparse(substitute(ydat))) {
  yname <- .npplregbw_sanitize_yname(yname)
  c(
    list(list(values = ydat, name = yname)),
    lapply(seq_len(ncol(xdat)), function(i) {
      list(values = xdat[, i], name = names(xdat)[i])
    })
  )
}

.npplregbw_child_degree_names <- function(xdat) {
  c("yzbw", names(toFrame(xdat)))
}

.npplregbw_validate_child_degree_entry <- function(value, ncon, child_name) {
  npValidateGlpDegree(
    regtype = "lp",
    degree = value,
    ncon = ncon,
    argname = sprintf("degree for npplreg child '%s'", child_name)
  )
}

.npplregbw_normalize_child_degree <- function(degree,
                                              regtype,
                                              ncon,
                                              child.names,
                                              degree.select) {
  out <- list(
    degree.for.spec = degree,
    child.degree = NULL,
    child.degree.common = TRUE
  )

  if (!identical(regtype, "lp") || is.null(degree))
    return(out)

  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  nchild <- length(child.names)

  if (is.list(degree) && !is.data.frame(degree)) {
    if (!identical(degree.select, "manual"))
      stop("child-specific npplreg degree lists are fixed-degree controls; use degree.select='manual'")
    if (length(degree) != nchild)
      stop(sprintf("npplreg degree list must have one entry per child regression (%d expected, got %d)",
                   nchild, length(degree)))
    if (!is.null(names(degree)) && all(nzchar(names(degree)))) {
      missing.names <- setdiff(child.names, names(degree))
      if (length(missing.names))
        stop(sprintf("npplreg degree list is missing child %s",
                     paste(sQuote(missing.names), collapse = ", ")))
      degree <- degree[child.names]
    }
    child.degree <- Map(
      function(value, child_name) .npplregbw_validate_child_degree_entry(value, ncon, child_name),
      degree,
      child.names
    )
  } else if (is.matrix(degree) || is.data.frame(degree)) {
    if (!identical(degree.select, "manual"))
      stop("child-specific npplreg degree matrices are fixed-degree controls; use degree.select='manual'")
    degree.matrix <- as.matrix(degree)
    if (nrow(degree.matrix) != nchild || ncol(degree.matrix) != ncon) {
      stop(sprintf(
        "npplreg degree matrix must have %d rows (yzbw plus x children) and %d columns (continuous z predictors)",
        nchild,
        ncon
      ))
    }
    if (!is.null(rownames(degree.matrix)) && all(nzchar(rownames(degree.matrix)))) {
      missing.names <- setdiff(child.names, rownames(degree.matrix))
      if (length(missing.names))
        stop(sprintf("npplreg degree matrix is missing child row %s",
                     paste(sQuote(missing.names), collapse = ", ")))
      degree.matrix <- degree.matrix[child.names, , drop = FALSE]
    }
    child.degree <- lapply(seq_len(nchild), function(i) {
      .npplregbw_validate_child_degree_entry(degree.matrix[i, ], ncon, child.names[[i]])
    })
    names(child.degree) <- child.names
  } else {
    common.degree <- .npplregbw_validate_child_degree_entry(degree, ncon, "all")
    child.degree <- replicate(nchild, common.degree, simplify = FALSE)
    names(child.degree) <- child.names
  }

  degree.keys <- vapply(child.degree, paste, collapse = ",", character(1L))
  out$degree.for.spec <- child.degree[[1L]]
  out$child.degree <- child.degree
  out$child.degree.common <- length(unique(degree.keys)) <= 1L
  out
}

.npplregbw_degree_search_controls <- function(regtype,
                                              regtype.named,
                                              bandwidth.compute,
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
  if (!isTRUE(bandwidth.compute))
    stop("automatic degree search requires bandwidth.compute=TRUE")
  if (ncon < 1L)
    stop("automatic degree search requires at least one continuous smoothing predictor")

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
      stop("degree.start must lie within the searched degree candidates for every continuous smoothing predictor")
    start.raw
  }

  list(
    method = if (identical(search.engine, "cell")) degree.select else search.engine,
    engine = search.engine,
    degree.select = degree.select,
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

.npplregbw_attach_degree_search <- function(bws, search_result) {
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
  bws <- .np_attach_nomad_restart_summary(bws, search_result)
  bws$degree.search <- metadata
  bws
}

.npplregbw_child_nomad_call_args <- function(zdat,
                                             response,
                                             bw.start,
                                             reg.args,
                                             opt.args,
                                             degree.search,
                                             nomad.inner.nmulti,
                                             random.seed) {
  kernel.arg.names <- c(
    "basis", "bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
    "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype",
    "scale.factor.search.lower"
  )
  optimizer.arg.names <- c(
    "ftol", "itmax", "nmulti", "nomad.remin", "powell.remin", "small", "tol",
    "scale.factor.search.lower"
  )
  child.args <- c(
    list(
      xdat = zdat,
      ydat = response,
      bws = bw.start,
      bandwidth.compute = TRUE,
      regtype = "lp",
      nomad = TRUE,
      search.engine = degree.search$engine,
      degree.select = degree.search$degree.select,
      bernstein.basis = degree.search$bernstein.basis,
      degree.min = degree.search$lower,
      degree.max = degree.search$upper,
      degree.start = degree.search$start.degree,
      degree.restarts = degree.search$restarts,
      degree.max.cycles = degree.search$max.cycles,
      degree.verify = degree.search$verify,
      nomad.nmulti = nomad.inner.nmulti,
      random.seed = random.seed
    ),
    reg.args[intersect(names(reg.args), kernel.arg.names)],
    opt.args[intersect(names(opt.args), optimizer.arg.names)]
  )
  child.args[c("degree", "degree.engine")] <- NULL
  child.args <- child.args[!duplicated(names(child.args), fromLast = TRUE)]
  child.args
}

.npplregbw_child_specific_nomad_search <- function(xdat,
                                                   ydat,
                                                   zdat,
                                                   bws,
                                                   reg.args,
                                                   outer.args,
                                                   opt.args,
                                                   degree.search,
                                                   nomad.inner.nmulti = 0L,
                                                   random.seed = 42L,
                                                   yname = deparse(substitute(ydat))) {
  if (isTRUE(degree.search$verify))
    stop("npplreg child-specific NOMAD does not support degree.verify")

  if (is.null(degree.search$degree.select) ||
      identical(degree.search$degree.select, "manual"))
    stop("npplreg child-specific NOMAD requires automatic child degree search")

  xdat <- toFrame(xdat)
  zdat <- toFrame(zdat)
  child.responses <- .npplregbw_child_responses(xdat = xdat, ydat = ydat, yname = yname)
  child.names <- c("yzbw", names(xdat))
  child.labels <- c("E[y|z]", sprintf("E[%s|z]", names(xdat)))
  child.list <- vector("list", length(child.responses))
  names(child.list) <- child.names

  total.start <- proc.time()[3L]
  .np_progress_bandwidth_set_coordinator(
    total_groups = length(child.responses),
    local_total = 1L
  )

  for (i in seq_along(child.responses)) {
    child.context <- sprintf("%s (%d/%d)", child.labels[[i]], i, length(child.responses))
    .np_progress_bandwidth_set_coordinator_group(i, child.context)
    old.context <- .np_progress_runtime$bandwidth_context_label
    .np_progress_bandwidth_set_context(child.context)
    child.list[[i]] <- tryCatch({
      child.args <- .npplregbw_child_nomad_call_args(
        zdat = zdat,
        response = child.responses[[i]]$values,
        bw.start = bws[i, ],
        reg.args = reg.args,
        opt.args = opt.args,
        degree.search = degree.search,
        nomad.inner.nmulti = nomad.inner.nmulti,
        random.seed = random.seed
      )
      do.call(npregbw, child.args)
    }, finally = {
      .np_progress_bandwidth_set_context(old.context)
    })
    child.list[[i]]$ynames <- child.responses[[i]]$name
  }

  total.time <- proc.time()[3L] - total.start
  child.fval <- vapply(child.list, function(bwi) {
    if (is.null(bwi$fval) || !length(bwi$fval)) NA_real_ else as.numeric(bwi$fval[1L])
  }, numeric(1L))
  child.baseline.fval <- vapply(child.list, function(bwi) {
    if (is.null(bwi$degree.search$baseline.fval) ||
        !length(bwi$degree.search$baseline.fval)) NA_real_ else as.numeric(bwi$degree.search$baseline.fval[1L])
  }, numeric(1L))
  child.feval <- vapply(child.list, function(bwi) {
    if (is.null(bwi$num.feval) || identical(bwi$num.feval, NA)) 0 else as.numeric(bwi$num.feval[1L])
  }, numeric(1L))
  child.feval.fast <- vapply(child.list, function(bwi) {
    if (is.null(bwi$num.feval.fast) || identical(bwi$num.feval.fast, NA)) 0 else as.numeric(bwi$num.feval.fast[1L])
  }, numeric(1L))
  child.nomad.time <- vapply(child.list, function(bwi) {
    if (is.null(bwi$nomad.time) || !length(bwi$nomad.time)) NA_real_ else as.numeric(bwi$nomad.time[1L])
  }, numeric(1L))
  child.powell.time <- vapply(child.list, function(bwi) {
    if (is.null(bwi$powell.time) || !length(bwi$powell.time)) NA_real_ else as.numeric(bwi$powell.time[1L])
  }, numeric(1L))

  common.degree <- {
    degree.keys <- vapply(child.list, function(bwi) paste(as.integer(bwi$degree), collapse = ","), character(1L))
    length(unique(degree.keys)) <= 1L
  }
  object.degree <- as.integer(child.list[[1L]]$degree)

  plbw.args <- c(
    list(
      bws = child.list,
      nobs = dim(xdat)[1],
      fval = if (all(!is.finite(child.fval))) NA_real_ else sum(child.fval[is.finite(child.fval)]),
      num.feval = sum(child.feval),
      num.feval.fast = sum(child.feval.fast),
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      zdati = untangle(zdat),
      xnames = names(xdat),
      ynames = child.responses[[1L]]$name,
      znames = names(zdat),
      bandwidth.compute = TRUE,
      total.time = total.time
    ),
    utils::modifyList(
      outer.args,
      list(
        regtype = "lp",
        degree = object.degree,
        bernstein.basis = degree.search$bernstein.basis
      )
    )
  )
  tbw <- do.call(plbandwidth, plbw.args)
  tbw$degree.policy <- if (isTRUE(common.degree)) "common-child-degree" else "child-specific"
  tbw$nomad.time <- if (all(!is.finite(child.nomad.time))) NA_real_ else sum(child.nomad.time[is.finite(child.nomad.time)])
  tbw$powell.time <- if (all(!is.finite(child.powell.time))) NA_real_ else sum(child.powell.time[is.finite(child.powell.time)])
  tbw$total.time <- sum(c(tbw$nomad.time, tbw$powell.time), na.rm = TRUE)
  tbw$child.fval <- child.fval
  tbw$child.baseline.fval <- child.baseline.fval
  tbw$child.num.feval <- child.feval
  tbw$child.num.feval.fast <- child.feval.fast
  tbw$child.nomad.time <- child.nomad.time
  tbw$child.powell.time <- child.powell.time

  selected.degrees <- lapply(child.list, function(bwi) as.integer(bwi$degree))
  names(selected.degrees) <- child.names
  baseline.degrees <- lapply(child.list, function(bwi) {
    if (!is.null(bwi$degree.search$baseline.degree)) as.integer(bwi$degree.search$baseline.degree) else integer(0L)
  })
  names(baseline.degrees) <- child.names

  list(
    method = degree.search$engine,
    direction = "min",
    verify = FALSE,
    completed = TRUE,
    certified = TRUE,
    interrupted = FALSE,
    baseline = list(
      degree = baseline.degrees,
      objective = if (all(!is.finite(child.baseline.fval))) NA_real_ else sum(child.baseline.fval[is.finite(child.baseline.fval)])
    ),
    best = list(
      degree = selected.degrees,
      objective = as.numeric(tbw$fval[1L]),
      child.fval = child.fval
    ),
    best_payload = tbw,
    nomad.time = tbw$nomad.time,
    powell.time = tbw$powell.time,
    optim.time = tbw$total.time,
    n.unique = NA_integer_,
    n.visits = NA_integer_,
    n.cached = NA_integer_,
    grid.size = NA_integer_,
    best.restart = NA_integer_,
    restart.starts = NULL,
    restart.degree.starts = NULL,
    restart.bandwidth.starts = NULL,
    restart.start.info = NULL,
    restart.results = NULL,
    trace = data.frame(
      child = child.names,
      label = child.labels,
      degree = vapply(selected.degrees, paste, collapse = ",", character(1L)),
      objective = child.fval,
      num.feval = child.feval,
      num.feval.fast = child.feval.fast,
      stringsAsFactors = FALSE
    )
  )
}

.npplregbw_nomad_search <- function(xdat,
                                    ydat,
                                    zdat,
                                    bws,
                                    reg.args,
                                    outer.args,
                                    opt.args,
                                    degree.search,
                                    nomad.inner.nmulti = 0L,
                                    random.seed = 42L) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template.outer.args <- outer.args
  template.outer.args$regtype <- "lp"
  template.outer.args$degree <- as.integer(degree.search$start.degree)
  template.outer.args$bernstein.basis <- degree.search$bernstein.basis

  template <- .npplregbw_build_plbandwidth(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args,
    outer.args = template.outer.args
  )

  child.templates <- template$bw
  if (!length(child.templates))
    stop("automatic degree search with search.engine='nomad' could not initialize child bandwidth templates")
  if (any(vapply(child.templates, function(tb) !identical(tb$type, "fixed"), logical(1))))
    stop("automatic degree search with search.engine='nomad' currently requires bwtype='fixed'")

  child.responses <- .npplregbw_child_responses(
    xdat = xdat,
    ydat = ydat,
    yname = template$ynames
  )
  child.setup <- lapply(child.templates, function(tb) .npregbw_nomad_bw_setup(xdat = zdat, template = tb))
  child.point.length <- vapply(seq_along(child.templates), function(i) {
    length(.npregbw_nomad_bw_to_point(child.templates[[i]]$bw, template = child.templates[[i]], setup = child.setup[[i]]))
  }, integer(1L))

  child_cont_lower <- function(i) {
    npGetScaleFactorSearchLower(
      child.templates[[i]],
      argname = "child scale.factor.search.lower"
    )
  }

  child.lower <- unlist(lapply(seq_along(child.setup), function(i) {
    c(rep.int(child_cont_lower(i), length(child.setup[[i]]$cont_idx)),
      rep.int(0, length(child.setup[[i]]$cat_idx)))
  }), use.names = FALSE)
  child.upper <- unlist(lapply(child.setup, function(setup) {
    c(rep.int(1e6, length(setup$cont_idx)), setup$cat_upper * setup$bandwidth.scale.categorical)
  }), use.names = FALSE)
  child.start <- unlist(lapply(seq_along(child.templates), function(i) {
    raw <- child.templates[[i]]$bw
    point.start <- if (all(raw == 0)) NULL else .npregbw_nomad_bw_to_point(raw, template = child.templates[[i]], setup = child.setup[[i]])
    .np_nomad_complete_start_point(
      point = point.start,
      lower = c(rep.int(child_cont_lower(i), length(child.setup[[i]]$cont_idx)),
                rep.int(0, length(child.setup[[i]]$cat_idx))),
      upper = c(rep.int(1e6, length(child.setup[[i]]$cont_idx)), child.setup[[i]]$cat_upper * child.setup[[i]]$bandwidth.scale.categorical),
      ncont = length(child.setup[[i]]$cont_idx)
    )
  }), use.names = FALSE)
  bwdim <- length(child.start)
  ndeg <- length(degree.search$start.degree)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(dim(zdat)[2]) else npValidateNmulti(opt.args$nmulti[1L])

  x0 <- c(child.start, as.integer(degree.search$start.degree))
  lb <- c(child.lower, degree.search$lower)
  ub <- c(child.upper, degree.search$upper)
  bbin <- c(rep.int(0L, bwdim), rep.int(1L, ndeg))
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  .np_nomad_baseline_note(degree.search$start.degree)

  point_to_matrix <- function(point) {
    cursor <- 1L
    out <- matrix(0, nrow = length(child.templates), ncol = ncol(zdat))
    for (i in seq_along(child.templates)) {
      seg <- point[cursor:(cursor + child.point.length[i] - 1L)]
      out[i, ] <- .npregbw_nomad_point_to_bw(seg, template = child.templates[[i]], setup = child.setup[[i]])
      cursor <- cursor + child.point.length[i]
    }
    out
  }

  child_eval_payload <- function(bw.matrix, degree, penalty.multiplier) {
    total.objective <- 0
    total.feval <- 0
    total.feval.fast <- 0
    child.payloads <- vector("list", length(child.templates))

    for (i in seq_along(child.templates)) {
      child.reg.args <- reg.args
      child.reg.args$regtype <- "lp"
      child.reg.args$degree <- degree
      child.reg.args$bernstein.basis <- degree.search$bernstein.basis
      child.reg.args$bandwidth.compute <- NULL
      resp <- child.responses[[i]]
      tbw <- .npregbw_build_rbandwidth(
        xdat = zdat,
        ydat = resp$values,
        bws = bw.matrix[i, ],
        bandwidth.compute = FALSE,
        reg.args = child.reg.args,
        yname = resp$name
      )
      out <- .npregbw_eval_only(
        xdat = zdat,
        ydat = resp$values,
        bws = tbw,
        invalid.penalty = "baseline",
        penalty.multiplier = penalty.multiplier
      )
      total.objective <- total.objective + out$objective
      total.feval <- total.feval + out$num.feval
      total.feval.fast <- total.feval.fast + out$num.feval.fast
      child.payloads[[i]] <- list(
        bws = tbw,
        objective = out$objective,
        response = resp,
        num.feval = out$num.feval,
        num.feval.fast = out$num.feval.fast
      )
    }

    list(
      objective = total.objective,
      num.feval = total.feval,
      num.feval.fast = total.feval.fast,
      child.payloads = child.payloads
    )
  }

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[bwdim + seq_len(ndeg)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw.matrix <- point_to_matrix(point[seq_len(bwdim)])
    out <- child_eval_payload(
      bw.matrix = bw.matrix,
      degree = degree,
      penalty.multiplier = 10
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
    bw.matrix <- point_to_matrix(point[seq_len(bwdim)])
    powell.elapsed <- NA_real_

    build_direct_payload <- function() {
      child.out <- child_eval_payload(
        bw.matrix = bw.matrix,
        degree = degree,
        penalty.multiplier = 10
      )
      child.list <- vector("list", length(child.templates))
      names(child.list) <- c("yzbw", names(xdat))

      for (i in seq_along(child.templates)) {
        resp <- child.responses[[i]]
        tbw <- child.out$child.payloads[[i]]$bws
        tbw$fval <- as.numeric(child.out$child.payloads[[i]]$objective)
        tbw$ifval <- as.numeric(child.out$child.payloads[[i]]$objective)
        tbw$num.feval <- as.numeric(child.out$child.payloads[[i]]$num.feval)
        tbw$num.feval.fast <- as.numeric(child.out$child.payloads[[i]]$num.feval.fast)
        tbw$fval.history <- as.numeric(child.out$child.payloads[[i]]$objective)
        tbw$eval.history <- 1
        tbw$invalid.history <- 0
        tbw$timing <- NA_real_
        tbw$total.time <- NA_real_
        child.list[[i]] <- npregbw.rbandwidth(
          xdat = zdat,
          ydat = resp$values,
          bws = tbw,
          bandwidth.compute = FALSE
        )
        if (!is.null(child.list[[i]]$method) && length(child.list[[i]]$method))
          child.list[[i]]$pmethod <- bwmToPrint(as.character(child.list[[i]]$method[1L]))
      }

      plbw.args <- c(
        list(
          bws = child.list,
          nobs = dim(xdat)[1],
          fval = child.out$objective,
          num.feval = as.numeric(nomad.num.feval.total),
          num.feval.fast = as.numeric(nomad.num.feval.fast.total),
          xdati = untangle(xdat),
          ydati = untangle(data.frame(ydat)),
          zdati = untangle(zdat),
          xnames = names(xdat),
          ynames = template$ynames,
          znames = names(zdat),
          bandwidth.compute = TRUE
        ),
        utils::modifyList(outer.args, list(regtype = "lp", degree = degree, bernstein.basis = degree.search$bernstein.basis))
      )
      do.call(plbandwidth, plbw.args)
    }

    direct.payload <- build_direct_payload()
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.outer.args <- outer.args
      hot.outer.args$regtype <- "lp"
      hot.outer.args$degree <- degree
      hot.outer.args$bernstein.basis <- degree.search$bernstein.basis
      hot.opt.args <- .np_nomad_powell_hotstart_opt_args(
        opt.args,
        strategy = "disable_multistart",
        remin = isTRUE(opt.args$powell.remin)
      )
      powell.start <- proc.time()[3L]
      hot.payload <- .np_nomad_with_powell_progress(
        degree = degree,
        best_record = best_record,
        expr = .npplregbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          zdat = zdat,
          bws = bw.matrix,
          reg.args = hot.reg.args,
          outer.args = hot.outer.args,
          opt.args = hot.opt.args
        )
      )
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

  .np_nomad_search(
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
    random.seed = random.seed,
    remin = isTRUE(opt.args$nomad.remin),
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
}

npplregbw.default = 
  function(xdat = stop("invoked without data `xdat'"),
           ydat = stop("invoked without data `ydat'"),
           zdat = stop("invoked without data `zdat'"),
           bandwidth.compute = TRUE,
           bws,
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
           scale.factor.search.lower = NULL,
           ftol, itmax, nmulti, nomad.remin, powell.remin, small, tol,
           ...){
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")

    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain z names and 'toFrame'
    zdat <- toFrame(zdat)

    dots <- list(...)
    dot.names <- names(dots)
    mc <- match.call(expand.dots = FALSE)
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
        regtype = if ("regtype" %in% dot.names) dots$regtype else NULL,
        search.engine = if ("search.engine" %in% mc.names) search.engine else NULL,
        degree.select = if ("degree.select" %in% mc.names) degree.select else NULL,
        bernstein.basis = if ("bernstein.basis" %in% dot.names) dots$bernstein.basis else NULL,
        degree.min = if ("degree.min" %in% mc.names) degree.min else NULL,
        degree.max = if ("degree.max" %in% mc.names) degree.max else NULL,
        degree.verify = if ("degree.verify" %in% mc.names) degree.verify else NULL,
        bwtype = if ("bwtype" %in% dot.names) dots$bwtype else NULL,
        degree = if (!missing(degree)) degree else NULL
      ),
      where = "npplregbw"
    )

    if (isTRUE(nomad.shortcut$enabled)) {
      if (!missing(degree))
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% dot.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% dot.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$bwtype, c("fixed", "generalized_nn", "adaptive_nn")))[1L], "fixed"))
        stop("nomad=TRUE currently requires bwtype='fixed'")
      if ("degree.select" %in% mc.names &&
          identical(as.character(match.arg(nomad.shortcut$values$degree.select, c("manual", "coordinate", "exhaustive")))[1L], "manual"))
        stop("nomad=TRUE requires automatic degree search; use degree.select='coordinate' or 'exhaustive'")
      if ("search.engine" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$search.engine, c("nomad+powell", "cell", "nomad")))[1L] %in%
              c("nomad", "nomad+powell")))
        stop("nomad=TRUE requires search.engine='nomad' or 'nomad+powell'")
      if ("bernstein.basis" %in% dot.names &&
          !isTRUE(npValidateGlpBernstein(regtype = "lp",
                                        bernstein.basis = nomad.shortcut$values$bernstein.basis)))
        stop("nomad=TRUE currently requires bernstein.basis=TRUE")
      if ("degree.verify" %in% mc.names &&
          isTRUE(npValidateScalarLogical(nomad.shortcut$values$degree.verify, "degree.verify")))
        stop("nomad=TRUE currently requires degree.verify=FALSE")
    }

    regtype.arg <- if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else "lc"
    basis.arg <- if ("basis" %in% dot.names) dots$basis else "glp"
    degree.arg <- if (!missing(degree)) degree else NULL
    bernstein.arg <- if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else FALSE
    ncon.z <- sum(untangle(zdat)$icon)
    child.degree.setup <- .npplregbw_normalize_child_degree(
      degree = degree.arg,
      regtype = regtype.arg,
      ncon = ncon.z,
      child.names = .npplregbw_child_degree_names(xdat),
      degree.select = if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    )
    degree.for.spec <- child.degree.setup$degree.for.spec

    spec.mc.names <- dot.names
    if (!missing(degree))
      spec.mc.names <- c(spec.mc.names, "degree")
    if (isTRUE(nomad.shortcut$enabled))
      spec.mc.names <- unique(c(spec.mc.names, "regtype", "bernstein.basis"))
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    degree.setup <- npSetupGlpDegree(
      regtype = regtype.arg,
      degree = degree.for.spec,
      ncon = ncon.z,
      degree.select = degree.select.value
    )
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)

    spec <- npResolveCanonicalConditionalRegSpec(
      mc.names = spec.mc.names,
      regtype = regtype.arg,
      basis = basis.arg,
      degree = degree.setup,
      bernstein.basis = bernstein.arg,
      ncon = ncon.z,
      where = "npplregbw"
    )

    random.seed.value <- .np_degree_extract_random_seed(dots)
    degree.search <- .npplregbw_degree_search_controls(
      regtype = regtype.arg,
      regtype.named = isTRUE(nomad.shortcut$enabled) || ("regtype" %in% dot.names),
      bandwidth.compute = bandwidth.compute,
      ncon = ncon.z,
      nobs = NROW(zdat),
      basis = basis.arg,
      degree.select = degree.select.value,
      search.engine = if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell",
      degree.min = nomad.shortcut$values$degree.min,
      degree.max = nomad.shortcut$values$degree.max,
      degree.start = if ("degree.start" %in% mc.names) degree.start else NULL,
      degree.restarts = if ("degree.restarts" %in% mc.names) degree.restarts else 0L,
      degree.max.cycles = if ("degree.max.cycles" %in% mc.names) degree.max.cycles else 20L,
      degree.verify = if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE,
      bernstein.basis = bernstein.arg,
      bernstein.named = isTRUE(nomad.shortcut$enabled) || ("bernstein.basis" %in% dot.names)
    )
    nomad.inner.named <- "nomad.nmulti" %in% mc.names
    nomad.inner.nmulti <- if (nomad.inner.named) {
      npValidateNonNegativeInteger(nomad.nmulti, "nomad.nmulti")
    } else {
      0L
    }
    if (nomad.inner.named &&
        (is.null(degree.search) || !(degree.search$engine %in% c("nomad", "nomad+powell")))) {
      stop("nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad+powell'")
    }

    reg.args <- list(
      regtype = spec$regtype.engine,
      basis = spec$basis.engine,
      degree = spec$degree.engine,
      bernstein.basis = spec$bernstein.basis.engine,
      bandwidth.compute = FALSE,
      scale.factor.search.lower = scale.factor.search.lower
    )
    if (!is.null(degree.search))
      reg.args$bernstein.basis <- degree.search$bernstein.basis
    kernel.arg.names <- intersect(
      c("bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
        "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype"),
      dot.names
    )
    if (length(kernel.arg.names))
      reg.args[kernel.arg.names] <- dots[kernel.arg.names]
    outer.args <- reg.args[setdiff(names(reg.args), "bandwidth.compute")]
    outer.args$regtype <- spec$regtype
    outer.args$basis <- spec$basis
    outer.args$degree <- spec$degree
    outer.args$bernstein.basis <- spec$bernstein.basis
    outer.args$scale.factor.search.lower <- scale.factor.search.lower

    opt.args <- list()
    margs <- c("regtype", "basis", "degree", "bernstein.basis",
               "bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype",
               "scale.factor.search.lower",
               "ftol", "itmax", "nmulti", "nomad.remin", "powell.remin", "small", "tol")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (bandwidth.compute) {
      if (!missing(nmulti)) nmulti <- npValidateNmulti(nmulti)
      if (!missing(powell.remin)) powell.remin <- npValidateScalarLogical(powell.remin, "powell.remin")
      if (!missing(nomad.remin)) nomad.remin <- npValidateScalarLogical(nomad.remin, "nomad.remin")
      if (!missing(itmax)) itmax <- npValidatePositiveInteger(itmax, "itmax")
      if (!missing(ftol)) ftol <- npValidatePositiveFiniteNumeric(ftol, "ftol")
      if (!missing(tol)) tol <- npValidatePositiveFiniteNumeric(tol, "tol")
      if (!missing(small)) small <- npValidatePositiveFiniteNumeric(small, "small")
      if (any.m) {
        nms <- mc.names[m]
        opt.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
      }
      opt.args$scale.factor.search.lower <- scale.factor.search.lower

      if (!is.null(degree.search)) {
        if (identical(degree.search$engine, "cell")) {
          eval_fun <- function(degree.vec) {
            cell.reg.args <- reg.args
            cell.reg.args$regtype <- "lp"
            cell.reg.args$degree <- as.integer(degree.vec)
            cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
            cell.outer.args <- outer.args
            cell.outer.args$regtype <- "lp"
            cell.outer.args$degree <- as.integer(degree.vec)
            cell.outer.args$bernstein.basis <- degree.search$bernstein.basis
            cell.bws <- .npplregbw_run_fixed_degree(
              xdat = xdat,
              ydat = ydat,
              zdat = zdat,
              bws = bws,
              reg.args = cell.reg.args,
              outer.args = cell.outer.args,
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
            direction = "min",
            trace_level = "full",
            objective_name = "fval"
          )
        } else {
          search.result <- .npplregbw_child_specific_nomad_search(
            xdat = xdat,
            ydat = ydat,
            zdat = zdat,
            bws = bws,
            reg.args = reg.args,
            outer.args = outer.args,
            opt.args = opt.args,
            degree.search = degree.search,
            nomad.inner.nmulti = nomad.inner.nmulti,
            random.seed = random.seed.value,
            yname = deparse(substitute(ydat))
          )
        }
        tbw <- .npplregbw_attach_degree_search(
          bws = search.result$best_payload,
          search_result = search.result
        )
        mc <- match.call(expand.dots = FALSE)
        environment(mc) <- parent.frame()
        tbw$call <- mc
        tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)
        return(tbw)
      }
    }

    tbw <- .npplregbw_build_plbandwidth(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      bws = bws,
      bandwidth.compute = bandwidth.compute,
      reg.args = reg.args,
      outer.args = outer.args,
      child.degree = child.degree.setup$child.degree
    )

    if (bandwidth.compute) {
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting partially linear regression bandwidth",
        do.call(npplregbw.plbandwidth, c(list(xdat = xdat, ydat = ydat, zdat = zdat, bws = tbw), opt.args))
      )
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)
  }
