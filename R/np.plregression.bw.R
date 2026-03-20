npplregbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
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

    tbw <- npplregbw(xdat = xdat, ydat = ydat, zdat = zdat, ...)

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
      updateBwNameMetadata(nameList = list(ynames = deparse(substitute(ydat))),
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
      nmulti <- min(5,dim(zdat)[2])
    }
    nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
    
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
      local_total = max(1L, nmulti)
    )
    total.time <-
      system.time({
        .np_progress_bandwidth_set_coordinator_group(1L, "y~z")
        bws$bw$yzbw  <- npregbw(xdat = zdat, ydat = ydat,
                                bws = bws$bw$yzbw, nmulti = nmulti, ...)
        
        ## x on z

        for (i in seq_len(ncol(xdat))) {
          .np_progress_bandwidth_set_coordinator_group(i + 1L, sprintf("x%s~z", i))
          bws$bw[[i+1]] <- npregbw(xdat=zdat, ydat=xdat[,i],
                  bws = bws$bw[[i+1]], nmulti = nmulti, ...)
        }
      })[["elapsed"]]
    num.feval <- sum(sapply(bws$bw, function(bwi) {
      if (is.null(bwi$num.feval) || identical(bwi$num.feval, NA)) 0 else bwi$num.feval
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
                       rows.omit = rows.omit,
                       total.time = total.time)

  }

.npplregbw_build_plbandwidth <- function(xdat,
                                         ydat,
                                         zdat,
                                         bws,
                                         bandwidth.compute,
                                         reg.args,
                                         outer.args) {
  yname <- deparse(substitute(ydat))

  plband <- list()
  plband$yzbw <- do.call(
    npregbw,
    c(list(xdat = zdat, ydat = ydat, bws = bws[1, ]), reg.args)
  )

  for (i in seq_len(dim(xdat)[2])) {
    plband[[i + 1L]] <- do.call(
      npregbw,
      c(list(xdat = zdat, ydat = xdat[, i], bws = bws[i + 1L, ]), reg.args)
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

  do.call(plbandwidth, plbw.args)
}

.npplregbw_run_fixed_degree <- function(xdat, ydat, zdat, bws, reg.args, outer.args, opt.args) {
  tbw <- .npplregbw_build_plbandwidth(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    bandwidth.compute = TRUE,
    reg.args = reg.args,
    outer.args = outer.args
  )

  do.call(
    npplregbw.plbandwidth,
    c(list(xdat = xdat, ydat = ydat, zdat = zdat, bws = tbw), opt.args)
  )
}

.npplregbw_child_responses <- function(xdat, ydat, yname = deparse(substitute(ydat))) {
  c(
    list(list(values = ydat, name = yname)),
    lapply(seq_len(ncol(xdat)), function(i) {
      list(values = xdat[, i], name = names(xdat)[i])
    })
  )
}

.npplregbw_degree_search_controls <- function(regtype,
                                              regtype.named,
                                              bandwidth.compute,
                                              ncon,
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

  if (!isTRUE(bern.auto) && any(bounds$upper > 3L))
    stop("automatic degree search with bernstein.basis=FALSE currently requires degree.max <= 3")

  baseline.degree <- rep.int(0L, ncon)
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, baseline.degree))
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
    candidates = bounds$candidates,
    lower = bounds$lower,
    upper = bounds$upper,
    baseline.degree = baseline.degree,
    start.degree = start.degree,
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
    restart.starts = lapply(search_result$restart.starts, as.integer),
    trace = search_result$trace
  )

  if (!is.null(search_result$nomad.time))
    bws$nomad.time <- as.numeric(search_result$nomad.time[1L])
  if (!is.null(search_result$powell.time))
    bws$powell.time <- as.numeric(search_result$powell.time[1L])
  if (!is.null(search_result$optim.time) && is.finite(search_result$optim.time))
    bws$total.time <- as.numeric(search_result$optim.time[1L])
  bws$degree.search <- metadata
  bws
}

.npplregbw_nomad_search <- function(xdat,
                                    ydat,
                                    zdat,
                                    bws,
                                    reg.args,
                                    outer.args,
                                    opt.args,
                                    degree.search) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  baseline.reg.args <- reg.args
  baseline.reg.args$regtype <- "lp"
  baseline.reg.args$degree <- as.integer(degree.search$baseline.degree)
  baseline.reg.args$bernstein.basis <- degree.search$bernstein.basis
  baseline.outer.args <- outer.args
  baseline.outer.args$regtype <- "lp"
  baseline.outer.args$degree <- as.integer(degree.search$baseline.degree)
  baseline.outer.args$bernstein.basis <- degree.search$bernstein.basis

  .np_nomad_baseline_note(degree.search$baseline.degree)
  baseline.bws <- .npplregbw_run_fixed_degree(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    reg.args = baseline.reg.args,
    outer.args = baseline.outer.args,
    opt.args = opt.args
  )

  child.templates <- baseline.bws$bw
  if (!length(child.templates))
    stop("automatic degree search with search.engine='nomad' could not initialize child bandwidth templates")
  if (any(vapply(child.templates, function(tb) !identical(tb$type, "fixed"), logical(1))))
    stop("automatic degree search with search.engine='nomad' currently requires bwtype='fixed'")

  child.responses <- .npplregbw_child_responses(
    xdat = xdat,
    ydat = ydat,
    yname = baseline.bws$ynames
  )
  child.setup <- lapply(child.templates, function(tb) .npregbw_nomad_bw_setup(xdat = zdat, template = tb))
  child.point.length <- vapply(seq_along(child.templates), function(i) {
    length(.npregbw_nomad_bw_to_point(child.templates[[i]]$bw, template = child.templates[[i]], setup = child.setup[[i]]))
  }, integer(1L))

  baseline.record <- list(
    eval_id = 0L,
    degree = as.integer(degree.search$baseline.degree),
    objective = as.numeric(baseline.bws$fval[1L]),
    status = "ok",
    cached = FALSE,
    message = NULL,
    elapsed = 0,
    num.feval = if (!is.null(baseline.bws$num.feval)) as.numeric(baseline.bws$num.feval[1L]) else NA_real_
  )

  child.lower <- unlist(lapply(child.setup, function(setup) {
    c(rep.int(1e-2, length(setup$cont_idx)), rep.int(0, length(setup$cat_idx)))
  }), use.names = FALSE)
  child.upper <- unlist(lapply(child.setup, function(setup) {
    c(rep.int(1e6, length(setup$cont_idx)), setup$cat_upper * setup$bandwidth.scale.categorical)
  }), use.names = FALSE)
  child.start <- unlist(lapply(seq_along(child.templates), function(i) {
    .npregbw_nomad_bw_to_point(child.templates[[i]]$bw, template = child.templates[[i]], setup = child.setup[[i]])
  }), use.names = FALSE)
  bwdim <- length(child.start)
  ndeg <- length(degree.search$baseline.degree)

  x0 <- c(child.start, as.integer(degree.search$start.degree))
  lb <- c(child.lower, degree.search$lower)
  ub <- c(child.upper, degree.search$upper)
  bbin <- c(rep.int(0L, bwdim), rep.int(1L, ndeg))

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
      child.payloads[[i]] <- list(bws = tbw, objective = out$objective, response = resp)
    }

    list(objective = total.objective, num.feval = total.feval, child.payloads = child.payloads)
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
    list(
      objective = out$objective,
      degree = degree,
      num.feval = out$num.feval
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

      for (i in seq_along(child.templates)) {
        resp <- child.responses[[i]]
        tbw <- child.out$child.payloads[[i]]$bws
        tbw$fval <- as.numeric(child.out$child.payloads[[i]]$objective)
        tbw$ifval <- as.numeric(child.out$child.payloads[[i]]$objective)
        tbw$num.feval <- 1
        tbw$num.feval.fast <- 0
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
      }

      plbw.args <- c(
        list(
          bws = child.list,
          nobs = dim(xdat)[1],
          fval = child.out$objective,
          num.feval = if (!is.null(solution$bbe)) as.numeric(solution$bbe) else child.out$num.feval,
          xdati = untangle(xdat),
          ydati = untangle(data.frame(ydat)),
          zdati = untangle(zdat),
          xnames = names(xdat),
          ynames = baseline.bws$ynames,
          znames = names(zdat),
          bandwidth.compute = TRUE
        ),
        utils::modifyList(outer.args, list(regtype = "lp", degree = degree, bernstein.basis = degree.search$bernstein.basis))
      )
      do.call(plbandwidth, plbw.args)
    }

    use.baseline.payload <- identical(as.integer(degree), as.integer(degree.search$baseline.degree))
    direct.payload <- if (isTRUE(use.baseline.payload)) baseline.bws else build_direct_payload()
    direct.objective <- if (isTRUE(use.baseline.payload)) {
      as.numeric(baseline.record$objective)
    } else {
      as.numeric(best_record$objective)
    }

    if (identical(degree.search$engine, "nomad+powell")) {
      .np_nomad_powell_note(degree)
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.outer.args <- outer.args
      hot.outer.args$regtype <- "lp"
      hot.outer.args$degree <- degree
      hot.outer.args$bernstein.basis <- degree.search$bernstein.basis
      hot.opt.args <- opt.args
      hot.opt.args$nmulti <- 0L
      powell.start <- proc.time()[3L]
      hot.payload <- .npplregbw_run_fixed_degree(
        xdat = xdat,
        ydat = ydat,
        zdat = zdat,
        bws = bw.matrix,
        reg.args = hot.reg.args,
        outer.args = hot.outer.args,
        opt.args = hot.opt.args
      )
      powell.elapsed <- proc.time()[3L] - powell.start
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
    x0 = x0,
    bbin = bbin,
    lb = lb,
    ub = ub,
    eval_fun = eval_fun,
    build_payload = build_payload,
    direction = "min",
    objective_name = "fval"
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
           degree.min = NULL,
           degree.max = NULL,
           degree.start = NULL,
           degree.restarts = 0L,
           degree.max.cycles = 20L,
           degree.verify = FALSE,
           ftol, itmax, nmulti, remin, small, tol,
           ...){
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")

    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain z names and 'toFrame'
    zdat <- toFrame(zdat)

    dots <- list(...)
    dot.names <- names(dots)
    regtype.arg <- if ("regtype" %in% dot.names) dots$regtype else "lc"
    basis.arg <- if ("basis" %in% dot.names) dots$basis else "glp"
    degree.arg <- if (!missing(degree)) degree else NULL
    bernstein.arg <- if ("bernstein.basis" %in% dot.names) dots$bernstein.basis else FALSE

    spec.mc.names <- dot.names
    if (!missing(degree))
      spec.mc.names <- c(spec.mc.names, "degree")

    spec <- npResolveCanonicalConditionalRegSpec(
      mc.names = spec.mc.names,
      regtype = regtype.arg,
      basis = basis.arg,
      degree = degree.arg,
      bernstein.basis = bernstein.arg,
      ncon = sum(untangle(zdat)$icon),
      where = "npplregbw"
    )

    mc.names <- names(match.call(expand.dots = FALSE))
    degree.search <- .npplregbw_degree_search_controls(
      regtype = regtype.arg,
      regtype.named = "regtype" %in% dot.names,
      bandwidth.compute = bandwidth.compute,
      ncon = sum(untangle(zdat)$icon),
      degree.select = if ("degree.select" %in% mc.names) degree.select else "manual",
      search.engine = if ("search.engine" %in% mc.names) search.engine else "nomad+powell",
      degree.min = if ("degree.min" %in% mc.names) degree.min else NULL,
      degree.max = if ("degree.max" %in% mc.names) degree.max else NULL,
      degree.start = if ("degree.start" %in% mc.names) degree.start else NULL,
      degree.restarts = if ("degree.restarts" %in% mc.names) degree.restarts else 0L,
      degree.max.cycles = if ("degree.max.cycles" %in% mc.names) degree.max.cycles else 20L,
      degree.verify = if ("degree.verify" %in% mc.names) degree.verify else FALSE,
      bernstein.basis = bernstein.arg,
      bernstein.named = "bernstein.basis" %in% dot.names
    )

    reg.args <- list(
      regtype = spec$regtype.engine,
      basis = spec$basis.engine,
      degree = spec$degree.engine,
      bernstein.basis = spec$bernstein.basis.engine,
      bandwidth.compute = FALSE
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

    opt.args <- list()
    margs <- c("regtype", "basis", "degree", "bernstein.basis",
               "bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype",
               "ftol", "itmax", "nmulti", "remin", "small", "tol")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (bandwidth.compute) {
      if (!missing(nmulti)) nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
      if (!missing(remin)) remin <- npValidateScalarLogical(remin, "remin")
      if (!missing(itmax)) itmax <- npValidatePositiveInteger(itmax, "itmax")
      if (!missing(ftol)) ftol <- npValidatePositiveFiniteNumeric(ftol, "ftol")
      if (!missing(tol)) tol <- npValidatePositiveFiniteNumeric(tol, "tol")
      if (!missing(small)) small <- npValidatePositiveFiniteNumeric(small, "small")
      if (any.m) {
        nms <- mc.names[m]
        opt.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
      }

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
          search.result <- .npplregbw_nomad_search(
            xdat = xdat,
            ydat = ydat,
            zdat = zdat,
            bws = bws,
            reg.args = reg.args,
            outer.args = outer.args,
            opt.args = opt.args,
            degree.search = degree.search
          )
        }
        tbw <- .npplregbw_attach_degree_search(
          bws = search.result$best_payload,
          search_result = search.result
        )
        mc <- match.call(expand.dots = FALSE)
        environment(mc) <- parent.frame()
        tbw$call <- mc
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
      outer.args = outer.args
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

    return(tbw)
  }
