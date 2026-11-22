.indextest_recipe <- function(bw) {
  # These three controls survive evaluation in the bandwidth object. Reading
  # the original call's symbols could pick up a value changed after the fit.
  starts <- length(bw[["fval.vector"]])
  optimizer <- bw[["pomethod"]]
  floor <- bw[["scale.factor.search.lower"]]
  if (!starts || length(optimizer) != 1L || !optimizer %in% c("Nelder-Mead", "BFGS", "CG") ||
      length(floor) != 1L || !is.finite(floor))
    stop("The fitted object lacks the joint search metadata. Refit from x and y.")
  if (!bw[["pmethod"]] %in% c("Ichimura", "Klein and Spady"))
    stop("Bootstrap refits require an original joint index and bandwidth search.")
  controls <- list(optim.method = optimizer, scale.factor.search.lower = floor)
  call <- bw[["call"]]
  if (!is.call(call)) stop("Original search call unavailable. Refit from x and y.")
  args <- as.list(call)[-1L]
  stored <- c("formula","data","subset","na.action","xdat","ydat","method","regtype",
              "bwtype","ckertype","ckerorder","ckerbound","ckerlb","ckerub",
              "nmulti","optim.method","scale.factor.search.lower","only.optimize.beta",
              "bandwidth.compute")
  replay <- c("bws","optim.maxit","optim.maxattempts","optim.reltol","optim.abstol",
              "random.seed","scale.factor.init.lower","scale.factor.init.upper",
              "scale.factor.init")
  unhandled <- setdiff(names(args), c(stored,replay))
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))) || length(unhandled)))
    stop("Cannot reconstruct this original search call. Use a supported, named npindexbw call.")
  literal <- function(expr) {
    if (is.atomic(expr)) return(TRUE)
    is.call(expr) && as.character(expr[[1L]])[1L] %in% c("c", "-", "+") &&
      all(vapply(as.list(expr)[-1L], literal, logical(1)))
  }
  for (key in intersect(replay,names(args))) {
    expr <- args[[key]]
    if (!literal(expr)) stop("Original ",key," value is unavailable. Refit with literal controls or do.call using evaluated controls.")
    controls[key] <- list(eval(expr, baseenv()))
  }
  c(list(method = bw[["method"]], nmulti = starts, regtype = "lc",
         bwtype = "fixed", ckertype = "gaussian"), controls)
}

.indextest_setup <- function(object, bw.method) {
  bw <- if (inherits(object, "singleindex")) object[["bws"]] else object
  if (!inherits(bw, "sibandwidth"))
    stop("object must be an npindex() fit or npindexbw() bandwidth object.", call. = FALSE)
  if (!identical(bw[["regtype"]], "lc") || !identical(bw[["type"]], "fixed") ||
      !identical(bw[["ckertype"]], "gaussian") || !identical(bw[["ckerbound"]], "none") ||
      !isTRUE(bw[["ckerorder"]] == 2))
    stop("The test requires a fixed, unbounded, second-order Gaussian local-constant fit.", call. = FALSE)
  if (length(bw[["fval"]]) != 1L || !is.finite(bw[["fval"]]))
    stop("The test requires a joint index and bandwidth search.", call. = FALSE)
  if (!is.null(bw[["index.design"]]))
    stop("The test requires continuous numeric covariates.", call. = FALSE)
  recipe <- .indextest_recipe(bw)
  training <- bw[[".np.native.training"]]
  if (is.null(training)) {
    frame <- bw[[".np.formula.training"]][["frame"]]
    if (is.null(frame))
      stop("Original training data unavailable. Refit the index model before testing.", call. = FALSE)
    training <- list(xdat = .np_index_formula_xdat(frame), ydat = stats::model.response(frame))
  }
  x <- as.data.frame(training$xdat)
  y <- training$ydat
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) != nrow(x) ||
      ncol(x) < 2L || !all(vapply(x, function(v) is.numeric(v) && is.null(dim(v)), logical(1))))
    stop("At least two continuous numeric covariates and an aligned numeric response are required.", call. = FALSE)
  keep <- complete.cases(x, y)
  x <- x[keep, , drop = FALSE]
  y <- as.double(y[keep])
  n <- length(y)
  if (n <= ncol(x)+2L || n != bw[["nobs"]] || any(!is.finite(as.matrix(x))) ||
      any(!is.finite(y)) || length(unique(y)) < 2L)
    stop("The joint search must have a matching, finite, nonconstant complete training sample.", call. = FALSE)
  if (recipe$method == "kleinspady" && !all(y %in% c(0, 1)))
    stop("Klein-Spady requires a numeric 0/1 response.", call. = FALSE)
  # A held-bandwidth npindex fit can replace its training data. Check the
  # original search criterion once before using that direction for inference.
  if (inherits(object, "singleindex")) {
    checked <- .np_plot_with_local_compiled_eval(.npindexbw_certify_selected_candidate(c(bw[["beta"]][-1L], bw[["bw"]]),
      as.matrix(x), y, bw, .npindex_resolve_spec(bw)))
    if (!isTRUE(all.equal(as.double(checked$objective), as.double(bw[["fval"]]),
                          tolerance = 1e-8)))
      stop("The stored index search does not match the fit's training sample. Refit npindexbw on these data.", call. = FALSE)
  }
  centered <- sweep(as.matrix(x), 2L, colMeans(x))
  scales <- sqrt(colMeans(centered^2))
  if (any(!is.finite(scales)) || any(scales <= 0))
    stop("Rescale the covariates to have finite positive variance.", call. = FALSE)
  s <- sweep(centered, 2L, scales, "/")
  if (qr(s)$rank != ncol(s))
    stop("The centered covariates must have full column rank.", call. = FALSE)
  h <- switch(bw.method, published = n^(-2/9), maximum = c(.5,1,2)*n^(-2/9),
              separate = n^(-17/36))
  # Recreate the specified starts, including automatic zero placeholders.
  # Reusing the observed winning coefficients would change the bootstrap.
  initial <- recipe[["bws"]]
  if (is.null(initial)) initial <- numeric(ncol(x) + 1L)
  template <- sibandwidth(beta = initial[seq_len(ncol(x))], h = tail(initial, 1L),
    method = recipe$method, nobs = n, xdati = bw$xdati, ydati = bw$ydati,
    xnames = names(x), ynames = bw$ynames, bandwidth = tail(initial, 1L),
    bandwidth.compute = TRUE)
  template <- npSetScaleFactorSearchLower(template, recipe$scale.factor.search.lower)
  descriptor <- .npindexbw_prepare_lp_regression_leaf_descriptor(
    x[[1L]], y, 1, template, .npindex_resolve_spec(template))
  cv <- .npregbw_call_fixed_degree_core(xdat = x[,1L,drop=FALSE], ydat = y,
    bws = descriptor$bws, invalid.penalty = "dbmax", eval.only = TRUE,
    objective = if (recipe$method == "ichimura") "ls" else "ks", .prepare.index = TRUE)
  list(x=x, y=y, s=s, scales=scales, template=template, cv=cv, recipe=recipe, h=h,
       g.test=if (bw.method=="separate") n^(-2/9), beta=bw[["beta"]], bw=bw[["bw"]])
}

.indextest_index <- function(data, beta, bw) {
  direction <- beta*data$scales
  scale <- sqrt(sum(direction^2))
  g <- as.double(bw)/scale
  z <- drop(data$s %*% (direction/scale))
  if (length(g)!=1L || !is.finite(g) || g<=0 || any(!is.finite(z)))
    stop("The fitted index or link bandwidth is invalid.", call. = FALSE)
  list(z=z, g=g)
}

.indextest_score <- function(y, index, data) {
  g <- if (is.null(data$g.test)) index$g else data$g.test
  if (is.null(.np_progress_runtime$fit_forward))
    return(.Call(C_np_indextest_score, y, index$z, data$s, g, data$h))
  .np_with_compiled_fit_progress("Computing the single-index statistic", length(y),
    expr = .Call(C_np_indextest_score, y, index$z, data$s, g, data$h))
}

.indextest_pilot <- function(data, index) {
  y <- data$y
  n <- length(y)
  pilot <- numeric(n)
  binary <- data$recipe$method == "kleinspady"
  # Retain the original rowwise normalization and matrix product. Only a
  # fixed number of target rows is live at any point.
  for (first in seq.int(1L, n, by = 64L)) {
    .np_fit_progress_step()
    rows <- seq.int(first, min(n, first + 63))
    logw <- -0.5*(outer(index$z[rows], index$z, "-")/index$g)^2
    if (binary) {
      weights <- exp(logw)
      success <- rowSums(weights[,y == 1,drop = FALSE])
      failure <- rowSums(weights[,y == 0,drop = FALSE])
      pilot[rows] <- success/(success + failure)
    } else {
      logw[cbind(seq_along(rows), rows)] <- -Inf
      weights <- exp(logw - apply(logw, 1L, max))
      pilot[rows] <- drop(weights %*% y)/rowSums(weights)
    }
  }
  pilot
}

.indextest_response <- function(data) {
  n <- length(data$y)
  if (data$recipe$method == "kleinspady")
    return(as.double(stats::rbinom(n, 1L, data$pilot)))
  eta <- ifelse(stats::runif(n) < (5+sqrt(5))/10, (1-sqrt(5))/2, (1+sqrt(5))/2)
  data$pilot + eta*(data$y-data$pilot)
}

.indextest_cv_leaf <- function(cv) {
  force(cv)
  heartbeat <- .np_progress_runtime$fit_forward
  function(index, ydat, h, bws, spec, leaf.descriptor) {
    if (!is.null(heartbeat)) heartbeat()
    out <- .Call(C_np_regression_bw_eval, numeric(), numeric(), as.double(index),
      ydat, as.double(EssDee(index)), cv$opti, cv$optd, as.double(h), 1L, 0L, 10,
      cv$degree, cv$bernstein, cv$basis, cv$lo, cv$hi, cv$support)
    list(objective = out$fval[1L], num.feval.fast = out$fast.history[1L])
  }
}

.indextest_refit <- function(y, data, draw) {
  # Generate this response in the bootstrap RNG scope, before the search resets its seed.
  force(y)
  tryCatch({
    controls <- data$recipe[setdiff(names(data$recipe),
        c("method","bwtype","regtype","ckertype","bws"))]
    fit <- do.call(npindexbw.sibandwidth, c(list(xdat = quote(data$x), ydat = quote(y),
      bws = quote(data$template), .indextest.leaf = .indextest_cv_leaf(data$cv)), controls))
    ix <- .indextest_index(data, fit[["beta"]], fit[["bw"]])
    max(.indextest_score(y, ix, data))
  }, error = function(e) stop("Single-index bootstrap draw ", draw, ": ",
                              conditionMessage(e), call. = FALSE))
}

npindextest <- function(object, bw.method = c("published", "maximum", "separate"),
                        B = 199L, random.seed = 42L) {
  started <- proc.time()[["elapsed"]]
  call <- match.call()
  absent <- missing(object)
  .pltest_local({
    if (absent) stop("Supply an npindex() fit or npindexbw() bandwidth object.", call. = FALSE)
    bw.method <- match.arg(bw.method, c("published", "maximum", "separate"))
    if (!is.null(dim(B)) || is.object(B) || (is.numeric(B) && any(B > .Machine$integer.max, na.rm=TRUE)))
      stop("B must be one representable positive integer.", call. = FALSE)
    B <- npValidatePositiveInteger(B, "B")
    if (B < 9L) stop("B must be at least 9.", call. = FALSE)
    if (!is.null(random.seed)) {
      if (!is.null(dim(random.seed)) || is.object(random.seed) ||
          (is.numeric(random.seed) && any(random.seed > .Machine$integer.max, na.rm=TRUE)))
        stop("random.seed must be NULL or one representable nonnegative integer.", call. = FALSE)
      random.seed <- npValidateNonNegativeInteger(random.seed, "random.seed")
    }
    data <- .np_progress_activity_run("Preparing the single-index test",
      .indextest_setup(object, bw.method))
    index <- .indextest_index(data, data$beta, data$bw)
    scores <- .np_progress_activity_run("Computing the single-index statistic",
      .indextest_score(data$y, index, data))
    observed <- max(scores)
    data$pilot <- .np_progress_activity_run("Computing the bootstrap pilot",
      .indextest_pilot(data, index))
    if (any(!is.finite(data$pilot))) stop("The bootstrap pilot is not finite.", call. = FALSE)
    .pltest_preflight(B, "npindextest")
  }, where="npindextest")
  .pltest_same(list(bw.method,B,random.seed,data$x,data$y,data$recipe,data$beta,data$bw))
  if (is.null(random.seed)) {
    .pltest_current_stream()
    on.exit(.pltest_current_stream(), add=TRUE)
  } else {
    state <- .pltest_rng_capture()
    on.exit(.np_rng_scope_restore(state), add=TRUE)
    .np_rng_scope_seed(state, random.seed)
  }
  boot <- .indextest_bootstrap(data, B)
  values <- boot$statistics
  if (length(values)!=B || any(!is.finite(values)))
    stop("Expected B finite single-index bootstrap statistics.", call. = FALSE)
  exceed <- sum(values >= observed)
  structure(list(call=call, model.call=object[["call"]], method=data$recipe$method,
    bw.method=bw.method, statistic=observed, components=scores, P=(1+exceed)/(1+as.double(B)),
    exceedances=exceed, critical.values=stats::quantile(values,c(.9,.95,.99),type=1),
    bootstrap=values, rank.counts=boot$rank.counts, boot.num=B, random.seed=random.seed,
    n=length(data$y), beta=data$beta, refit.controls=data$recipe,
    g.cv=index$g, g.test=if (is.null(data$g.test)) index$g else data$g.test, h=data$h,
    total.time = proc.time()[["elapsed"]] - started),
    class="indextest")
}

print.indextest <- function(x, ...) {
  cat("Maistre-Patilea single-index specification test\n")
  .np_spec_test_print_result(x)
  invisible(x)
}

summary.indextest <- function(object, ...) {
  print(object)
  cat("Index estimator:", object$method, "   Bandwidth method:", object$bw.method,
      "\nCV link bandwidth:", format(object$g.cv), "   Residual bandwidth:", format(object$g.test),
      "\nPairwise bandwidths:", format(object$h),
      "\nBootstrap draws:", object$boot.num, "   Search starts per refit:",
      object$refit.controls$nmulti, "\n")
  if (!is.null(object[["total.time"]]))
    cat("Computation Time: ", format(object[["total.time"]]), " seconds\n", sep = "")
  invisible(object)
}

.indextest_local_refits <- function(data, first, size, responses = NULL) {
  .npRmpi_with_local_bootstrap(.np_plot_with_local_compiled_eval({
    vapply(seq_len(size), function(j) {
      y <- if (is.null(responses)) .indextest_response(data) else responses[,j]
      .indextest_refit(y, data, first+j-1)
    }, numeric(1))
  }))
}

.indextest_bootstrap <- function(data, B) {
  if (.npRmpi_cms_collective_context()) {
    rank <- mpi.comm.rank(1L)
    size <- mpi.comm.size(1L)
    values <- numeric(B)
    counts <- integer(size)
    progress <- .pltest_local(if (rank==0L) .np_progress_activity_begin(B, "Bootstrapping the single-index test", detail = NULL))
    complete <- FALSE
    on.exit(if (rank == 0L) .np_progress_activity_end(progress, complete), add=TRUE)
    for (first in seq.int(1, B, by=32)) {
      ids <- seq.int(first, min(as.double(B),first+31))
      responses <- .pltest_local(if (rank==0L)
        matrix(unlist(lapply(ids, function(i) .indextest_response(data)),use.names=FALSE),
               nrow=length(data$y)))
      responses <- .pltest_root_value(responses)
      local <- numeric(length(ids))
      pieces <- .npRmpi_bootstrap_collective_apply(length(ids),length(ids),function(j,pos) {
        result <- vapply(j,function(k) .indextest_local_refits(data,ids[k],1L,
                            responses[,k,drop=FALSE]),numeric(1))
        list(ids=j,values=result)
      })
      for (piece in pieces) {
        local[piece$ids] <- piece$values
        counts[rank+1L] <- counts[rank+1L]+length(piece$ids)
      }
      values[ids] <- mpi.allreduce(local,type=2L,op="sum",comm=1L)
      .pltest_local(if (rank==0L) .np_progress_activity_step(progress, done = tail(ids,1L)))
    }

    complete <- TRUE
    counts <- mpi.allreduce(counts,type=1L,op="sum",comm=1L)
    names(counts) <- seq_len(size)-1L
    return(list(statistics=values,rank.counts=counts[counts>0L]))
  }
  size <- min(8L,max(1L,ceiling(B/(.npRmpi_bootstrap_worker_count()+1L))))
  tasks <- .npRmpi_bootstrap_chunk_tasks(B,size,with.seeds=FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks,function(task) {
    for (i in seq_len(task$bsz)) .indextest_response(data)
  })
  worker <- function(task, data) {
    .npRmpi_bootstrap_use_task_rng(task)
    v <- .indextest_local_refits(data,task$start,task$bsz)
    cbind(statistic=v,rank=mpi.comm.rank(1L))
  }
  out <- .np_spec_test_fanout(B,tasks,worker,ncol.out=2L,what="single-index test",
    progress.label="Bootstrapping the single-index test",profile.where="mpi.applyLB:npindextest",
    required.bindings=list(.indextest_local_refits=.indextest_local_refits),data=data)
  counts <- table(out[,2L])
  list(statistics=as.numeric(out[,1L]),rank.counts=setNames(as.integer(counts),names(counts)))
}
