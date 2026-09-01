# This function implements an individual test of significance for both
# discrete (Racine, hart, Li, 2006, ER) and continuous variables
# (Racine, 1997, JBES). It accepts a data frame for explanatory data
# (mixed datatypes allowed), a vector for y for a regression model, an
# npregbw object, and a set of indices for the columns of X for which
# the test is to be run (default = all).

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

.np_npsig_bootstrap_bw_reselect <- function(xdat,
                                            ydat,
                                            bws.seed,
                                            extra.args = list(),
                                            bootstrap.iter,
                                            bw.fun = npregbw) {
  bw.args <- if (length(extra.args)) extra.args else list()
  bw.args[c("xdat", "ydat", "bws")] <- NULL

  user.nmulti <- !is.null(names(bw.args)) &&
    "nmulti" %in% names(bw.args) &&
    !is.null(bw.args$nmulti)

  if (!user.nmulti && bootstrap.iter > 1L)
    bw.args$nmulti <- 1L

  .np_progress_with_legacy_suppressed(
    do.call(bw.fun, c(list(xdat = xdat, ydat = ydat, bws = bws.seed), bw.args))
  )
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

  if (is.factor(ydat))
    stop("dependent variable must be continuous.")

  ## Test for valid entries in index before RNG or estimator work.

  if(anyNA(index)) stop("index must not contain missing values")
  if(any(index < 1 | index > NCOL(xdat), na.rm = TRUE)) stop(paste("invalid index provided: index entries must lie between 1 and ",NCOL(xdat),sep=""))
  if(length(unique(index)) < length(index)) stop("index contains repeated values (must be unique)")

  pivot.plan <- .np_npsig_pivot_plan(
    pivot = pivot,
    xdat = xdat,
    index = index,
    joint = joint
  )
  pivot <- pivot.plan$requested
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
  on.exit({
    if (isTRUE(progress.active))
      .np_progress_abort(progress)
  }, add = TRUE)

  ## Save seed prior to setting.  The computational owner is selected before
  ## this boundary so unsupported calls enter the incumbent without RNG work.

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  if(boot.type=="II") {
    ## Store a copy of the bandwidths passed in
    bws.original <- bws
  }

  num.obs <- nrow(xdat)

  if(!joint) {

    In <- numeric(length(index))
    P <- numeric(length(index))

  }

  ## Some constants for the wild bootstrap

  a <- -0.6180339887499  # (1-sqrt(5))/2
  b <- 1.6180339887499   # (1+sqrt(5))/2
  P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))

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

    if(boot.type=="II") {

      ## Reset bw vals to original as the ith component of bws gets
      ## overwritten when index changes so needs to be set to its
      ## original value

      bws <- bws.original

    }

    ## Note - xdat must be a data frame

    ## Construct In, the average value of the squared derivatives of
    ## the jth element, discrete or continuous

    progress <- .np_progress_step(progress)
    npreg.out <- npreg(txdat = xdat,
                       tydat = ydat,
                       bws = bws,
                       gradients = TRUE,
                       se = pivot.use,
                       ...)

    In <- .np_npsig_statistic(npreg.out, index = index, pivot = pivot.use)
    progress <- .np_progress_step(progress)

    if(boot.method != "pairwise") {

      ## Compute scale and mean of unrestricted residuals

      progress <- .np_progress_step(progress)
      ei.unres <- scale(residuals(npreg(bws=bws)))
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
      mhat.xi <-  npreg(txdat = xdat,
                        tydat = ydat,
                        exdat = xdat.eval,
                        bws = bws,
                        ...)$mean
      progress <- .np_progress_step(progress)

      ## Rescale and recenter the residuals under the null to those
      ## under the alternative
      
      ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
      
      ## Recenter the residuals to have mean zero

      ei <- ei - mean(ei)
      
    }
    
    if(boot.type=="II")
      bws.boot.prev <- bws.original

    if (streamed.iid) {
      tile.width <- 8L
      for (tile.start in seq.int(1L, B, by = tile.width)) {
        tile.count <- min(tile.width, B - tile.start + 1L)
        if (identical(boot.method, "iid")) {
          donor.index <- vapply(
            seq_len(tile.count),
            function(unused) sample.int(num.obs, replace = TRUE),
            integer(num.obs)
          )
          response.matrix <- matrix(
            mhat.xi + ei[donor.index],
            nrow = num.obs,
            ncol = tile.count
          )
        } else {
          wild.values <- if (identical(boot.method, "wild"))
            c(a, b, P.a) else c(-1, 1, P.a)
          response.matrix <- vapply(
            seq_len(tile.count),
            function(unused)
              mhat.xi + ei * draw.wild.mult(
                num.obs, wild.values[[1L]], wild.values[[2L]],
                wild.values[[3L]]
              ),
            numeric(num.obs)
          )
        }
        tile.statistic <- .np_npsig_streamed_response_statistic(
          bws = bws,
          xdat = xdat,
          index = index,
          response.matrix = response.matrix,
          pivotal = pivot.use
        )
        tile.rows <- tile.start:(tile.start + tile.count - 1L)
        In.vec[tile.rows] <- tile.statistic
        tile.done <- tile.rows[[length(tile.rows)]]
        if (!isTRUE(progress$known_total)) {
          progress <- .np_npsig_progress_promote(
            progress, total = B, done = tile.start
          )
          progress <- .np_progress_step(progress, done = tile.done)
        } else {
          progress <- .np_progress_step(progress, done = tile.done)
        }
      }
    } else for (i.star in seq_len(B)) {
      if(boot.method == "iid") {

        ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]

      } else if(boot.method == "wild") {

        ## Conduct a wild bootstrap. We generate a sample for ydat
        ## (ydat.star) drawn from the conditional mean evaluated
        ## holding the variable tested at its median, and add to that
        ## a wild bootstrap draw from the original disturbance vector

        ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)

      } else if(boot.method == "wild-rademacher") {

        ## Conduct a wild bootstrap. We generate a sample for ydat
        ## (ydat.star) drawn from the conditional mean evaluated
        ## holding the variable tested at its median, and add to that
        ## a wild bootstrap draw from the original disturbance vector

        ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, -1, 1, P.a)

      } else if(boot.method =="pairwise") {

        ## Leave variable being tested untouched, resample remaining
        ## pairs of y,X thereby breaking any systematic relationship
        ## between variable being tested in y
        boot.index <- sample.int(num.obs, replace = TRUE)
        ydat.star <- ydat[boot.index]
        xdat.star <- xdat[boot.index,]
        for(i in index) xdat.star[,i] <- xdat[,i]

      }

      if(boot.type=="II") {

        ## Bootstrap II reuses the previous bootstrap optimum as a
        ## hot start and drops to nmulti=1 after the first
        ## re-selection unless the user explicitly supplied nmulti.

        if(boot.method == "pairwise") {

          bws.boot <- .np_npsig_bootstrap_bw_reselect(
            xdat = xdat.star,
            ydat = ydat.star,
            bws.seed = bws.boot.prev,
            extra.args = extra.args,
            bootstrap.iter = i.star
          )

        } else {

          bws.boot <- .np_npsig_bootstrap_bw_reselect(
            xdat = xdat,
            ydat = ydat.star,
            bws.seed = bws.boot.prev,
            extra.args = extra.args,
            bootstrap.iter = i.star
          )

        }

        bws.boot.prev <- bws.boot

        ## Copy the new cross-validated bandwidth for variable i into
        ## bw.original and use this below.

        bws <- bws.original

        bws$bw[index] <- bws.boot$bw[index]

      }

      if (direct.pairwise) {

        In.vec[i.star] <- .np_npsig_streamed_response_statistic(
          bws = bws,
          xdat = xdat.star,
          index = index,
          response.matrix = matrix(ydat.star, ncol = 1L),
          pivotal = pivot.use
        )

      } else if(boot.method == "pairwise") {

        npreg.boot <- npreg(txdat = xdat.star,
                            tydat = ydat.star,
                            bws = bws,
                            gradients = TRUE,
                            se = pivot.use,
                            ...)

      } else {

        npreg.boot <- npreg(txdat = xdat,
                            tydat = ydat.star,
                            bws = bws,
                            gradients = TRUE,
                            se = pivot.use,
                            ...)

      }

      if (!direct.pairwise)
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
      streamed.unrestricted <- npreg(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = TRUE,
        se = any(pivot.plan$effective),
        ...
      )
      progress <- .np_progress_step(progress)
      progress <- .np_progress_step(progress)
      streamed.ei.unres <- scale(residuals(npreg(bws = bws)))
      streamed.ei.unres.scale <- attr(streamed.ei.unres, "scaled:scale")
      streamed.ei.unres.center <- attr(streamed.ei.unres, "scaled:center")
      streamed.ei.unres <- NULL
      progress <- .np_progress_step(progress)
    }

    for(i in index) {
      
      ## Increment counter...
      
      ii <- ii + 1
      progress$label <- paste("Testing", tested.names[[ii]])
      pivot.use <- pivot.plan$effective[[ii]]
      
      if(boot.type=="II") {
        
        ## Reset bw vals to original as the ith component of bws gets
        ## overwritten when index changes so needs to be set to its
        ## original value
        
        bws <- bws.original
        
      }
      
      ## Note - xdat must be a data frame
      
      ## Construct In, the average value of the squared derivatives of
      ## the jth element, discrete or continuous
      
      if (streamed.iid) {
        npreg.out <- streamed.unrestricted
      } else {
        progress <- .np_progress_step(progress)
        npreg.out <- npreg(txdat = xdat,
                           tydat = ydat,
                           bws = bws,
                           gradients = TRUE,
                           se = pivot.use,
                           ...)
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
          ei.unres <- scale(residuals(npreg(bws=bws)))
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
        mhat.xi <-  npreg(txdat = xdat,
                          tydat = ydat,
                          exdat = xdat.eval,
                          bws = bws,
                          ...)$mean
        progress <- .np_progress_step(progress)
        
        ## Rescale and recenter the residuals under the null to those
        ## under the alternative
        
        ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
        
        ## Recenter the residuals to have mean zero
        
        ei <- ei - mean(ei)
        
      }
      
      if(boot.type=="II")
        bws.boot.prev <- bws.original

      if (streamed.iid) {
        tile.width <- 8L
        for (tile.start in seq.int(1L, B, by = tile.width)) {
          tile.count <- min(tile.width, B - tile.start + 1L)
          donor.index <- NULL
          response.matrix <- NULL
          if (identical(boot.method, "iid")) {
            donor.index <- vapply(
              seq_len(tile.count),
              function(unused) sample.int(num.obs, replace = TRUE),
              integer(num.obs)
            )
          } else {
            wild.values <- if (identical(boot.method, "wild"))
              c(a, b, P.a) else c(-1, 1, P.a)
            response.matrix <- vapply(
              seq_len(tile.count),
              function(unused)
                mhat.xi + ei * draw.wild.mult(
                  num.obs, wild.values[[1L]], wild.values[[2L]],
                  wild.values[[3L]]
                ),
              numeric(num.obs)
            )
          }
          tile.statistic <- .np_npsig_streamed_iid_tile(
            bws = bws,
            xdat = xdat,
            tested.index = i,
            donor.index = donor.index,
            response.matrix = response.matrix,
            null.mean = mhat.xi,
            residual.pool = ei
          )
          tile.rows <- tile.start:(tile.start + tile.count - 1L)
          In.vec[tile.rows] <- tile.statistic
          tile.done <- tile.rows[[length(tile.rows)]]
          if (length(index) == 1L && !isTRUE(progress$known_total)) {
            progress <- .np_npsig_progress_promote(
              progress, total = B, done = tile.done
            )
          } else if (length(index) == 1L) {
            progress <- .np_progress_step(progress, done = tile.done)
          } else {
            progress <- .np_progress_step(progress)
          }
        }
      } else for (i.star in seq_len(B)) {
        if(boot.method == "iid") {
          
          ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
          
        } else if(boot.method == "wild") {
          
          ## Conduct a wild bootstrap. We generate a sample for ydat
          ## (ydat.star) drawn from the conditional mean evaluated
          ## holding the variable tested at its median, and add to that
          ## a wild bootstrap draw from the original disturbance vector
          
          ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
          
        } else if(boot.method == "wild-rademacher") {
          
          ## Conduct a wild bootstrap. We generate a sample for ydat
          ## (ydat.star) drawn from the conditional mean evaluated
          ## holding the variable tested at its median, and add to that
          ## a wild bootstrap draw from the original disturbance vector
          
          ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, -1, 1, P.a)
          
        } else if(boot.method =="pairwise") {
          
          ## Leave variable being tested untouched, resample remaining
          ## pairs of y,X thereby breaking any systematic relationship
          ## between variable being tested in y
          boot.index <- sample.int(num.obs, replace = TRUE)
          ydat.star <- ydat[boot.index]
          xdat.star <- xdat
          xdat.star[,-i] <- xdat[boot.index,-i]
          
        }
        
        if(boot.type=="II") {
          
          ## Bootstrap II reuses the previous bootstrap optimum as a
          ## hot start and drops to nmulti=1 after the first
          ## re-selection unless the user explicitly supplied nmulti.
          
          if(boot.method == "pairwise") {
            
            bws.boot <- .np_npsig_bootstrap_bw_reselect(
              xdat = xdat.star,
              ydat = ydat.star,
              bws.seed = bws.boot.prev,
              extra.args = extra.args,
              bootstrap.iter = i.star
            )
            
          } else {
            
            bws.boot <- .np_npsig_bootstrap_bw_reselect(
              xdat = xdat,
              ydat = ydat.star,
              bws.seed = bws.boot.prev,
              extra.args = extra.args,
              bootstrap.iter = i.star
            )
            
          }

          bws.boot.prev <- bws.boot
          
          ## Copy the new cross-validated bandwidth for variable i into
          ## bw.original and use this below.
          
          bws <- bws.original
          
          bws$bw[i] <- bws.boot$bw[i]
          
        }
        
        if (direct.pairwise) {

          In.vec[i.star] <- .np_npsig_streamed_response_statistic(
            bws = bws,
            xdat = xdat.star,
            index = i,
            response.matrix = matrix(ydat.star, ncol = 1L),
            pivotal = pivot.use
          )

        } else if(boot.method == "pairwise") {
          
          npreg.boot <- npreg(txdat = xdat.star,
                              tydat = ydat.star,
                              bws = bws,
                              gradients = TRUE,
                              se = pivot.use,
                              ...)
          
        } else {
          
          npreg.boot <- npreg(txdat = xdat,
                              tydat = ydat.star,
                              bws = bws,
                              gradients = TRUE,
                              se = pivot.use,
                              ...)
          
        }
        
        if (!direct.pairwise)
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

npsigtest.default <- function(bws, xdat, ydat, ...){
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

  if(xdat.named)
    xdat <- toFrame(xdat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npregbw)

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
