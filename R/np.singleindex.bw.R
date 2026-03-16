
npindexbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
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

    tbw <- npindexbw(xdat = xdat, ydat = ydat, ...)

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

    xdat <- toFrame(xdat)

    bws <- double(ncol(xdat)+1)

    tbw <- npindexbw.default(xdat = xdat,
                             ydat = ydat,
                             bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
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

.npindex_nn_candidate_bandwidth <- function(h, bwtype, nobs) {
  if (identical(bwtype, "fixed")) {
    return(list(ok = is.finite(h) && (h > 0), value = as.double(h)))
  }

  if (!is.finite(h)) {
    return(list(ok = FALSE, value = NA_real_))
  }

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  k <- .np_round_half_to_even(h)
  list(ok = (k >= lower) && (k <= upper), value = as.double(k))
}

.npindex_default_start_bandwidth <- function(fit, bwtype, nobs) {
  if (identical(bwtype, "fixed")) {
    return(1.059224 * EssDee(fit) * nobs^(-1 / 5))
  }

  lower <- 2L
  max(lower, min(max(1L, as.integer(nobs) - 1L), .np_round_half_to_even(sqrt(nobs))))
}

.npindex_random_start_bandwidth <- function(fit, bwtype, nobs) {
  if (identical(bwtype, "fixed")) {
    return(runif(1, min = 0.5, max = 1.5) * EssDee(fit) * nobs^(-1 / 5))
  }

  upper <- max(1L, as.integer(nobs) - 1L)
  runif(1, min = 2, max = max(2L, upper))
}

.npindex_finalize_bandwidth <- function(h, bwtype, nobs, where = "npindexbw") {
  candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bwtype, nobs = nobs)
  if (!candidate$ok) {
    if (identical(bwtype, "fixed")) {
      stop(sprintf("%s: bandwidth must be positive and finite", where), call. = FALSE)
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

  candidate$value
}

.npindex_lp_loo_fit <- function(index,
                                ydat,
                                h,
                                bws,
                                spec,
                                chunk.size = 128L) {
  index <- as.double(index)
  ydat <- as.double(ydat)
  n <- length(index)
  if (length(ydat) != n)
    stop("index and ydat must have the same length")

  index.df <- data.frame(index = index)
  kbw <- kbandwidth(
    bw = c(h),
    bwtype = bws$type,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder,
    ckerbound = bws$ckerbound,
    ckerlb = bws$ckerlb,
    ckerub = bws$ckerub,
    nobs = n,
    xdati = untangle(index.df),
    ydati = untangle(data.frame(ydat)),
    xnames = "index",
    ynames = bws$ynames
  )
  degree <- as.integer(spec$degree.engine)
  W <- W.lp(
    xdat = index.df,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  W.eval <- W.lp(
    xdat = index.df,
    exdat = index.df,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  if (!is.matrix(W))
    W <- matrix(W, nrow = n)
  if (!is.matrix(W.eval))
    W.eval <- matrix(W.eval, nrow = n)

  ridge.grid <- npRidgeSequenceFromBase(n.train = n, ridge.base = 0.0, cap = 1.0)
  fit.loo <- double(n)
  chunk.size <- max(1L, npValidatePositiveInteger(chunk.size, "chunk.size"))

  for (start in seq.int(1L, n, by = chunk.size)) {
    idx <- seq.int(start, min(n, start + chunk.size - 1L))
    kw <- .np_kernel_weights_direct(
      bws = kbw,
      txdat = index.df,
      exdat = index.df[idx, , drop = FALSE],
      bandwidth.divide = TRUE,
      kernel.pow = 1.0
    )

    kw[cbind(idx, seq_along(idx))] <- 0.0
    W.eval.block <- W.eval[idx, , drop = FALSE]

    for (jj in seq_along(idx)) {
      w <- kw[, jj]
      XtWy0 <- sum(w * ydat)
      lc.mean <- XtWy0 / NZD(sum(w))
      A <- crossprod(W, W * w)
      z <- crossprod(W, ydat * w)
      solved <- FALSE

      for (ridge in ridge.grid) {
        Ar <- A
        zr <- z
        if (ridge > 0) {
          diag(Ar) <- diag(Ar) + ridge
          zr[1L] <- XtWy0 + ridge * XtWy0 / NZD(A[1L, 1L])
        }

        coef <- tryCatch(chol2inv(chol(Ar)) %*% zr, error = function(e) NULL)
        if (!is.null(coef) && all(is.finite(coef))) {
          alpha <- min(1.0, ridge)
          fit.loo[idx[jj]] <- (1.0 - alpha) * drop(W.eval.block[jj, , drop = FALSE] %*% coef) +
            alpha * lc.mean
          solved <- TRUE
          break
        }
      }

      if (!solved)
        fit.loo[idx[jj]] <- lc.mean
    }
  }

  fit.loo
}

npindexbw.default <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws,
           bandwidth.compute = TRUE,
           basis = c("glp", "additive", "tensor"),
           bernstein.basis = FALSE,
           degree = NULL,
           nmulti,
           only.optimize.beta,
           optim.abstol,
           optim.maxattempts,
           optim.maxit,
           optim.method,
           optim.reltol,
           random.seed,
           regtype = c("lc", "ll", "lp"),
           ...){

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
    spec <- npResolveCanonicalConditionalRegSpec(
      mc.names = names(match.call(expand.dots = FALSE)),
      regtype = regtype,
      basis = basis,
      degree = degree,
      bernstein.basis = bernstein.basis,
      ncon = 1L,
      where = "npindexbw"
    )
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

    if (tbw$method == "kleinspady" && !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("nmulti","random.seed", "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit", "only.optimize.beta")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (bandwidth.compute) {
      bwsel.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
      if (any.m) {
        nms <- mc.names[m]
        bwsel.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
      }
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting single-index bandwidth",
        do.call(npindexbw.sibandwidth, bwsel.args)
      )
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

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
           ...){

    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)


    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- min(5,ncol(xdat))
    }
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    only.optimize.beta <- npValidateScalarLogical(only.optimize.beta, "only.optimize.beta")
    nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
    .np_progress_bandwidth_set_total(nmulti)
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")

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

    total.time <-
      system.time({

        if(bandwidth.compute){

          ## Invariant objects used by objective evaluations.
          xmat <- xdat
          wmat <- cbind(ydat, 1.0)
          bandwidth_eval_count <- 0L

          bandwidth_progress_step <- function() {
            bandwidth_eval_count <<- bandwidth_eval_count + 1L
            .np_progress_bandwidth_activity_step(done = bandwidth_eval_count)
            invisible(NULL)
          }

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

            ## Next we define the sum of squared leave-one-out residuals

            sum.squares.leave.one.out <- function(beta,h) {

              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)

              index <- xmat %*% c(1, beta)

              fit.loo <- if (identical(spec$regtype.engine, "lc")) {
                ## One call to npksum to avoid repeated computation of the
                ## product kernel (the expensive part)
                tww <- tryCatch(
                  npksum(txdat=index,
                         tydat=wmat,
                         weights=wmat,
                         leave.one.out=TRUE,
                         bandwidth.divide=TRUE,
                         bws=c(h),
                         bwtype = bws$type,
                         ckertype = bws$ckertype,
                         ckerorder = bws$ckerorder,
                         ckerbound = bws$ckerbound,
                         ckerlb = bws$ckerlb,
                         ckerub = bws$ckerub)$ksum,
                  error = function(e) NULL
                )
              if (is.null(tww))
                  return(rep.int(NA_real_, length(ydat)))
                tww[1,2,]/NZD(tww[2,2,])
              } else {
                ok.design <- tryCatch({
                  npCheckRegressionDesignCondition(
                    reg.code = REGTYPE_LP,
                    xcon = data.frame(index = index),
                    basis = spec$basis.engine,
                    degree = spec$degree.engine,
                    bernstein.basis = spec$bernstein.basis.engine,
                    where = "npindexbw"
                  )
                  TRUE
                }, error = function(e) FALSE)
                if (!ok.design)
                  return(ichimuraMaxPenalty)

                .npindex_lp_loo_fit(
                  index = index,
                  ydat = ydat,
                  h = h,
                  bws = bws,
                  spec = spec
                )
              }

              if (any(!is.finite(fit.loo)))
                return(ichimuraMaxPenalty)

              t.ret <- mean((ydat-fit.loo)^2)
              return(t.ret)

            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              return(sum.squares.leave.one.out(beta,h))
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

          kleinspadyFloor <- sqrt(.Machine$double.eps)

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

            ## Next we define the sum of logs

            sum.log.leave.one.out <- function(beta,h) {

              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)

              index <- xmat %*% c(1, beta)

              ks.loo <- if (identical(spec$regtype.engine, "lc")) {
                ## One call to npksum to avoid repeated computation of the
                ## product kernel (the expensive part)
                tww <- tryCatch(
                  npksum(txdat=index,
                         tydat=wmat,
                         weights=wmat,
                         leave.one.out=TRUE,
                         bandwidth.divide=TRUE,
                         bws=c(h),
                         bwtype = bws$type,
                         ckertype = bws$ckertype,
                         ckerorder = bws$ckerorder,
                         ckerbound = bws$ckerbound,
                         ckerlb = bws$ckerlb,
                         ckerub = bws$ckerub)$ksum,
                  error = function(e) NULL
                )
              if (is.null(tww))
                  return(rep.int(NA_real_, length(ydat)))
                tww[1,2,]/NZD(tww[2,2,])
              } else {
                ok.design <- tryCatch({
                  npCheckRegressionDesignCondition(
                    reg.code = REGTYPE_LP,
                    xcon = data.frame(index = index),
                    basis = spec$basis.engine,
                    degree = spec$degree.engine,
                    bernstein.basis = spec$bernstein.basis.engine,
                    where = "npindexbw"
                  )
                  TRUE
                }, error = function(e) FALSE)
                if (!ok.design)
                  return(sqrt(.Machine$double.xmax))

                .npindex_lp_loo_fit(
                  index = index,
                  ydat = ydat,
                  h = h,
                  bws = bws,
                  spec = spec
                )
              }

              if (any(!is.finite(ks.loo)))
                return(sqrt(.Machine$double.xmax))

              ## Avoid taking log of zero (ks.loo = 0 or 1 since we take
              ## the log of ks.loo and the log of 1-ks.loo)

              ks.loo[ks.loo < kleinspadyFloor] <- kleinspadyFloor
              ks.loo[ks.loo > 1 - kleinspadyFloor] <- 1 - kleinspadyFloor

              ## Maximize the log likelihood, therefore minimize minus.
              ## Here ydat is binary (0/1), so this is equivalent to
              ## ifelse(ydat==1, log(ks.loo), log(1-ks.loo)) but avoids
              ## branchy ifelse and uses stable log1p for the second term.
              t.ret <- -mean(ydat * log(ks.loo) + (1.0 - ydat) * log1p(-ks.loo))
              return(t.ret)

            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              return(sum.log.leave.one.out(beta,h))
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
          num.feval.overall <- 0

          if(bws$method == "ichimura"){
            optim.fn <- if(only.optimize.beta) ichimura.nobw else ichimura
            optim.control <- list(abstol=optim.abstol,
                                  reltol=optim.reltol,
                                  maxit=optim.maxit)
          } else if(bws$method == "kleinspady"){
            optim.fn <- if(only.optimize.beta) kleinspady.nobw else  kleinspady
            optim.control <- list(reltol=optim.reltol,maxit=optim.maxit)
          }

          for (i in seq_len(nmulti)) {
            bandwidth_eval_count <- 0L
            ## We use the nlm command to minimize the objective function using
            ## starting values. Note that since we normalize beta_1=1 here beta
            ## is the k-1 vector containing beta_2...beta_k

            if(i == 1) {

              ## Initial values taken from OLS fit with a constant used for
              ## multistart 1
              ols.fit <- lm(ydat~xdat,x=TRUE)
              fit <- fitted(ols.fit)

              if (p != 1L){
                if (setequal(bws$beta[2:p], c(0)))
                  beta <- coef(ols.fit)[3:ncol(ols.fit$x)]
                else
                  beta = bws$beta[2:p]
              } else { beta = numeric(0) }

              if (bws$bw == 0)
                h <- .npindex_default_start_bandwidth(fit = fit, bwtype = bws$type, nobs = nobs)
              else
                h <- .npindex_finalize_bandwidth(h = bws$bw, bwtype = bws$type, nobs = nobs)
            } else {
              ## Random initialization used for remaining multistarts

              beta.length <- length(coef(ols.fit)[3:ncol(ols.fit$x)])
              beta <- runif(beta.length,min=0.5,max=1.5)*coef(ols.fit)[3:ncol(ols.fit$x)]
              if (!only.optimize.beta)
                h <- .npindex_random_start_bandwidth(fit = fit, bwtype = bws$type, nobs = nobs)
            }

            optim.parm <- if(only.optimize.beta) beta else c(beta,h)

            optim.base.args <- list(
              par = optim.parm,
              fn = optim.fn,
              gr = NULL,
              method = optim.method,
              control = optim.control
            )
            if (only.optimize.beta) {
              optim.base.args$h <- h
            }
            suppressWarnings(optim.return <- do.call(optim, optim.base.args))
            if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
              num.feval.overall <- num.feval.overall + optim.return$counts[1]
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              beta.length <- length(coef(ols.fit)[3:ncol(ols.fit$x)])
              beta <- runif(beta.length,min=0.5,max=1.5)*coef(ols.fit)[3:ncol(ols.fit$x)]
              if(!only.optimize.beta)
                h <- .npindex_random_start_bandwidth(fit = fit, bwtype = bws$type, nobs = nobs)

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
              
              optim.base.args$par <- optim.parm
              optim.base.args$control <- optim.control
              if (!only.optimize.beta && ("h" %in% names(optim.base.args))) {
                optim.base.args$h <- NULL
              }
              if (only.optimize.beta) {
                optim.base.args$h <- h
              }
              suppressWarnings(optim.return <- do.call(optim, optim.base.args))
              if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                num.feval.overall <- num.feval.overall + optim.return$counts[1]
            }

            if(optim.return$convergence != 0)
              stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))

            fval.value[i] <- optim.return$value
            if(optim.return$value < fval.min) {
              param <- if(only.optimize.beta) c(optim.return$par,h) else optim.return$par
              fval.min <- optim.return$value
              numimp <- numimp + 1
              best <- i
            }

            .np_progress_bandwidth_multistart_step(done = i, total = nmulti)

          }

          bws$beta <- c(1.0, param[beta.idx])
          bws$bw <- .npindex_finalize_bandwidth(h = param[p], bwtype = bws$type, nobs = nobs)
          bws$fval <- fval.min
          bws$ifval <- best
          bws$num.feval <- num.feval.overall
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

    .np_seed_exit(seed.state)

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

    bws

  }
