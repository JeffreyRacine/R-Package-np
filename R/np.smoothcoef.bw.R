npscoefbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat", "zdat"),
                                     eval_env = parent.frame())
    UseMethod("npscoefbw", target)
  }

npscoefbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = terms(formula),
                        data = environment(formula),
                        eval_env = environment(formula))
    else .np_terms_ts_mask(terms_obj = terms(formula, data = data),
                           data = data,
                           eval_env = environment(formula))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    formula.call <- .np_bw_formula_from_call(call_obj = call, eval_env = parent.frame())
    if (!is.null(formula.call))
      mf[[2]] <- formula.call

    mf[[1]] <- as.name("model.frame")

    formula.obj <- .np_bw_resolve_formula(formula_obj = formula,
                                        formula_call = formula.call,
                                        eval_env = parent.frame())
    chromoly <- explodePipe(formula.obj, env = environment(formula))

    bronze <- sapply(chromoly, paste, collapse = " + ")
    mf[["formula"]] <-
      as.formula(paste(bronze[1]," ~ ",
                       paste(bronze[2:length(bronze)],
                             collapse =" + ")),
                 env = environment(formula))

    mf[["formula"]] <- terms(mf[["formula"]])
    if(all(orig.ts)){
      args <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      attr(mf[["formula"]], "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
    }
    
    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())
    
    ydat <- model.response(mf)
    xdat <- mf[, chromoly[[2]], drop = FALSE]
    miss.z <- !(length(chromoly) == 3)
    if (!miss.z)
      zdat <- mf[, chromoly[[3]], drop = FALSE]
    
    bw.args <- list(xdat = xdat, ydat = ydat)
    if (!miss.z)
      bw.args$zdat <- zdat
    tbw <- do.call(npscoefbw, c(bw.args, list(...)))

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$chromoly <- chromoly

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)
    
    tbw
  }

npscoefbw.NULL <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    miss.z <- missing(zdat)

    xdat <- toFrame(xdat)

    if(!miss.z)
      zdat <- toFrame(zdat)

    n.bw <- if (miss.z) ncol(xdat) else ncol(zdat)
    bws <- double(n.bw)

    bw.args <- list(xdat = xdat, ydat = ydat, bws = bws)
    if (!miss.z)
      bw.args$zdat <- zdat
    tbw <- do.call(npscoefbw.default, c(bw.args, list(...)))

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <-
      updateBwNameMetadata(nameList = list(ynames = deparse(substitute(ydat))),
                           bws = tbw)

    tbw
  }

.npscoef_nn_candidate_bandwidth <- function(param, bwtype, nobs) {
  if (identical(bwtype, "fixed"))
    return(as.double(param))

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  vapply(param, function(h) {
    if (!is.finite(h))
      return(NA_real_)
    as.double(max(lower, min(upper, .np_round_half_to_even(h))))
  }, numeric(1))
}

.npscoef_default_start_bandwidth <- function(param, bwtype, nobs) {
  if (identical(bwtype, "fixed"))
    return(as.double(param))

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  start <- max(lower, min(upper, .np_round_half_to_even(sqrt(nobs))))
  rep.int(as.double(start), length(param))
}

.npscoef_random_start_bandwidth <- function(param, bwtype, nobs) {
  if (identical(bwtype, "fixed"))
    return(as.double(runif(length(param), min = 0.5, max = 1.5) * param))

  upper <- max(1L, as.integer(nobs) - 1L)
  .npscoef_nn_candidate_bandwidth(
    param = runif(length(param), min = 2, max = max(2L, upper)),
    bwtype = bwtype,
    nobs = nobs
  )
}

.npscoef_candidate_is_admissible <- function(param, bwtype, nobs) {
  candidate <- .npscoef_nn_candidate_bandwidth(
    param = param,
    bwtype = bwtype,
    nobs = nobs
  )
  if (any(!is.finite(candidate)))
    return(FALSE)
  if (identical(bwtype, "fixed"))
    return(all(candidate > 0))
  TRUE
}

.npscoef_finalize_bandwidth <- function(param, bwtype, nobs, where = "npscoefbw") {
  candidate <- .npscoef_nn_candidate_bandwidth(param = param, bwtype = bwtype, nobs = nobs)
  if (any(!is.finite(candidate))) {
    if (identical(bwtype, "fixed")) {
      stop(sprintf("%s: bandwidth must be finite", where), call. = FALSE)
    }
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth must be an integer vector in [2, %d]",
        where,
        max(2L, as.integer(nobs) - 1L)
      ),
      call. = FALSE
    )
  }
  if (identical(bwtype, "fixed") && any(candidate <= 0)) {
    stop(sprintf("%s: bandwidth must be strictly positive", where), call. = FALSE)
  }
  as.double(candidate)
}

npscoefbw.scbandwidth <- 
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws,
           backfit.iterate = FALSE,
           backfit.maxiter = 100,
           backfit.tol = .Machine$double.eps,
           bandwidth.compute = TRUE,
           cv.iterate = FALSE,
           cv.num.iterations = 1,
           nmulti,
           optim.abstol = .Machine$double.eps,
           optim.maxattempts = 10,
           optim.maxit = 500,
           optim.method = c("Nelder-Mead", "BFGS", "CG"),
           optim.reltol = sqrt(.Machine$double.eps),
           random.seed = 42,
           ...){
    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)


    miss.z <- missing(zdat)
    
    xdat <- toFrame(xdat)

    if (!miss.z)
      zdat <- toFrame(zdat)
    
    if (missing(nmulti)){
      nmulti <- min(5,length(bws$bw))
    }
    regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
    cv.iterate <- npValidateScalarLogical(cv.iterate, "cv.iterate")
    backfit.iterate <- npValidateScalarLogical(backfit.iterate, "backfit.iterate")
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
    .np_progress_bandwidth_set_total(nmulti)
    backfit.maxiter <- npValidatePositiveInteger(backfit.maxiter, "backfit.maxiter")
    backfit.tol <- npValidatePositiveFiniteNumeric(backfit.tol, "backfit.tol")
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")
    if (cv.iterate)
      cv.num.iterations <- npValidatePositiveInteger(cv.num.iterations, "cv.num.iterations")
    spec <- .npscoef_canonical_spec(
      source = bws,
      zdat = if (miss.z) xdat else zdat,
      where = "npscoefbw"
    )
    reg.engine <- spec$regtype.engine
    if (!identical(reg.engine, "lc") && cv.iterate)
      stop("cv.iterate currently supports regtype='lc' for npscoefbw")
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector or a factor")

    if (miss.z) {
      bwMatch(xdat, bws$xdati)
    } else {
      bwMatch(zdat, bws$zdati)
    }
    
    if (dim(xdat)[1] != length(ydat))
      stop("number of regression data and response data do not match")

    if (ncol(xdat) == 1 && missing(cv.iterate))
      cv.iterate = FALSE

    if (!all(bws$xdati$icon))
      stop("Only continuous 'x' regressors are supported in this version.")

    optim.method <- match.arg(optim.method)
      
    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    train.df <- data.frame(xdat, ydat)
    if (!miss.z)
      train.df <- data.frame(train.df, zdat)
    rows.omit <- attr(na.omit(train.df), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]

    if(!miss.z)
      zdat <- zdat[keep.rows,, drop = FALSE]
    
    nrow = dim(xdat)[1]
    ncol = dim(xdat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat <- toMatrix(xdat)

    ## if (!miss.z)
    ##  zdat <- toMatrix(zdat)
    
    ## bad data
    if (qr(xdat)$rank < ncol(xdat)){
      stop("columns of the independent variable (xdat) are linearly dependent") 
    }

    n <- nrow(xdat)
    
    ## ... do bandwidth selection
    
    ## construct 'W' matrix
    ## in the future one will be able to use a switch to npksum
    ## to emulate W

    W <- cbind(1.0, xdat)
    yW <- cbind(ydat, W)
    
    if (miss.z){
      zdat <- xdat
      dati <- bws$xdati
    }
    else
      dati <- bws$zdati
    zdat.df <- if (is.data.frame(zdat)) zdat else as.data.frame(zdat)

    mysd <- EssDee(zdat[, dati$icon, drop = FALSE])
    nconfac <- n^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- n^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    bws$sdev <- mysd
    bws$nconfac <- nconfac
    bws$ncatfac <- ncatfac
    bw.scale.multiplier <- NULL
    if (bws$scaling) {
      bw.scale.multiplier <- rep(ncatfac, bws$ndim)
      if (any(bws$icon)) {
        icon.cumsum <- cumsum(dati$icon)
        bw.scale.multiplier[bws$icon] <- nconfac * bws$sdev[icon.cumsum[bws$icon]]
      }
    }
    apply_bw_to_scbw <- function(scbw, param) {
      param <- .npscoef_nn_candidate_bandwidth(param = param, bwtype = scbw$type, nobs = n)
      scbw$bw <- param
      if (scbw$scaling)
        scbw$bandwidth[[1]] <- scbw$bw * bw.scale.multiplier
      else
        scbw$bandwidth[[1]] <- scbw$bw
      scbw
    }

    fast_largeh_tol <- getOption("np.largeh.rel.tol", 1e-3)
    if (!is.numeric(fast_largeh_tol) || length(fast_largeh_tol) != 1L ||
        is.na(fast_largeh_tol) || !is.finite(fast_largeh_tol) ||
        fast_largeh_tol <= 0 || fast_largeh_tol >= 0.1)
      fast_largeh_tol <- 1e-3

    fast_disc_tol <- getOption("np.disc.upper.rel.tol", 1e-2)
    if (!is.numeric(fast_disc_tol) || length(fast_disc_tol) != 1L ||
        is.na(fast_disc_tol) || !is.finite(fast_disc_tol) ||
        fast_disc_tol <= 0 || fast_disc_tol >= 0.5)
      fast_disc_tol <- 1e-2

    cont_utol <- switch(
      bws$ckertype,
      gaussian = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
      "truncated gaussian" = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
      epanechnikov = sqrt(fast_largeh_tol),
      uniform = 1.0 - 32.0 * .Machine$double.eps,
      0.0
    )

    cont_hmin <- numeric(0)
    if (any(dati$icon) && is.finite(cont_utol) && cont_utol > 0) {
      zcon <- zdat.df[, dati$icon, drop = FALSE]
      cont_hmin <- vapply(zcon, function(col) {
        vals <- as.double(col)
        vals <- vals[is.finite(vals)]
        if (!length(vals))
          return(Inf)
        diff(range(vals)) / cont_utol
      }, numeric(1))
    }

    disc_upper_tol <- function(upper) {
      max(fast_disc_tol * abs(upper),
          16.0 * .Machine$double.eps * max(1.0, abs(upper)))
    }

    uno_upper <- numeric(0)
    if (any(dati$iuno)) {
      uno_idx <- which(dati$iuno)
      uno_upper <- vapply(uno_idx, function(i) {
        uMaxL(dati$all.nlev[[i]], kertype = bws$ukertype)
      }, numeric(1))
    }

    ord_upper <- numeric(0)
    if (any(dati$iord)) {
      ord_idx <- which(dati$iord)
      ord_upper <- vapply(ord_idx, function(i) {
        oMaxL(dati$all.nlev[[i]], kertype = bws$okertype)
      }, numeric(1))
    }

    npscoef_fast_eligible <- function(sbw) {
      if (!identical(sbw$type, "fixed"))
        return(FALSE)

      bwv <- sbw$bandwidth[[1L]]
      if (!length(bwv) || length(bwv) != length(dati$icon))
        return(FALSE)

      if (any(dati$icon)) {
        bw_cont <- bwv[dati$icon]
        if (any(!is.finite(bw_cont)) || any(bw_cont <= 0) ||
            any(bw_cont < cont_hmin))
          return(FALSE)
      }

      if (any(dati$iuno)) {
        bw_uno <- bwv[dati$iuno]
        ok_uno <- mapply(function(bw, upper) {
          is.finite(bw) && abs(bw - upper) <= disc_upper_tol(upper)
        }, bw = bw_uno, upper = uno_upper, SIMPLIFY = TRUE, USE.NAMES = FALSE)
        if (!all(ok_uno))
          return(FALSE)
      }

      if (any(dati$iord)) {
        bw_ord <- bwv[dati$iord]
        ok_ord <- mapply(function(bw, upper) {
          is.finite(bw) && abs(bw - upper) <= disc_upper_tol(upper)
        }, bw = bw_ord, upper = ord_upper, SIMPLIFY = TRUE, USE.NAMES = FALSE)
        if (!all(ok_ord))
          return(FALSE)
      }

      TRUE
    }

    solve_cv_moment_system <- function(tyw, tww, W.eval.design, maxPenalty, Wz.eval = NULL) {
      neval <- ncol(tyw)
      ncoef <- nrow(tyw)
      pcoef <- ncol(W.eval.design)
      coef.out <- matrix(maxPenalty, nrow = pcoef, ncol = neval)
      ridge.grid <- npRidgeSequenceAdditive(n.train = n, cap = 1.0)
      ridge <- rep.int(ridge.grid[1L], neval)
      ridge.idx <- rep.int(1L, neval)
      doridge <- rep.int(TRUE, neval)

      while(any(doridge)){
        iloo <- seq_len(neval)[doridge]
        for (ii in iloo) {
          doridge[ii] <- FALSE
          ridge.val <- ridge[ii]*tyw[,ii][1]/NZD(tww[,,ii][1,1])
          theta.ii <- tryCatch(
            solve(tww[,,ii] + diag(rep(ridge[ii], ncoef)),
                  tyw[,ii] + c(ridge.val, rep(0, ncoef - 1))),
            error = function(e) e
          )
          if (inherits(theta.ii, "error")) {
            ridge.idx[ii] <- ridge.idx[ii] + 1L
            if (ridge.idx[ii] <= length(ridge.grid)) {
              ridge[ii] <- ridge.grid[ridge.idx[ii]]
              doridge[ii] <- TRUE
            }
            theta.ii <- rep(maxPenalty, ncoef)
          }

          if (is.null(Wz.eval)) {
            coef.out[,ii] <- theta.ii
          } else {
            coef.out[,ii] <- as.vector(crossprod(
              Wz.eval[ii,],
              matrix(theta.ii, nrow = ncol(Wz.eval), ncol = pcoef)
            ))
          }
        }
      }

      coef.out
    }

    lp_full_coef <- function(sbw, leave.one.out.eval) {
      lp_state <- .npscoef_lp_state(
        bws = sbw,
        tzdat = zdat.df,
        ezdat = zdat.df,
        leave.one.out = leave.one.out.eval,
        where = "npscoefbw"
      )
      tensor.train <- .npscoef_row_tensor_design(W, lp_state$W.train)
      ytensor <- cbind(ydat, tensor.train)
      ksum.args <- list(
        txdat = lp_state$z.train,
        tydat = ytensor,
        weights = ytensor,
        bws = lp_state$rbw,
        leave.one.out = leave.one.out.eval,
        bandwidth.divide = TRUE
      )
      main.ks <- do.call(npksum, ksum.args)$ksum
      tyw <- main.ks[-1L, 1L, , drop = FALSE]
      if (length(dim(tyw)) == 3L)
        dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
      tww <- main.ks[-1L, -1L, , drop = FALSE]
      solve_cv_moment_system(
        tyw = tyw,
        tww = tww,
        W.eval.design = W,
        Wz.eval = lp_state$W.eval,
        maxPenalty = maxPenalty
      )
    }

    lp_partial_coef <- function(sbw, wj, partial.y, leave.one.out.eval) {
      lp_state <- .npscoef_lp_state(
        bws = sbw,
        tzdat = zdat.df,
        ezdat = zdat.df,
        leave.one.out = leave.one.out.eval,
        where = "npscoefbw"
      )
      U <- lp_state$W.train * wj
      yU <- cbind(partial.y, U)
      ksum.args <- list(
        txdat = lp_state$z.train,
        tydat = yU,
        weights = yU,
        bws = lp_state$rbw,
        leave.one.out = leave.one.out.eval,
        bandwidth.divide = TRUE
      )
      main.ks <- do.call(npksum, ksum.args)$ksum
      tyw <- main.ks[-1L, 1L, , drop = FALSE]
      if (length(dim(tyw)) == 3L)
        dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
      tww <- main.ks[-1L, -1L, , drop = FALSE]
      as.vector(solve_cv_moment_system(
        tyw = tyw,
        tww = tww,
        W.eval.design = matrix(1.0, nrow = n, ncol = 1L),
        Wz.eval = lp_state$W.eval,
        maxPenalty = maxPenalty
      ))
    }

    total.time <-
      system.time({
        if (bandwidth.compute){
          maxPenalty <- sqrt(.Machine$double.xmax)
          cv_state <- new.env(parent = emptyenv())
          cv_state$fast_total <- 0L
          cv_state$objective_fast <- FALSE
          cv_state$optim_progress <- NULL
          cv_state$optim_eval <- 0L
          cv_state$multistart_index <- NA_integer_
          cv_state$partial_progress <- NULL
          cv_state$partial_eval <- 0L
          cv_state$backfit_iteration <- NA_integer_
          cv_state$partial_index <- NA_integer_

          cv_progress_detail <- function(ridging = FALSE) {
            detail <- sprintf("multistart %d", cv_state$multistart_index)
            if (isTRUE(ridging)) {
              paste(detail, "near-singular system encountered, ridging", sep = ", ")
            } else {
              detail
            }
          }

          cv_progress_begin <- function() {
            cv_state$optim_eval <- 0L
            cv_state$optim_progress <- .np_progress_begin("Optimizing smooth coefficient bandwidth")
            invisible(NULL)
          }

          cv_progress_step <- function(ridging = FALSE) {
            cv_state$optim_eval <- cv_state$optim_eval + 1L
            .np_progress_bandwidth_activity_step(done = cv_state$optim_eval)
            cv_state$optim_progress <- .np_progress_step(
              state = cv_state$optim_progress,
              done = cv_state$optim_eval,
              detail = cv_progress_detail(ridging = ridging)
            )
            invisible(NULL)
          }

          cv_progress_end <- function(state) {
            if (is.null(state))
              return(invisible(NULL))

            if (isTRUE(state$known_total) && identical(state$last_done, state$total))
              return(invisible(NULL))

            state$last_emit <- -Inf
            .np_progress_end(state)
            invisible(NULL)
          }

          cv_progress_finish <- function(ridging = FALSE) {
            if (is.null(cv_state$optim_progress))
              return(invisible(NULL))

            cv_state$optim_progress$last_emit <- -Inf
            cv_state$optim_progress <- .np_progress_end(
              cv_state$optim_progress,
              detail = cv_progress_detail(ridging = ridging)
            )
            cv_state$optim_progress <- NULL
            invisible(NULL)
          }

          partial_progress_detail <- function(fv = NULL) {
            detail <- sprintf(
              "backfitting iteration %d of %d, partial residual %d of %d",
              cv_state$backfit_iteration,
              cv.num.iterations,
              cv_state$partial_index,
              ncol(W)
            )
            if (!is.null(fv)) {
              detail <- paste(
                detail,
                sprintf(
                  "fval %s",
                  format(signif(fv, digits = getOption("digits", 7L)), trim = TRUE)
                ),
                sep = ", "
              )
            }
            detail
          }

          partial_progress_begin <- function(iteration, partial.index) {
            cv_state$backfit_iteration <- iteration
            cv_state$partial_index <- partial.index
            cv_state$partial_eval <- 0L
            cv_state$partial_progress <- .np_progress_begin("Optimizing partial residual bandwidth")
            invisible(NULL)
          }

          partial_progress_step <- function(fv) {
            cv_state$partial_eval <- cv_state$partial_eval + 1L
            cv_state$partial_progress <- .np_progress_step(
              state = cv_state$partial_progress,
              done = cv_state$partial_eval,
              detail = partial_progress_detail(fv = fv)
            )
            invisible(NULL)
          }

          partial_progress_finish <- function(fv = NULL) {
            if (is.null(cv_state$partial_progress))
              return(invisible(NULL))

            cv_state$partial_progress$last_emit <- -Inf
            cv_state$partial_progress <- .np_progress_end(
              cv_state$partial_progress,
              detail = partial_progress_detail(fv = fv)
            )
            cv_state$partial_progress <- NULL
            invisible(NULL)
          }

          overall.cv.ls <- function(param) {
            cv_state$objective_fast <- FALSE
            sbw <- apply_bw_to_scbw(bws, param)
            if (!validateBandwidthTF(sbw) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            cv_state$objective_fast <- npscoef_fast_eligible(sbw)

            if (identical(reg.engine, "lc")) {
              tww <- npksum(txdat = zdat, tydat = yW, weights = yW, bws = sbw,
                            leave.one.out = TRUE)$ksum

              mean.loo <- rep(maxPenalty,n)
              ridge.grid <- npRidgeSequenceAdditive(n.train = n, cap = 1.0)
              ridge <- rep.int(ridge.grid[1L], n)
              ridge.idx <- rep.int(1L, n)
              doridge <- rep.int(TRUE, n)

              nc <- ncol(tww[-1,-1,1])

              while(any(doridge)){
                iloo <- which(doridge)
                for (ii in iloo) {
                  doridge[ii] <- FALSE
                  ridge.val <- ridge[ii]*tww[-1,1,ii][1]/NZD(tww[-1,-1,ii][1,1])
                  beta.ii <- tryCatch(
                    solve(tww[-1,-1,ii] + diag(rep(ridge[ii], nc)),
                          tww[-1,1,ii] + c(ridge.val, rep(0, nc - 1))),
                    error = function(e) e
                  )
                  if (inherits(beta.ii, "error")) {
                    ridge.idx[ii] <- ridge.idx[ii] + 1L
                    if (ridge.idx[ii] <= length(ridge.grid)) {
                      ridge[ii] <- ridge.grid[ridge.idx[ii]]
                      doridge[ii] <- TRUE
                    }
                    beta.ii <- rep(maxPenalty, nc)
                  }
                  mean.loo[ii] <- W[ii,, drop = FALSE] %*% beta.ii
                }
              }
            } else {
              coef.loo <- lp_full_coef(sbw = sbw, leave.one.out.eval = TRUE)
              mean.loo <- rowSums(W * t(coef.loo))
            }

            stopifnot(all(is.finite(mean.loo)))

            if(!any(mean.loo == maxPenalty)){
              fv <- sum((ydat-mean.loo)^2)/n
              cv_progress_step()
            } else {
              cv_progress_step(ridging = TRUE)
              fv <- maxPenalty
            }

            if (isTRUE(cv_state$objective_fast))
              cv_state$fast_total <- cv_state$fast_total + 1L

            return((if (is.finite(fv)) fv else maxPenalty))

          }

          scoef.loo.args <- list(
            bws = bws, txdat = xdat, tydat = ydat,
            leave.one.out = TRUE, iterate = TRUE,
            maxiter = backfit.maxiter, tol = backfit.tol,
            betas = TRUE
          )
          if (!miss.z)
            scoef.loo.args$tzdat <- zdat
          
          partial.cv.ls <- function(param, partial.index) {
            cv_state$objective_fast <- FALSE
            sbw <- apply_bw_to_scbw(bws, param)

            if (!validateBandwidthTF(sbw) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            cv_state$objective_fast <- npscoef_fast_eligible(sbw)
            
            if (backfit.iterate){
              bws$bw.fitted[,partial.index] <- param
              scoef.loo <- do.call(npscoef, scoef.loo.args)
              partial.loo <- W[,partial.index]*scoef.loo$beta[,partial.index]
            } else {
              wj <- W[,partial.index]
              if (identical(reg.engine, "lc")) {
                tww <- npksum(txdat=zdat,
                              tydat=cbind(partial.orig * wj, wj * wj),
                              weights=cbind(partial.orig * wj, 1),
                              bws=sbw,
                              leave.one.out=TRUE)$ksum

                partial.loo <- wj * tww[1,2,]/NZD(tww[2,2,])
              } else {
                partial.loo <- wj * lp_partial_coef(
                  sbw = sbw,
                  wj = wj,
                  partial.y = partial.orig,
                  leave.one.out.eval = TRUE
                )
              }
            }
            

            fv <- sum((partial.orig - partial.loo)^2)/n

            if (isTRUE(cv_state$objective_fast))
              cv_state$fast_total <- cv_state$fast_total + 1L
            
            partial_progress_step(fv = fv)
            return((if (is.finite(fv)) fv else maxPenalty))
          }

          ## Now we implement multistarting

          fval.min <- .Machine$double.xmax
          have_best <- FALSE
          numimp <- 0
          value.overall <- numeric(nmulti)
          num.feval.overall <- 0

          x.scale <- sapply(seq_len(bws$ndim), function(i){
            if (dati$icon[i]){
              return(1.059224*((if (bws$scaling) 1.0 else mysd[sum(dati$icon[seq_len(i)])]*nconfac)))
            }
            
            if (dati$iord[i])
              return(0.5*oMaxL(dati$all.nlev[[i]], kertype = bws$okertype)*
                     (if (bws$scaling) ncatfac else 1.0))
            
            if (dati$iuno[i])
              return(0.5*uMaxL(dati$all.nlev[[i]], kertype = bws$ukertype)*
                     (if (bws$scaling) ncatfac else 1.0))       
          })

          optim.control <- list(abstol = optim.abstol,
                                reltol = optim.reltol,
                                maxit = optim.maxit)

          for (i in seq_len(nmulti)) {

            cv_state$multistart_index <- i
            cv_progress_begin()
            
            if (i == 1) {
              tbw <- .npscoef_default_start_bandwidth(param = x.scale, bwtype = bws$type, nobs = n)
              if (all(bws$bw != 0) &&
                  .npscoef_candidate_is_admissible(param = bws$bw, bwtype = bws$type, nobs = n)) {
                tbw <- .npscoef_finalize_bandwidth(
                  param = bws$bw,
                  bwtype = bws$type,
                  nobs = n,
                  where = "npscoefbw"
                )
              }
            } else {
              tbw <- .npscoef_random_start_bandwidth(param = x.scale, bwtype = bws$type, nobs = n)
            }

            suppressWarnings(optim.return <- optim(tbw,
                                                   fn = overall.cv.ls,
                                                   control = optim.control))
            if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
              num.feval.overall <- num.feval.overall + optim.return$counts[1]
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              tbw <- .npscoef_random_start_bandwidth(param = x.scale, bwtype = bws$type, nobs = n)
              optim.control <- lapply(optim.control, '*', 10.0)
              suppressWarnings(optim.return <- optim(tbw,
                                                     fn = overall.cv.ls,
                                                     control = optim.control))
              if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                num.feval.overall <- num.feval.overall + optim.return$counts[1]

            }

            cv_progress_finish()

            value.overall[i] <- optim.return$value

            if (.npscoef_candidate_is_admissible(
              param = optim.return$par,
              bwtype = bws$type,
              nobs = n
            ) && (!have_best || optim.return$value < fval.min)) {
              param <- .npscoef_finalize_bandwidth(
                param = optim.return$par,
                bwtype = bws$type,
                nobs = n,
                where = "npscoefbw"
              )
              min.overall <- optim.return$value
              fval.min <- min.overall ## Added by jracine Jul 22 2010
              numimp.overall <- numimp + 1
              best.overall <- i
              have_best <- TRUE
            }

            .np_progress_bandwidth_multistart_step(done = i, total = nmulti)
          }

          if (!have_best) {
            if (identical(bws$type, "fixed")) {
              stop("npscoefbw: no feasible fixed bandwidths found", call. = FALSE)
            }
            stop("npscoefbw: no feasible bandwidths found", call. = FALSE)
          }

          param.overall <- bws$bw <- .npscoef_finalize_bandwidth(
            param = param,
            bwtype = bws$type,
            nobs = n,
            where = "npscoefbw"
          )
          bws <- apply_bw_to_scbw(bws, bws$bw)

          if(cv.iterate){
            n.part <- (ncol(xdat)+1)
            backfit.progress <- .np_progress_begin(
              "Backfitting smooth coefficient bandwidth",
              total = cv.num.iterations
            )
            on.exit(cv_progress_end(backfit.progress), add = TRUE)
            bws$bw.fitted <- matrix(data = bws$bw, nrow = length(bws$bw), ncol = n.part)
            ## obtain matrix of alpha.hat | h0 and beta.hat | h0

            scoef.args <- list(bws = bws, txdat = xdat, tydat = ydat, iterate = FALSE, betas = TRUE)
            if (!miss.z)
              scoef.args$tzdat <- zdat
            scoef <- do.call(npscoef, scoef.args)
            
            resid.full <- ydat - scoef$mean

            
            for (i in seq_len(cv.num.iterations)) {
              backfit.progress <- .np_progress_step(
                state = backfit.progress,
                done = i,
                detail = sprintf("iteration %d of %d", i, cv.num.iterations)
              )

              for (j in seq_len(n.part)) {
                ## estimate partial residuals
                partial.orig <- W[,j] * scoef$beta[,j] + resid.full
                partial_progress_begin(iteration = i, partial.index = j)
                
                ## minimise
                suppressWarnings(optim.return <-
                                 optim(tbw, fn = partial.cv.ls,
                                       control = optim.control,
                                       partial.index = j))
                if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                  num.feval.overall <- num.feval.overall + optim.return$counts[1]
                partial_progress_finish(fv = optim.return$value)
                
                ## grab parameter
                bws$bw.fitted[,j] <- optim.return$par

                if (backfit.iterate){
                  ## re-estimate all betas
                  scoef.args <- list(
                    bws = bws, txdat = xdat, tydat = ydat,
                    iterate = TRUE, maxiter = backfit.maxiter,
                    tol = backfit.tol, betas = TRUE
                  )
                  if (!miss.z)
                    scoef.args$tzdat <- zdat
                  scoef <- do.call(npscoef, scoef.args)
                  resid.full <- ydat - scoef$mean
                } else {
                  bws$bw <- bws$bw.fitted[,j]
                  ## estimate new beta.hats

                  bws <- apply_bw_to_scbw(bws, bws$bw)

                  if (identical(reg.engine, "lc")) {
                    wj <- W[,j]
                    tww <- npksum(txdat=zdat,
                                  tydat=cbind(partial.orig * wj, wj * wj),
                                  weights=cbind(partial.orig * wj, 1),
                                  bws=bws)$ksum
                    scoef$beta[,j] <- tww[1,2,]/NZD(tww[2,2,])
                  } else {
                    wj <- W[,j]
                    scoef$beta[,j] <- lp_partial_coef(
                      sbw = bws,
                      wj = wj,
                      partial.y = partial.orig,
                      leave.one.out.eval = FALSE
                    )
                  }
                  
                  bws$bw <- param.overall
                  bws <- apply_bw_to_scbw(bws, bws$bw)
                  ## estimate new full residuals 
                  resid.full <- partial.orig - W[,j] * scoef$beta[,j]
                }
              }
            }
            scoef.loo.args <- list(
              bws = bws, txdat = xdat, tydat = ydat,
              iterate = TRUE, maxiter = backfit.maxiter,
              tol = backfit.tol, leave.one.out = TRUE
            )
            if (!miss.z)
              scoef.loo.args$tzdat <- zdat
            scoef.loo <- do.call(npscoef, scoef.loo.args)$mean
            bws$fval.fitted <- sum((ydat - scoef.loo)^2)/n
          }

          bws$fval = min.overall
          bws$ifval = best.overall
          bws$num.feval = num.feval.overall
          bws$num.feval.fast = cv_state$fast_total
          bws$numimp = numimp.overall
          bws$fval.vector = value.overall
        }
      })[["elapsed"]]
    
    bws$sfactor <- bws$bandwidth <- bws$bw
    nfactor <- nrow^(-2.0/(2.0*bws$ckerorder+bws$ncon))
    dfactor <- EssDee(zdat[, dati$icon, drop = FALSE])*nrow^(-1.0/(2.0*bws$ckerorder+sum(dati$icon)))

    if (bws$scaling) {
      bws$bandwidth[dati$icon] <- bws$bandwidth[dati$icon]*dfactor

      if(bws$nuno > 0)
        bws$bandwidth[dati$iuno] <- bws$bandwidth[dati$iuno]*nfactor

      if(bws$nord > 0)
        bws$bandwidth[dati$iord] <- bws$bandwidth[dati$iord]*nfactor
      
    } else {
      bws$sfactor[dati$icon] <- bws$sfactor[dati$icon]/dfactor

      if(bws$nuno > 0)
        bws$sfactor[dati$iuno] <- bws$sfactor[dati$iuno]/nfactor

      if(bws$nord > 0)
        bws$sfactor[dati$iord] <- bws$sfactor[dati$iord]/nfactor
    }

    ## Restore seed

    .np_seed_exit(seed.state)
    
    bws <- scbandwidth(bw = bws$bw,
                       regtype = regtype,
                       basis = if (is.null(bws$basis)) "glp" else bws$basis,
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
                       fval = bws$fval,
                       ifval = bws$ifval,
                       num.feval = bws$num.feval,
                       num.feval.fast = bws$num.feval.fast,
                       numimp = bws$numimp,
                       fval.vector = bws$fval.vector,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       zdati = bws$zdati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       znames = bws$znames,
                       sfactor = bws$sfactor,
                       bandwidth = bws$bandwidth,
                       rows.omit = rows.omit,
                       bandwidth.compute = bandwidth.compute,
                       optim.method = optim.method,
                       total.time = total.time)

    bws
  }

npscoefbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws,
           backfit.iterate,
           backfit.maxiter,
           backfit.tol,
           bandwidth.compute = TRUE,
           basis,
           bernstein.basis,
           bwmethod,
           bwscaling,
           bwtype,
           ckerbound,
           ckerlb,
           ckerorder,
           ckertype,
           ckerub,
           cv.iterate,
           cv.num.iterations,
           degree,
           nmulti,
           okertype,
           optim.abstol,
           optim.maxattempts,
           optim.maxit,
           optim.method,
           optim.reltol,
           random.seed,
           regtype,
           ukertype,
           ...){
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (!missing(bwmethod) && identical(match.arg(bwmethod, c("cv.ls", "manual")), "manual") &&
        missing(bws))
      stop("bwmethod='manual' requires argument 'bws'")


    miss.z <- missing(zdat)
    xdat <- toFrame(xdat)
    
    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector or a factor")

    if(!miss.z)
      zdat <- toFrame(zdat)

    ## first grab dummy args for scbandwidth() and perform 'bootstrap'
    ## bandwidth call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("regtype", "basis", "degree", "bernstein.basis",
               "bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    sbw.args <- list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      zdati = untangle(zdat),
      xnames = names(xdat),
      ynames = deparse(substitute(ydat)),
      znames = names(zdat),
      bandwidth.compute = bandwidth.compute
    )
    if (any.m) {
      nms <- mc.names[m]
      sbw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    tbw <- do.call(scbandwidth, sbw.args)

    ## next grab dummies for actual bandwidth selection and perform call
    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("zdat", "bandwidth.compute",
               "nmulti",
               "random.seed",
               "cv.iterate",
               "cv.num.iterations",
               "backfit.iterate",
               "backfit.maxiter",
               "backfit.tol",
               "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)


    scbw.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
    if (any.m) {
      nms <- mc.names[m]
      scbw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    tbw <- .np_progress_select_bandwidth_enhanced(
      "Selecting smooth coefficient bandwidth",
      do.call(npscoefbw.scbandwidth, scbw.args)
    )

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
    
  }
