npscoef <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npscoef",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npscoef",bws$call)
        else if (!is.call(bws))
          UseMethod("npscoef",bws)
        else
          UseMethod("npscoef",NULL)
      } else {
        UseMethod("npscoef", NULL)
      }
    } else {
      UseMethod("npscoef", NULL)
    }
  }

npscoef.formula <-
  function(bws, data = NULL, newdata = NULL, y.eval = FALSE, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    response.name <- attr(tmf, "names")[attr(attr(tmf, "terms"), "response")]
    tydat <- model.response(tmf)
    txdat <- tmf[, bws$chromoly[[2]], drop = FALSE]
    miss.z <- !(length(bws$chromoly) == 3)
    if (!miss.z)
      tzdat <- tmf[, bws$chromoly[[3]], drop = FALSE]

    has.eval <- !is.null(newdata)
    if (has.eval) {
      if (!y.eval){
        tt <- delete.response(tt)

        orig.ts <- .np_terms_ts_mask(terms_obj = tt, data = newdata)
        
        ## delete.response clobbers predvars, which is used for timeseries objects
        ## so we need to reconstruct it

        if(all(orig.ts)){
          args <- (as.list(attr(tt, "variables"))[-1])
          attr(tt, "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
        }else if(any(orig.ts)){
          arguments <- (as.list(attr(tt, "variables"))[-1])
          arguments.normal <- arguments[which(!orig.ts)]
          arguments.timeseries <- arguments[which(orig.ts)]

          ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
          attr(tt, "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
        }else{
          attr(tt, "predvars") <- attr(tt, "variables")
        }
      }
      
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      emf <- umf

      if (y.eval)
        eydat <- model.response(emf)
      
      exdat <- emf[, bws$chromoly[[2]], drop = FALSE]
      if (!miss.z)
        ezdat <- emf[, bws$chromoly[[3]], drop = FALSE]
    }


    sc.args <- list(txdat = txdat, tydat = tydat)
    if (!miss.z)
      sc.args$tzdat <- tzdat
    if (has.eval) {
      sc.args$exdat <- exdat
      if (y.eval)
        sc.args$eydat <- eydat
      if (!miss.z)
        sc.args$ezdat <- ezdat
    }
    sc.args$bws <- bws
    ev <- do.call(npscoef, c(sc.args, list(...)))

    if (length(response.name) == 1L && !is.na(response.name) && nzchar(response.name)) {
      if (!is.null(ev$bws))
        ev$bws$ynames <- response.name
    }

    ev$omit <- attr(umf,"na.action")
    ev$rows.omit <- as.vector(ev$omit)
    ev$nobs.omit <- length(ev$rows.omit)

    ev$mean <- napredict(ev$omit, ev$mean)
    ev$merr <- napredict(ev$omit, ev$merr)

    if(ev$residuals){
        ev$resid <- naresid(ev$omit, ev$resid)
    }    
    return(ev)
  }

npscoef.call <-
  function(bws, ...) {
    call.args <- list(
      txdat = .np_eval_bws_call_arg(bws, "xdat"),
      tydat = .np_eval_bws_call_arg(bws, "ydat")
    )
    if (!is.null(bws$zdati))
      call.args$tzdat <- .np_eval_bws_call_arg(bws, "zdat")
    call.args$bws <- bws
    do.call(npscoef, c(call.args, list(...)))
  }

.np_scoef_fit_progress_begin <- function(handoff = FALSE, detail = NULL) {
  state <- .np_progress_begin(
    "Fitting smooth coefficient model",
    surface = "bandwidth"
  )

  if (isTRUE(handoff)) {
    state <- .np_progress_show_now(
      state = state,
      detail = detail
    )
  }

  state
}

.npRmpi_npscoef_should_localize <- function(bws) {
  isa(bws, "scbandwidth") && identical(bws$type, "adaptive_nn")
}

npscoef.default <- function(bws, txdat, tydat, tzdat, nomad = FALSE, ...) {
  .npRmpi_require_active_slave_pool(where = "npscoef()")
  explicit.scbandwidth <- (!missing(bws)) && inherits(bws, "scbandwidth")
  nomad <- npValidateScalarLogical(nomad, "nomad")
  degree.select.value <- if (isTRUE(nomad)) {
    "coordinate"
  } else if ("degree.select" %in% names(list(...))) {
    match.arg(list(...)$degree.select, c("manual", "coordinate", "exhaustive"))
  } else {
    "manual"
  }
  if (!missing(bws) &&
      .npRmpi_npscoef_should_localize(bws) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
  if (.npRmpi_autodispatch_active() &&
      (explicit.scbandwidth || identical(degree.select.value, "manual")))
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")
  tzdat.named <- any(sc.names == "tzdat")

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  no.tzdat <- missing(tzdat)
  has.explicit.bws <- (!no.bws) && isa(bws, "scbandwidth")

  ## if bws was passed in explicitly, do not compute bandwidths

  if(txdat.named)
    txdat <- toFrame(txdat)

  ## if(tydat.named)
  ## tydat <- toFrame(tydat)

  if(tzdat.named)
    tzdat <- toFrame(tzdat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npscoefbw)

  bws.formula <- (!no.bws) && inherits(bws, "formula")
  if (bws.formula) {
    ib <- match("bws", names(sc.bw), nomatch = 0L)
    if (ib > 0L) names(sc.bw)[ib] <- "formula"
  }

  if(bws.named && !bws.formula){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat','tzdat')
  nstxy <- c('xdat','ydat','zdat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }
    
  use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
    sc.bw,
    caller_env = parent.frame()
  )

  tbw <- if (!has.explicit.bws) {
    if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting smooth coefficient bandwidth",
        .np_eval_bw_call(sc.bw, caller_env = parent.frame())
      )
    } else {
      .np_eval_bw_call(sc.bw, caller_env = parent.frame())
    }
  } else {
    .np_eval_bw_call(sc.bw, caller_env = parent.frame())
  }

  ## because of some ambiguities in how the function might be called
  ## we only drop up to two unnamed arguments, when sometimes dropping
  ## three would be appropriate.  also, for simplicity, we don't allow
  ## for inconsistent mixes of named/unnamed arguments so bws is named
  ## or unnamed, and t[xyz]dat collectively either named or unnamed

  call.args <- list(bws = tbw)
  if (no.bws) {
    call.args$txdat <- txdat
    call.args$tydat <- tydat
    if (!no.tzdat) call.args$tzdat <- tzdat
  } else {
    if (txdat.named) call.args$txdat <- txdat
    if (tydat.named) call.args$tydat <- tydat
    if (tzdat.named) call.args$tzdat <- tzdat
    if ((!bws.named) && (!txdat.named) && (!no.tydat) && (!tydat.named)) {
      call.args <- c(call.args, list(tydat))
    }
    if ((!bws.named) && (!txdat.named) && (!no.tzdat) && (!tzdat.named)) {
      call.args <- c(call.args, list(tzdat))
    }
  }
  if (!has.explicit.bws)
    call.args$.np_fit_progress_handoff <- TRUE
  do.call(npscoef, c(call.args, list(...)))

}

.np_scoef_fit_internal <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           tzdat = NULL,
           exdat,
           eydat,
           ezdat,
           betas = FALSE,
           errors = TRUE,
           iterate = TRUE,
           leave.one.out = FALSE,
           maxiter = 100,
           residuals = FALSE,
           tol = .Machine$double.eps,
           ...){

    fit.start <- proc.time()[3]
    residuals <- npValidateScalarLogical(residuals, "residuals")
    errors <- npValidateScalarLogical(errors, "errors")
    iterate <- npValidateScalarLogical(iterate, "iterate")
    leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
    betas <- npValidateScalarLogical(betas, "betas")
    if (!is.numeric(maxiter) || length(maxiter) != 1L || is.na(maxiter) ||
        !is.finite(maxiter) || maxiter < 1 || maxiter != floor(maxiter))
      stop("'maxiter' must be a positive integer")
    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol < 0)
      stop("'tol' must be a finite numeric scalar >= 0")
    maxiter <- as.integer(maxiter)
    tol <- as.double(tol)
    regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
    dots <- list(...)
    fit.progress.allow <- !isFALSE(dots$.np_fit_progress_allow) &&
      isTRUE(.np_progress_enabled(domain = "bandwidth"))
    fit.progress.handoff <- fit.progress.allow && isTRUE(dots$.np_fit_progress_handoff)
    fit.progress <- NULL
    fit.progress.active <- FALSE
    fit.progress.step <- NULL
    if (isTRUE(fit.progress.allow)) {
      fit.progress <- .np_scoef_fit_progress_begin(
        handoff = fit.progress.handoff,
        detail = if (fit.progress.handoff) "building moments" else NULL
      )
      fit.progress.active <- TRUE
      fit.progress.counter <- 0L
      fit.progress.step <- function(detail) {
        fit.progress.counter <<- fit.progress.counter + 1L
        fit.progress <<- .np_progress_step(
          fit.progress,
          done = fit.progress.counter,
          detail = detail
        )
        invisible(NULL)
      }
      on.exit({
        if (isTRUE(fit.progress.active))
          .np_progress_abort(fit.progress)
      }, add = TRUE)
    }

    miss.z <- missing(tzdat) || is.null(tzdat)

    miss.ex = missing(exdat)
    miss.ey = missing(eydat)

    ## if miss.ex then if !miss.ey then ey and tx must match, to get
    ## oos errors alternatively if miss.ey you get is errors if
    ## !miss.ex then if !miss.ey then ey and ex must match, to get oos
    ## errors alternatively if miss.ey you get NO errors since we
    ## don't evaluate on the training data

    txdat <- toFrame(txdat)

    if (!(is.vector(tydat) || is.factor(tydat)))
      stop("'tydat' must be a vector or a factor")

    if (!miss.z)
      tzdat <- toFrame(tzdat)

    if (!miss.ex){
      exdat <- toFrame(exdat)

      if (!miss.z)
        ezdat <- toFrame(ezdat)

      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")

      if (!miss.ey){
        if (dim(exdat)[1] != length(eydat))
          stop("number of evaluation data 'exdat' and dependent data 'eydat' do not match")
      }

    } else if(!miss.ey) {
      if (dim(txdat)[1] != length(eydat))
        stop("number of training data 'txdat' and dependent data 'eydat' do not match")
    }

    if(iterate && !is.null(bws$bw.fitted) && !miss.ex){
      .np_warning("iteration is not supported for out of sample evaluations; using overall bandwidths")
      iterate = FALSE
    }

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(txdat))
    train.df <- data.frame(txdat, tydat)
    if (!miss.z)
      train.df <- data.frame(train.df, tzdat)
    rows.omit <- attr(na.omit(train.df), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows]
    if (!miss.z)
      tzdat <- tzdat[keep.rows,, drop = FALSE]

    if (!miss.ex){
      keep.eval <- rep_len(TRUE, nrow(exdat))
      eval.df <- data.frame(exdat)
      if (!miss.ey)
        eval.df <- data.frame(eval.df, eydat)
      if (!miss.z)
        eval.df <- data.frame(eval.df, ezdat)
      rows.omit <- attr(na.omit(eval.df), "na.action")
      if (length(rows.omit) > 0L)
        keep.eval[as.integer(rows.omit)] <- FALSE

      exdat <- exdat[keep.eval,,drop = FALSE]
      if (!miss.ey)
        eydat <- eydat[keep.eval]
      if (!miss.z)
        ezdat <- ezdat[keep.eval,, drop = FALSE]

      if (!any(keep.eval))
        stop("Evaluation data has no rows without NAs")
    }

    ## convert tydat, eydat to numeric, from a factor with levels from the y-data
    ## used during bandwidth selection.

    if (is.factor(tydat)){
      tydat <- adjustLevels(as.data.frame(tydat), bws$ydati)[,1]
      tydat <- (bws$ydati$all.dlev[[1]])[as.integer(tydat)]
    }
    else
      tydat <- as.double(tydat)

    if (miss.ey)
      eydat <- double()
    else {
      if (is.factor(eydat)){
        eydat <- adjustLevels(as.data.frame(eydat), bws$ydati)[,1]
        eydat <- (bws$ydati$all.dlev[[1]])[as.integer(eydat)]
      }
      else
        eydat <- as.double(eydat)
    }

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.

    txdat <- adjustLevels(txdat, bws$xdati)
    if (!miss.z)
      tzdat <- adjustLevels(tzdat, bws$zdati)

    if (!miss.ex){
      exdat <- adjustLevels(exdat, bws$xdati)
      if (!miss.z)
        ezdat <- adjustLevels(ezdat, bws$zdati)
    }

    ## grab the evaluation data before it is converted to numeric
    if(miss.ex){
      teval <- txdat
      if (!miss.z)
        teval <- list(exdat = txdat, ezdat = tzdat)
    } else {
      teval <- exdat
      if (!miss.z)
        teval <- list(exdat = exdat, ezdat = ezdat)
    }

    ## put the unordered, ordered, and continuous data in their own objects
    ## data that is not a factor is continuous.

    txdat <- toMatrix(txdat)

    if (!miss.ex){
      exdat <- toMatrix(exdat)
    }

    if (miss.z){
      tzdat <- txdat
      if (!miss.ex)
        ezdat <- exdat
    }
    ## from this point on txdat and exdat have been recast as matrices
    ## construct 'W' matrix

    spec <- .npscoef_canonical_spec(source = bws, zdat = tzdat, where = "npscoef")
    reg.engine <- spec$regtype.engine
    W.train <- W <- as.matrix(data.frame(1,txdat))
    maxPenalty <- sqrt(.Machine$double.xmax)
    tnrow <- nrow(txdat)
    enrow <- (if (miss.ex) nrow(txdat) else nrow(exdat))

    if (!miss.ex)
      W <- as.matrix(data.frame(1,exdat))

    safe_chol2inv <- function(a, ridge0, eps, maxiter = 1000L){
      nc.local <- ncol(a)
      I.local <- diag(rep(1.0, nc.local))
      ridge.local <- max(as.double(ridge0), 0.0)
      for (iter in seq_len(maxiter)) {
        cm <- tryCatch(
          chol2inv(chol(a + ridge.local * I.local)),
          error = function(e) NULL
        )
        if (!is.null(cm))
          return(cm)
        ridge.local <- ridge.local + eps
      }
      NULL
    }

    fast_moment_solve <- function(tww.slice, tyw.slice, ridge.add, ridge.val) {
      ncoef.local <- nrow(tww.slice)

      if (ncoef.local == 1L) {
        denom <- tww.slice[1L, 1L] + ridge.add
        if (!is.finite(denom) || abs(denom) <= .Machine$double.eps)
          return(NULL)
        return((tyw.slice[1L] + ridge.val) / denom)
      }

      if (ncoef.local == 2L) {
        a11 <- tww.slice[1L, 1L] + ridge.add
        a12 <- tww.slice[1L, 2L]
        a21 <- tww.slice[2L, 1L]
        a22 <- tww.slice[2L, 2L] + ridge.add
        b1 <- tyw.slice[1L] + ridge.val
        b2 <- tyw.slice[2L]
        detA <- a11 * a22 - a12 * a21
        scale <- max(abs(c(a11, a12, a21, a22, b1, b2)), 1.0)
        if (!is.finite(detA) || abs(detA) <= .Machine$double.eps * scale)
          return(NULL)
        return(c(
          (b1 * a22 - a12 * b2) / detA,
          (a11 * b2 - b1 * a21) / detA
        ))
      }

      NULL
    }

    solve_moment_system <- function(tyw, tww, W.eval.design, Wz.eval = NULL, progress_detail = NULL) {
      neval.local <- ncol(tyw)
      ncoef <- nrow(tyw)
      pcoef <- ncol(W.eval.design)
      coef.out <- matrix(maxPenalty, nrow = pcoef, ncol = neval.local)
      theta.out <- if (is.null(Wz.eval)) NULL else matrix(NA_real_, nrow = ncoef, ncol = neval.local)
      ridge.grid <- npRidgeSequenceAdditive(n.train = tnrow, cap = 1.0)
      ridge <- rep.int(ridge.grid[1L], neval.local)
      ridge.idx <- rep.int(1L, neval.local)
      doridge <- rep.int(TRUE, neval.local)

      while(any(doridge)){
        iloo <- seq_len(neval.local)[doridge]
        for (ii in iloo) {
          doridge[ii] <- FALSE
          ridge.val <- ridge[ii]*tyw[,ii][1]/NZD(tww[,,ii][1,1])
          theta.ii <- fast_moment_solve(
            tww.slice = tww[, , ii],
            tyw.slice = tyw[, ii],
            ridge.add = ridge[ii],
            ridge.val = ridge.val
          )
          if (is.null(theta.ii)) {
            theta.ii <- tryCatch(
              solve(tww[,,ii] + diag(rep(ridge[ii], ncoef)),
                    tyw[,ii] + c(ridge.val, rep(0, ncoef - 1))),
              error = function(e) e
            )
          }
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
            theta.out[,ii] <- theta.ii
            coef.out[,ii] <- as.vector(crossprod(
              Wz.eval[ii,],
              matrix(theta.ii, nrow = ncol(Wz.eval), ncol = pcoef)
            ))
          }
          if (!is.null(fit.progress.step))
            fit.progress.step(progress_detail)
        }
      }

      list(coef = coef.out, theta = theta.out, ridge = ridge)
    }

    solve_single_moment_system <- function(tyw.vec, tww.mat) {
      ncoef <- length(tyw.vec)
      ridge.grid <- npRidgeSequenceAdditive(n.train = tnrow, cap = 1.0)
      ridge <- ridge.grid[1L]
      ridge.idx <- 1L

      repeat {
        ridge.val <- ridge * tyw.vec[1L] / NZD(tww.mat[1L, 1L])
        theta <- fast_moment_solve(
          tww.slice = tww.mat,
          tyw.slice = tyw.vec,
          ridge.add = ridge,
          ridge.val = ridge.val
        )
        if (is.null(theta)) {
          theta <- tryCatch(
            solve(
              tww.mat + diag(rep(ridge, ncoef)),
              tyw.vec + c(ridge.val, rep(0, ncoef - 1L))
            ),
            error = function(e) e
          )
        }
        if (!inherits(theta, "error"))
          return(list(coef = as.double(theta), ridge = ridge))

        ridge.idx <- ridge.idx + 1L
        if (ridge.idx > length(ridge.grid))
          return(list(coef = rep(maxPenalty, ncoef), ridge = ridge))
        ridge <- ridge.grid[ridge.idx]
      }
    }

    local_npksum <- function(args) {
      .npRmpi_with_local_regression(do.call(npksum, args))
    }
    lp1_local_npksum <- function(...) {
      local_npksum(list(...))
    }

    lc_moments <- function(z.eval, leave.one.out.eval, u2 = NULL) {
      yW.local <- cbind(tydat, W.train)
      ksum.args <- list(
        txdat = tzdat,
        tydat = yW.local,
        weights = yW.local,
        bws = bws,
        leave.one.out = leave.one.out.eval,
        bandwidth.divide = TRUE
      )
      if (!leave.one.out.eval && !is.null(z.eval))
        ksum.args$exdat <- z.eval
      main.ks <- local_npksum(ksum.args)$ksum
      tyw.out <- main.ks[-1L, 1L, , drop = FALSE]
      if (length(dim(tyw.out)) == 3L)
        dim(tyw.out) <- c(dim(tyw.out)[1L], dim(tyw.out)[3L])
      tww.out <- main.ks[-1L, -1L, , drop = FALSE]

      s.out <- NULL
      if (!is.null(u2)) {
        cov.args <- list(
          txdat = tzdat,
          tydat = W.train,
          weights = W.train * as.double(u2),
          bws = bws,
          leave.one.out = leave.one.out.eval,
          bandwidth.divide = TRUE,
          kernel.pow = 2
        )
        if (!leave.one.out.eval && !is.null(z.eval))
          cov.args$exdat <- z.eval
        s.out <- local_npksum(cov.args)$ksum
      }

      list(tyw = tyw.out, tww = tww.out, s = s.out)
    }

    lp_state <- if (identical(reg.engine, "lp")) {
      .npscoef_lp_state(
        bws = bws,
        tzdat = tzdat,
        ezdat = if (miss.ex) tzdat else ezdat,
        leave.one.out = leave.one.out,
        where = "npscoef"
      )
    } else {
      NULL
    }

    do.iterate <- (iterate && !is.null(bws$bw.fitted) && miss.ex && identical(reg.engine, "lc"))
    gate.zdat <- if (miss.ex) tzdat else rbind(tzdat, ezdat)
    fast.largeh.lc <- identical(reg.engine, "lc") &&
      !leave.one.out &&
      identical(bws$type, "fixed") &&
      !do.iterate &&
      .npscoefbw_fast_eligible(bws, eval.zdat = gate.zdat)
    fast.largeh.lp1 <- identical(reg.engine, "lp") &&
      !leave.one.out &&
      identical(bws$type, "fixed") &&
      identical(spec$basis.engine, "glp") &&
      !isTRUE(spec$bernstein.basis.engine) &&
      all(as.integer(spec$degree.engine) == 1L) &&
      .npscoefbw_fast_eligible(bws, eval.zdat = gate.zdat)
    fast.largeh <- fast.largeh.lc || fast.largeh.lp1

    lc_fast_global_moments <- function(z.eval.one, u2 = NULL) {
      yW.local <- cbind(tydat, W.train)
      ksum.args <- list(
        txdat = tzdat,
        tydat = yW.local,
        weights = yW.local,
        exdat = z.eval.one,
        bws = bws,
        bandwidth.divide = TRUE
      )
      main.ks <- local_npksum(ksum.args)$ksum
      out <- list(
        tyw = as.double(main.ks[-1L, 1L, 1L]),
        tww = main.ks[-1L, -1L, 1L, drop = TRUE]
      )

      if (!is.null(u2)) {
        cov.args <- list(
          txdat = tzdat,
          tydat = W.train,
          weights = W.train * as.double(u2),
          exdat = z.eval.one,
          bws = bws,
          bandwidth.divide = TRUE,
          kernel.pow = 2
        )
        out$s <- local_npksum(cov.args)$ksum[, , 1L, drop = TRUE]
      } else {
        out$s <- NULL
      }

      out
    }

    lp_tensor_moments <- function(state, u2 = NULL) {
      tensor.train <- .npscoef_row_tensor_design(W.train, state$W.train)
      ytensor <- cbind(tydat, tensor.train)
      ksum.args <- list(
        txdat = state$z.train,
        tydat = ytensor,
        weights = ytensor,
        bws = state$rbw,
        leave.one.out = state$leave.one.out,
        bandwidth.divide = TRUE
      )
      if (!state$leave.one.out)
        ksum.args$exdat <- state$z.eval
      main.ks <- local_npksum(ksum.args)$ksum
      tyw.out <- main.ks[-1L, 1L, , drop = FALSE]
      if (length(dim(tyw.out)) == 3L)
        dim(tyw.out) <- c(dim(tyw.out)[1L], dim(tyw.out)[3L])
      tww.out <- main.ks[-1L, -1L, , drop = FALSE]

      s.out <- NULL
      if (!is.null(u2)) {
        cov.args <- list(
          txdat = state$z.train,
          tydat = tensor.train,
          weights = tensor.train * as.double(u2),
          bws = state$rbw,
          leave.one.out = state$leave.one.out,
          bandwidth.divide = TRUE,
          kernel.pow = 2
        )
        if (!state$leave.one.out)
          cov.args$exdat <- state$z.eval
        s.out <- local_npksum(cov.args)$ksum
      }

      list(tyw = tyw.out, tww = tww.out, s = s.out)
    }

    moments <- NULL
    if (fast.largeh.lc) {
      eval.z.one <- if (miss.ex) tzdat[1L, , drop = FALSE] else ezdat[1L, , drop = FALSE]
      fast.eval <- lc_fast_global_moments(z.eval.one = eval.z.one)
      if (!is.null(fit.progress.step))
        fit.progress.step("solving global coefficients")
      fast.solve <- solve_single_moment_system(
        tyw.vec = fast.eval$tyw,
        tww.mat = fast.eval$tww
      )
      if (!is.null(fit.progress.step))
        fit.progress.step("assembling fitted values")
      coef.vec <- fast.solve$coef
      coef.mat <- matrix(coef.vec, nrow = length(coef.vec), ncol = enrow)
      ridge <- rep.int(fast.solve$ridge, enrow)
    } else if (fast.largeh.lp1) {
      fast.eval <- .npscoef_lp1_largeh_global_fit(
        bws = bws,
        tzdat = tzdat,
        ezdat = if (miss.ex) tzdat else ezdat,
        W.train = W.train,
        tydat = tydat,
        leave.one.out = leave.one.out,
        where = "npscoef",
        solver = solve_single_moment_system,
        ksum_fun = lp1_local_npksum
      )
      if (!is.null(fit.progress.step))
        fit.progress.step("solving global coefficients")
      coef.mat <- fast.eval$coef
      ridge <- rep.int(fast.eval$ridge, enrow)
    } else if (identical(reg.engine, "lc")) {
      moments <- lc_moments(
        z.eval = if (miss.ex) NULL else ezdat,
        leave.one.out.eval = leave.one.out
      )
      solver <- solve_moment_system(
        tyw = moments$tyw,
        tww = moments$tww,
        W.eval.design = W,
        progress_detail = "solving coefficient rows"
      )
    } else {
      moments <- lp_tensor_moments(lp_state)
      solver <- solve_moment_system(
        tyw = moments$tyw,
        tww = moments$tww,
        W.eval.design = W,
        Wz.eval = lp_state$W.eval,
        progress_detail = "solving coefficient rows"
      )
    }

    if (!fast.largeh) {
      coef.mat <- solver$coef
      ridge <- solver$ridge
    }

    if (iterate && !is.null(bws$bw.fitted) && miss.ex && !identical(reg.engine, "lc"))
      .np_warning("iterate=TRUE currently supports regtype='lc' for npscoef; using iterate=FALSE")
    if (do.iterate){
      resid <- tydat - sapply(seq_len(enrow), function(i) { W[i,, drop = FALSE] %*% coef.mat[,i] })

      i = 0
      max.err <- .Machine$double.xmax
      aydat <- abs(tydat) + .Machine$double.eps

      n.part <- (ncol(txdat)+1)

      while((max.err > tol) && ((i <- i + 1) <= maxiter)){
        resid.old <- resid
        for (j in seq_len(n.part)) {
          partial <- W[,j] * coef.mat[j,] + resid

          twww <- .npRmpi_with_local_regression(npksum(
            txdat=tzdat,
            tydat=cbind(partial * W[,j],W[,j]^2),
            weights=cbind(partial * W[,j],1),
            bws=bws,
            leave.one.out=leave.one.out
          ))$ksum

          coef.mat[j,] <- twww[1,2,]/NZD(twww[2,2,])
          resid <- partial - W[,j] * coef.mat[j,]
          if (!is.null(fit.progress.step))
            fit.progress.step(sprintf("backfit cycle %d partial %d/%d", i, j, n.part))
        }
        max.err <- max(abs(resid.old - resid)/aydat)
      }
      if (max.err > tol)
        .np_warning(paste("backfit iterations did not converge. max err= ", max.err,", tol= ", tol,", maxiter= ", maxiter, sep=''))
      mean <- tydat - resid
    } else {
      mean <- sapply(seq_len(enrow), function(i) { W[i,, drop = FALSE] %*% coef.mat[,i] })
    }

    if (!miss.ey) {
      RSQ = RSQfunc(eydat, mean)
      MSE = MSEfunc(eydat, mean)
      MAE = MAEfunc(eydat, mean)
      MAPE = MAPEfunc(eydat, mean)
      CORR = CORRfunc(eydat, mean)
      SIGN = SIGNfunc(eydat, mean)
    } else if(miss.ex) {
      RSQ = RSQfunc(tydat, mean)
      MSE = MSEfunc(tydat, mean)
      MAE = MAEfunc(tydat, mean)
      MAPE = MAPEfunc(tydat, mean)
      CORR = CORRfunc(tydat, mean)
      SIGN = SIGNfunc(tydat, mean)
    }

    if (errors || (residuals && miss.ex)) {
      if (errors) {
        if (fast.largeh.lc) {
          if (miss.ex) {
            train.solve <- fast.solve
            mean.fit <- mean
          } else {
            train.fast <- lc_fast_global_moments(z.eval.one = tzdat[1L, , drop = FALSE])
            train.solve <- solve_single_moment_system(
              tyw.vec = train.fast$tyw,
              tww.mat = train.fast$tww
            )
            mean.fit <- as.vector(W.train %*% train.solve$coef)
          }
          resid <- tydat - mean.fit
          u2.W <- resid^2
          s.fast <- lc_fast_global_moments(
            z.eval.one = if (miss.ex) tzdat[1L, , drop = FALSE] else ezdat[1L, , drop = FALSE],
            u2 = u2.W
          )$s
          if (!is.null(fit.progress.step))
            fit.progress.step("estimating standard errors")
          cm.fast <- safe_chol2inv(fast.eval$tww, fast.solve$ridge, 1.0 / nrow(txdat))
          merr <- rep(NA_real_, enrow)
          beta.se <- matrix(NA_real_, nrow = enrow, ncol = nrow(coef.mat))
          if (!is.null(cm.fast)) {
            vcv.beta.fast <- cm.fast %*% s.fast %*% cm.fast
            merr <- sqrt(pmax(rowSums((W %*% vcv.beta.fast) * W), 0.0))
            beta.se[] <- rep(sqrt(pmax(diag(vcv.beta.fast), 0.0)), each = enrow)
          }
        } else if (fast.largeh.lp1) {
          if (miss.ex) {
            train.fast <- fast.eval
            mean.fit <- mean
          } else {
            train.fast <- .npscoef_lp1_largeh_global_fit(
              bws = bws,
              tzdat = tzdat,
              ezdat = tzdat,
              W.train = W.train,
              tydat = tydat,
              leave.one.out = leave.one.out,
              where = "npscoef",
              solver = solve_single_moment_system,
              ksum_fun = lp1_local_npksum
            )
            mean.fit <- sapply(seq_len(tnrow), function(i) { W.train[i,, drop = FALSE] %*% train.fast$coef[,i] })
          }
          resid <- tydat - mean.fit
          u2.W <- resid^2
          s.fast <- .npscoef_lp1_largeh_global_fit(
            bws = bws,
            tzdat = tzdat,
            ezdat = if (miss.ex) tzdat else ezdat,
            W.train = W.train,
            tydat = tydat,
            u2 = u2.W,
            leave.one.out = leave.one.out,
            where = "npscoef",
            solver = solve_single_moment_system,
            ksum_fun = lp1_local_npksum
          )$s
          if (!is.null(fit.progress.step))
            fit.progress.step("estimating standard errors")
          cm.fast <- safe_chol2inv(fast.eval$tww, fast.eval$ridge, 1.0 / nrow(txdat))
          merr <- rep(NA_real_, enrow)
          beta.se <- matrix(NA_real_, nrow = enrow, ncol = nrow(coef.mat))
          if (!is.null(cm.fast)) {
            vcv.theta.fast <- cm.fast %*% s.fast %*% cm.fast
            for (i in seq_len(enrow)) {
              trans.i <- kronecker(diag(ncol(W)), matrix(fast.eval$lp_state$W.eval[i,], nrow = 1L))
              vcv.beta.fast <- trans.i %*% vcv.theta.fast %*% t(trans.i)
              w.i <- W[i,, drop = FALSE]
              merr[i] <- sqrt(max(drop(w.i %*% vcv.beta.fast %*% t(w.i)), 0.0))
              beta.se[i,] <- sqrt(pmax(diag(vcv.beta.fast), 0.0))
            }
          }
        } else if (miss.ex && !do.iterate) {
          mean.fit <- mean
          resid <- tydat - mean.fit
          u2.W <- resid^2
          if (identical(reg.engine, "lc")) {
            moments$s <- lc_moments(
              z.eval = NULL,
              leave.one.out.eval = leave.one.out,
              u2 = u2.W
            )$s
          } else {
            moments$s <- lp_tensor_moments(lp_state, u2 = u2.W)$s
          }
        } else if (identical(reg.engine, "lc")) {
          train.moments <- lc_moments(z.eval = NULL, leave.one.out.eval = leave.one.out)
          train.solve <- solve_moment_system(
            tyw = train.moments$tyw,
            tww = train.moments$tww,
            W.eval.design = W.train
          )
          mean.fit <- sapply(seq_len(tnrow), function(i) { W.train[i,, drop = FALSE] %*% train.solve$coef[,i] })
          u2.W <- (resid <- tydat - mean.fit)^2
          moments$s <- lc_moments(
            z.eval = if (miss.ex) NULL else ezdat,
            leave.one.out.eval = leave.one.out,
            u2 = u2.W
          )$s
        } else {
          lp_state.err <- .npscoef_lp_state(
            bws = bws,
            tzdat = tzdat,
            ezdat = tzdat,
            leave.one.out = leave.one.out,
            where = "npscoef"
          )
          train.moments <- lp_tensor_moments(lp_state.err)
          train.solve <- solve_moment_system(
            tyw = train.moments$tyw,
            tww = train.moments$tww,
            W.eval.design = W.train,
            Wz.eval = lp_state.err$W.eval
          )
          mean.fit <- sapply(seq_len(tnrow), function(i) { W.train[i,, drop = FALSE] %*% train.solve$coef[,i] })
          u2.W <- (resid <- tydat - mean.fit)^2
          moments$s <- lp_tensor_moments(lp_state, u2 = u2.W)$s
        }
      } else if (residuals && miss.ex) {
        resid <- tydat - mean
      }
    }

    if (!errors)
      beta.se <- NULL
    if(errors && !fast.largeh){
      u2 <- as.double(u2.W)
      merr <- rep(NA_real_, enrow)
      beta.se <- matrix(NA_real_, nrow = enrow, ncol = nrow(coef.mat))
      for (i in seq_len(enrow)) {
        cm <- safe_chol2inv(moments$tww[,,i], ridge[i], 1.0 / nrow(txdat))
        if (is.null(cm))
          next

        s.mat <- moments$s[,,i]

        if (identical(reg.engine, "lc")) {
          vcv.beta <- cm %*% s.mat %*% cm
        } else {
          vcv.theta <- cm %*% s.mat %*% cm
          trans.i <- kronecker(diag(ncol(W)), matrix(lp_state$W.eval[i,], nrow = 1L))
          vcv.beta <- trans.i %*% vcv.theta %*% t(trans.i)
        }

        w.i <- W[i,,drop=FALSE]
        merr[i] <- sqrt(max(drop(w.i %*% vcv.beta %*% t(w.i)), 0.0))
        beta.se[i,] <- sqrt(pmax(diag(vcv.beta), 0.0))
        if (!is.null(fit.progress.step))
          fit.progress.step("estimating standard errors")
      }

    }

    sc.obj.args <- list(
      bws = bws,
      eval = teval,
      mean = mean,
      residuals = residuals,
      betas = betas,
      ntrain = nrow(txdat),
      trainiseval = miss.ex
    )
    if (errors && !do.iterate)
      sc.obj.args$merr <- merr
    if (ncol(txdat) > 0L) {
      sc.obj.args$grad <- t(coef.mat[-1,,drop = FALSE])
      if (errors && !do.iterate && !is.null(beta.se))
        sc.obj.args$gerr <- beta.se[, -1, drop = FALSE]
    }
    if (betas)
      sc.obj.args$beta <- t(coef.mat)
    if (residuals)
      sc.obj.args$resid <- resid
    if (!(miss.ey && !miss.ex))
      sc.obj.args$xtra <- c(RSQ, MSE, MAE, MAPE, CORR, SIGN)
    ev <- do.call(smoothcoefficient, sc.obj.args)
    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)
    ev$timing <- bws$timing
    ev$total.time <- total.time
    ev$optim.time <- optim.time
    ev$fit.time <- fit.elapsed
    ev$nomad.time <- if (!is.null(bws$nomad.time) && is.finite(bws$nomad.time)) as.double(bws$nomad.time) else NA_real_
    ev$powell.time <- if (!is.null(bws$powell.time) && is.finite(bws$powell.time)) as.double(bws$powell.time) else NA_real_
    if (isTRUE(fit.progress.active)) {
      fit.progress <- .np_progress_end(fit.progress)
      fit.progress.active <- FALSE
    }
    ev

  }

npscoef.scbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           tzdat = NULL,
           exdat,
           eydat,
           ezdat,
           betas = FALSE,
           errors = TRUE,
           iterate = TRUE,
           leave.one.out = FALSE,
           maxiter = 100,
           residuals = FALSE,
           tol = .Machine$double.eps,
           ...){
    residuals <- npValidateScalarLogical(residuals, "residuals")
    errors <- npValidateScalarLogical(errors, "errors")
    iterate <- npValidateScalarLogical(iterate, "iterate")
    leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
    betas <- npValidateScalarLogical(betas, "betas")
    if (!is.numeric(maxiter) || length(maxiter) != 1L || is.na(maxiter) ||
        !is.finite(maxiter) || maxiter < 1 || maxiter != floor(maxiter))
      stop("'maxiter' must be a positive integer")
    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol < 0)
      stop("'tol' must be a finite numeric scalar >= 0")
    maxiter <- as.integer(maxiter)
    tol <- as.double(tol)

    .npRmpi_require_active_slave_pool(where = "npscoef()")
    if (.npRmpi_npscoef_should_localize(bws) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
      return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
    dots <- list(...)
    fit.progress.allow <- !isFALSE(dots$.np_fit_progress_allow)
    use.master.fit.progress <- .npRmpi_autodispatch_active() &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      isTRUE(fit.progress.allow) &&
      isTRUE(.np_progress_enabled(domain = "bandwidth"))
    if (!use.master.fit.progress && .npRmpi_autodispatch_active()) {
      result <- .npRmpi_autodispatch_call(match.call(), parent.frame())
      return(.npRmpi_restore_nomad_fit_bws_metadata(result, bws))
    }

    .np_scoef_fit_internal(
      bws = bws,
      txdat = txdat,
      tydat = tydat,
      tzdat = tzdat,
      exdat = exdat,
      eydat = eydat,
      ezdat = ezdat,
      betas = betas,
      errors = errors,
      iterate = iterate,
      leave.one.out = leave.one.out,
      maxiter = maxiter,
      residuals = residuals,
      tol = tol,
      ...
    )
  }
