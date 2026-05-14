npqreg <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npqreg",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npqreg",bws$call)
        else if (!is.call(bws))
          UseMethod("npqreg",bws)
        else
          UseMethod("npqreg",NULL)
      } else {
        UseMethod("npqreg", NULL)
      }
    } else {
      UseMethod("npqreg", NULL)
    }
  }

.npqreg.fit.control.names <- c("data", "newdata", "exdat", "tau", "gradients", "tol", "small", "itmax")
.npqreg.removed.solver.controls <- c("ftol",
                                     "lbc.dir", "dfc.dir", "cfac.dir", "initc.dir",
                                     "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir")

.npqreg_validate_tau <- function(tau) {
  if (!is.numeric(tau) || !length(tau) || anyNA(tau) ||
      any(!is.finite(tau)) || any(tau <= 0) || any(tau >= 1))
    stop("'tau' must contain numeric values in (0,1)")
  as.double(tau)
}

.npqreg_tau_labels <- function(tau) {
  paste0("tau=", format(tau, trim = TRUE, scientific = FALSE))
}

.npqreg_napredict_eval <- function(omit, x) {
  if (!length(omit))
    return(x)
  if (length(dim(x)) <= 2L)
    return(napredict(omit, x))
  d <- dim(x)
  dn <- dimnames(x)
  if (!is.null(dn))
    dn[[1L]] <- NULL
  out <- array(NA_real_,
               dim = c(d[1L] + length(omit), d[-1L]),
               dimnames = dn)
  keep <- seq_len(dim(out)[1L])[-as.integer(omit)]
  out[keep, , ] <- x
  out
}

.npqreg_validate_newdata_terms <- function(newdata, xnames) {
  nd <- toFrame(newdata)
  missing.names <- setdiff(xnames, names(nd))
  if (length(missing.names))
    stop(sprintf(
      "newdata must contain columns: %s",
      paste(shQuote(xnames), collapse = ", ")
    ), call. = FALSE)
  invisible(TRUE)
}

.npqreg_fit_dots <- function(dots, allow.bandwidth.controls = FALSE) {
  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep("", length(dots))

  stale <- intersect(dot.names[nzchar(dot.names)], .npqreg.removed.solver.controls)
  if (length(stale) && !allow.bandwidth.controls) {
    stop(sprintf(
      "'%s' %s no longer accepted by npqreg; the canonical one-dimensional quantile extractor is controlled by 'tol', 'small', and 'itmax'",
      paste(stale, collapse = "', '"),
      if (length(stale) == 1L) "is" else "are"
    ))
  }

  if (!allow.bandwidth.controls) {
    bad <- dot.names == "" | !(dot.names %in% .npqreg.fit.control.names)
    if (any(bad))
      .np_reject_unused_dots(dots[bad], "npqreg")
  }

  keep <- (!nzchar(dot.names)) | (dot.names %in% .npqreg.fit.control.names)
  dots[keep]
}

.npqreg_strip_fit_controls_from_bw_call <- function(call) {
  for (nm in c("tau", "gradients", "tol", "small", "itmax", "newdata", "exdat")) {
    if (nm %in% names(call))
      call[[nm]] <- NULL
  }
  call
}

.npqreg_quantile_delta_from_conditional <- function(bws,
                                                    xdat,
                                                    ydat,
                                                    exdat,
                                                    quantile,
                                                    gradients = FALSE) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  quantile <- as.double(quantile)
  gradients <- npValidateScalarLogical(gradients, "gradients")

  if (length(quantile) != nrow(exdat))
    stop("quantile delta helper requires one quantile per evaluation row")
  if (ncol(ydat) != 1L)
    stop("quantile delta helper requires a single response")

  eydat <- stats::setNames(data.frame(quantile), names(ydat)[1L])
  cdf.obj <- .np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = TRUE,
    gradients = gradients
  )
  dens.obj <- .np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = FALSE,
    gradients = FALSE
  )

  dens <- as.double(dens.obj$condens)
  quanterr <- as.double(cdf.obj$conderr) / NZD(dens)
  quanterr[!is.finite(quanterr) | quanterr < 0.0] <- NA_real_

  if (!gradients) {
    return(list(
      quanterr = quanterr,
      quantgrad = NA,
      quantgerr = NA,
      cdf = cdf.obj,
      dens = dens.obj
    ))
  }

  dens.mat <- matrix(NZD(dens),
                     nrow = nrow(cdf.obj$congrad),
                     ncol = ncol(cdf.obj$congrad))
  grad <- -cdf.obj$congrad / dens.mat
  grad[!is.finite(grad)] <- NA_real_

  gerr <- cdf.obj$congerr / dens.mat
  gerr[!is.finite(gerr) | gerr < 0.0] <- NA_real_

  list(
    quanterr = quanterr,
    quantgrad = grad,
    quantgerr = gerr,
    cdf = cdf.obj,
    dens = dens.obj
  )
}

.npqreg_selected_cdf_values <- function(bws,
                                        xdat,
                                        ydat,
                                        exdat,
                                        ycand) {
  ydat <- toFrame(ydat)
  yname <- names(ydat)[1L]
  eydat <- stats::setNames(data.frame(as.double(ycand)), yname)
  as.double(.np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = TRUE,
    gradients = FALSE
  )$condist)
}

.npqreg_assert_selected_cdf_metadata <- function(bws) {
  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  if (!identical(reg.engine, "lp"))
    return(invisible(TRUE))

  if (is.null(bws$degree.engine) && is.null(bws$degree))
    stop("selected LP conditional distribution metadata missing from bandwidth object: degree")
  invisible(TRUE)
}

.npqreg_invert_selected_cdf <- function(bws,
                                        xdat,
                                        ydat,
                                        exdat,
                                        tau,
                                        tol,
                                        small,
                                        itmax) {
  .npqreg_assert_selected_cdf_metadata(bws)

  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  y <- as.double(ydat[[1L]])
  y <- y[is.finite(y)]
  if (!length(y))
    stop("npqreg selected-CDF inversion requires finite response support")

  n.eval <- nrow(exdat)
  y.min <- min(y)
  y.max <- max(y)
  if (!is.finite(y.min) || !is.finite(y.max))
    stop("npqreg selected-CDF inversion found non-finite response support")
  if (identical(y.min, y.max))
    return(rep.int(y.min, n.eval))

  lo <- rep.int(y.min, n.eval)
  hi <- rep.int(y.max, n.eval)

  flo <- .npqreg_selected_cdf_values(bws, xdat, ydat, exdat, lo)
  fhi <- .npqreg_selected_cdf_values(bws, xdat, ydat, exdat, hi)
  if (any(!is.finite(flo)) || any(!is.finite(fhi)))
    stop("npqreg selected-CDF inversion encountered non-finite bracket values")

  done.low <- flo >= tau
  done.high <- fhi < tau
  active <- !(done.low | done.high)

  maxiter <- min(as.integer(itmax), 1000L)
  iter <- 0L
  while (any(active) && iter < maxiter) {
    iter <- iter + 1L
    mid <- (lo[active] + hi[active]) / 2.0
    fmid <- .npqreg_selected_cdf_values(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat[active, , drop = FALSE],
      ycand = mid
    )
    if (any(!is.finite(fmid)))
      stop("npqreg selected-CDF inversion encountered non-finite refinement values")

    active.idx <- which(active)
    upper <- fmid >= tau
    hi[active.idx[upper]] <- mid[upper]
    lo[active.idx[!upper]] <- mid[!upper]

    width <- hi[active.idx] - lo[active.idx]
    scale <- pmax(abs(hi[active.idx]), abs(lo[active.idx]), 1.0)
    active[active.idx] <- width > (tol * scale + small)
  }

  if (any(active))
    stop("npqreg selected-CDF inversion failed to converge within 'itmax'")

  out <- hi
  out[done.low] <- y.min
  out[done.high] <- y.max
  out
}

npqreg.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    tydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
    txdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]

    has.eval <- !is.null(newdata)
    if (has.eval) {
      .npqreg_validate_newdata_terms(newdata, bws$variableNames[["terms"]])
      tt <- drop.terms(tt, match(bws$variableNames$response, attr(tt, 'term.labels')))
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      emf <- umf
      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }

    q.args <- list(txdat = txdat, tydat = tydat)
    if (has.eval)
      q.args$exdat <- exdat
    q.args$bws <- bws
    tbw <- do.call(npqreg, c(q.args, .npqreg_fit_dots(list(...))))

    tbw$omit <- attr(umf,"na.action")
    tbw$rows.omit <- as.vector(tbw$omit)
    tbw$nobs.omit <- length(tbw$rows.omit)

    tbw$quantile <- .npqreg_napredict_eval(tbw$omit, tbw$quantile)
    tbw$quanterr <- .npqreg_napredict_eval(tbw$omit, tbw$quanterr)

    if(tbw$gradients){
        tbw$quantgrad <- .npqreg_napredict_eval(tbw$omit, tbw$quantgrad)
        tbw$quantgerr <- .npqreg_napredict_eval(tbw$omit, tbw$quantgerr)
    }

    return(tbw)
  }

npqreg.call <-
  function(bws, ...) {
    npqreg(txdat = .np_eval_bws_call_arg(bws, "xdat"),
           tydat = .np_eval_bws_call_arg(bws, "ydat"),
           bws = bws, ...)
  }

npqreg.conbandwidth <-
  function(bws, ...){
    stop("incorrect bandwidth type: expected conditional distribution bandwidths instead of conditional density bandwidths")
  }

npqreg.condbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           tau = 0.5,
           gradients = FALSE,
           tol = 1.490116e-04,
           small = 1.490116e-05, itmax = 10000,
           ...){

    fit.start <- proc.time()[3]
    fit.dots <- .npqreg_fit_dots(list(...))
    if (length(fit.dots))
      stop(sprintf("unused npqreg fit argument '%s'", names(fit.dots)[1L]))
    gradients <- npValidateScalarLogical(gradients, "gradients")
    if (!is.numeric(itmax) || length(itmax) != 1L || is.na(itmax) ||
        !is.finite(itmax) || itmax < 1 || itmax != floor(itmax))
      stop("'itmax' must be a positive integer")
    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol <= 0)
      stop("'tol' must be a positive finite numeric scalar")
    if (!is.numeric(small) || length(small) != 1L || is.na(small) ||
        !is.finite(small) || small <= 0)
      stop("'small' must be a positive finite numeric scalar")
    itmax <- as.integer(itmax)
    tol <- as.double(tol)
    small <- as.double(small)

    no.ex = missing(exdat)

    txdat = toFrame(txdat)
    tydat = toFrame(tydat)

    tau <- .npqreg_validate_tau(tau)
    ntau <- length(tau)
    tau.labels <- .npqreg_tau_labels(tau)

    if (dim(tydat)[2] != 1)
      stop("'tydat' has more than one column")

    if (!no.ex){
      exdat = toFrame(exdat)
      
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
    }

    if (length(bws$xbw) != length(txdat))
      stop("length of bandwidth vector does not match number of columns of 'txdat'")

    if (length(bws$ybw) != 1)
      stop("length of bandwidth vector does not match number of columns of 'tydat'")

    if (any(bws$iyord) || any(bws$iyuno) || coarseclass(tydat[,1]) != "numeric")
      stop("'tydat' is not continuous")

    if ((any(bws$ixcon) &&
         !all(vapply(txdat[, bws$ixcon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$ixord) &&
         !all(vapply(txdat[, bws$ixord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$ixuno) &&
         !all(vapply(txdat[, bws$ixuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'txdat' in type")

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(data.frame(txdat, tydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows,,drop = FALSE]

    if (!no.ex){
      keep.eval <- rep_len(TRUE, nrow(exdat))
      rows.omit <- attr(na.omit(exdat), "na.action")
      if (length(rows.omit) > 0L)
        keep.eval[as.integer(rows.omit)] <- FALSE
      exdat <- exdat[keep.eval,,drop = FALSE]
    }
    
    tnrow = dim(txdat)[1]
    enrow = (if (no.ex) tnrow else dim(exdat)[1])

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati)
    tydat <- adjustLevels(tydat, bws$ydati)
    
    if (!no.ex){
      exdat <- adjustLevels(exdat, bws$xdati)
    }

    ## grab the evaluation data before it is converted to numeric
    if(no.ex){
      txeval <- txdat
    } else {
      txeval <- exdat
    }
    txdat.df <- txdat
    tydat.df <- tydat
    if (!no.ex)
      exdat.df <- exdat

    fit_one_tau <- function(tau_i) {
      yq <- .npqreg_invert_selected_cdf(
        bws = bws,
        xdat = txdat.df,
        ydat = tydat.df,
        exdat = txeval,
        tau = tau_i,
        tol = tol,
        small = small,
        itmax = itmax
      )
      qdelta <- .npqreg_quantile_delta_from_conditional(
        bws = bws,
        xdat = txdat.df,
        ydat = tydat.df,
        exdat = txeval,
        quantile = yq,
        gradients = gradients
      )
      list(
        yq = yq,
        yqerr = qdelta$quanterr,
        yqgrad = if (gradients) qdelta$quantgrad else NA,
        yqgerr = if (gradients) qdelta$quantgerr else NA
      )
    }

    tau.out <- lapply(tau, fit_one_tau)

    if (ntau == 1L) {
      myout <- tau.out[[1L]]
    } else {
      myout <- list(
        yq = do.call(cbind, lapply(tau.out, `[[`, "yq")),
        yqerr = do.call(cbind, lapply(tau.out, `[[`, "yqerr")),
        yqgrad = NA,
        yqgerr = NA
      )
      colnames(myout$yq) <- tau.labels
      colnames(myout$yqerr) <- tau.labels
      if (gradients) {
        p <- ncol(tau.out[[1L]]$yqgrad)
        grad.names <- colnames(tau.out[[1L]]$yqgrad)
        myout$yqgrad <- array(NA_real_,
                              dim = c(enrow, p, ntau),
                              dimnames = list(NULL, grad.names, tau.labels))
        myout$yqgerr <- array(NA_real_,
                              dim = c(enrow, p, ntau),
                              dimnames = list(NULL, grad.names, tau.labels))
        for (j in seq_len(ntau)) {
          myout$yqgrad[, , j] <- tau.out[[j]]$yqgrad
          myout$yqgerr[, , j] <- tau.out[[j]]$yqgerr
        }
      }
    }


    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)

    qregression(bws = bws,
                xeval = txeval,
                tau = tau,
                quantile = myout$yq,
                quanterr = myout$yqerr,
                quantgrad = myout$yqgrad,
                quantgerr = myout$yqgerr,
                ntrain = tnrow,
                trainiseval = no.ex,
                gradients = gradients,
                timing = bws$timing, total.time = total.time,
                optim.time = optim.time, fit.time = fit.elapsed)
  }


npqreg.default <- function(bws, txdat, tydat, nomad = FALSE, ...){
  nomad <- npValidateScalarLogical(nomad, "nomad")
  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  has.explicit.bws <- (!no.bws) && isa(bws, "condbandwidth")
  bws.formula <- (!no.bws) && inherits(bws, "formula")

  if (bws.named && no.txdat && no.tydat && bws.formula) {
    sc$`bws` <- NULL
    sc$formula <- bws
    sc.bw <- sc
    sc.bw[[1]] <- quote(npcdistbw)
    sc.bw <- .npqreg_strip_fit_controls_from_bw_call(sc.bw)
    bws.named <- FALSE
  } else if (bws.formula) {
    sc.bw <- sc
    sc.bw$`bws` <- NULL
    sc.bw[[1]] <- quote(npcdistbw)
    sc.bw <- .npqreg_strip_fit_controls_from_bw_call(sc.bw)
    bws.named <- FALSE
  } else {
    sc.bw <- sc
    sc.bw[[1]] <- quote(npcdistbw)
    sc.bw <- .npqreg_strip_fit_controls_from_bw_call(sc.bw)
  }

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat')
  nstxy <- c('xdat','ydat')
  
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
        "Selecting conditional distribution bandwidth",
        .np_eval_bw_call(sc.bw, caller_env = parent.frame())
      )
    } else {
      .np_eval_bw_call(sc.bw, caller_env = parent.frame())
    }
  } else {
    .np_eval_bw_call(sc.bw, caller_env = parent.frame())
  }

  call.args <- list(bws = tbw)
  if (no.bws) {
    call.args$txdat <- txdat
    call.args$tydat <- tydat
  } else {
    if (txdat.named) call.args$txdat <- txdat
    if (tydat.named) call.args$tydat <- tydat
    if ((!bws.named) && (!txdat.named) && (!no.tydat) && (!tydat.named)) {
      call.args <- c(call.args, list(tydat))
    }
  }
  dots <- list(...)
  if (has.explicit.bws)
    fit.dots <- .npqreg_fit_dots(dots)
  else
    fit.dots <- .npqreg_fit_dots(dots, allow.bandwidth.controls = TRUE)
  do.call(npqreg, c(call.args, fit.dots))
}
