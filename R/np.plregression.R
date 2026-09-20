npplreg <-
  function(bws, ...){
    mc <- match.call(expand.dots = FALSE)
    .np_validate_public_dots(mc[["..."]], "npplreg")
    args <- .np_formula_dispatch_args(
      NULL, substitute(list(...))[-1L], environment())

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npplreg",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npplreg",bws$call)
        else if (!is.call(bws))
          UseMethod("npplreg",bws)
        else
          UseMethod("npplreg",NULL)
      } else {
        UseMethod("npplreg", NULL)
      }
    } else {
      UseMethod("npplreg", NULL)
    }
  }

npplreg.formula <-
  function(bws, data = NULL, newdata = NULL, y.eval = FALSE, ..., se = FALSE){
    se <- npValidateScalarLogical(se, "se")
    dots <- list(...)
    frame.state <- dots[[".np.formula.state", exact = TRUE]]
    dots$.np.formula.state <- NULL
    frame <- if (is.null(frame.state)) .np_plreg_formula_frame(bws, data, overrides = dots)
             else .np_formula_frame_take(frame.state)
    if ("tydat" %in% names(dots)) {
      frame <- .np_formula_replace_response(frame, bws, dots[["tydat"]],
                                             dots, data.override = !is.null(data))
      dots$tydat <- NULL
      dots$na.action <- NULL
    }
    bws <- .np_bws_retain_fit_frame(bws, frame)
    roles <- .np_plreg_formula_split(frame, bws$terms, bws$xterms)
    response.name <- names(roles$yz)[attr(bws$terms, "response")]
    pl.args <- list(txdat = roles$x, tydat = model.response(roles$yz),
                    tzdat = roles$yz[, .np_formula_term_names(bws$chromoly[[3L]]), drop = FALSE],
                    se = se)
    umf <- frame
    native.eval <- !is.null(dots[["exdat", exact = TRUE]]) ||
      !is.null(dots[["ezdat", exact = TRUE]])
    has.eval <- !is.null(newdata) && !native.eval
    if (has.eval) {
      tt <- attr(frame, "terms")
      response.eval <- y.eval && is.null(dots[["eydat", exact = TRUE]])
      npValidateNewdataFormula(newdata, tt, include.response = response.eval)
      yzterms <- bws$terms
      if (!response.eval) {
        tt <- delete.response(tt)
        yzterms <- delete.response(yzterms)
      }
      umf <- .np_formula_model_frame(tt, data = newdata)
      evaluation <- .np_plreg_formula_split(umf, yzterms, bws$xterms)
      pl.args$exdat <- evaluation$x
      pl.args$ezdat <- evaluation$yz[, .np_formula_term_names(bws$chromoly[[3L]]), drop = FALSE]
      if (response.eval)
        pl.args$eydat <- model.response(evaluation$yz)
    }
    pl.args$bws <- bws
    ev <- do.call(npplreg, c(pl.args, dots))

    if (length(response.name) == 1L && !is.na(response.name) && nzchar(response.name)) {
      if (!is.null(ev$bws))
        ev$bws$ynames <- response.name
    }

    # Residuals always index the training sample, independently of evaluation.
    if (ev$residuals)
      ev$resid <- naresid(attr(frame, "na.action"), ev$resid)
    # Native evaluation already owns its row omissions and padding.
    if (native.eval)
      return(ev)
    ev$omit <- attr(umf, "na.action")
    ev$rows.omit <- as.vector(ev$omit)
    ev$nobs.omit <- length(ev$rows.omit)
    ev$mean <- napredict(ev$omit, ev$mean)
    ev$merr <- napredict(ev$omit, ev$merr)
    ev
  }

npplreg.call <-
  function(bws, ...) {
    do.call(npplreg, .np_retained_training_args(
      bws, c(txdat = "xdat", tydat = "ydat", tzdat = "zdat"), list(...)))
  }

.np_plreg_fit_progress_targets <- function(xnames) {
  c("E[y|z]", sprintf("E[%s|z]", xnames), "final partially linear solve")
}

.np_plreg_fit_progress_begin <- function(xnames, handoff = FALSE) {
  state <- .np_progress_begin(
    "Fitting partially linear regression",
    total = length(.np_plreg_fit_progress_targets(xnames)),
    surface = "bandwidth"
  )

  if (isTRUE(handoff)) {
    state <- .np_progress_show_now(
      state = state,
      done = 0L,
      detail = paste("starting", .np_plreg_fit_progress_targets(xnames)[1L])
    )
  }

  state
}

.np_plreg_numeric_response <- function(y, ydati) {
  if (is.factor(y)) {
    yy <- adjustLevels(data.frame(y), ydati)
    return(ydati$all.dlev[[1L]][as.integer(yy[, 1L])])
  }

  as.double(y)
}

.np_plreg_residual_formation_error <- function(x, xhat) {
  .Machine$double.eps * sqrt(length(x)) * (max(abs(x)) + max(abs(xhat)))
}

.np_plreg_check_residualized_rank <- function(qrX, p, where,
                                              formation.error = NULL) {
  deficient <- qrX[["rank"]] < p
  if (!deficient) {
    # Inspect only the existing small factor; do not refactor or change the
    # pivot/solve. Column scaling makes the guard independent of other units.
    R <- qrX[["qr"]][seq_len(p), seq_len(p), drop = FALSE]
    R[lower.tri(R)] <- 0
    scale <- vapply(seq_len(p), function(j) norm(R[, j, drop = FALSE], "F"),
                    numeric(1L))
    precision <- min(max(dim(qrX[["qr"]])) * .Machine$double.eps,
                     sqrt(.Machine$double.eps))
    tolerance <- precision * scale
    if (!is.null(formation.error))
      tolerance <- tolerance + formation.error[qrX[["pivot"]]]
    deficient <- any(abs(diag(R)) <= tolerance)
  }
  if (deficient) {
    stop(sprintf(
      "%s: residualized linear regressors are rank deficient after smoothing on z; the parametric component is not identified",
      where
    ), call. = FALSE)
  }
}

.np_plreg_linear_solve <- function(resy,
                                   resx,
                                   response,
                                   yhat.train,
                                   zdim,
                                   where = "npplreg",
                                   se = TRUE,
                                   formation.error = NULL) {
  X <- as.matrix(resx)
  y <- as.double(resy)
  response <- as.double(response)
  yhat.train <- as.double(yhat.train)
  p <- ncol(X)
  qrX <- qr(X, tol = .Machine$double.eps)

  .np_plreg_check_residualized_rank(qrX = qrX, p = p, where = where,
                                    formation.error = formation.error)

  beta <- as.double(qr.coef(qrX, y))
  linear.fit <- as.vector(X %*% beta)
  train.fit <- as.vector(yhat.train + linear.fit)

  vcov <- stderr <- NULL
  if (se) {
    R <- qr.R(qrX)[seq_len(p), seq_len(p), drop = FALSE]
    XtX.inv.pivoted <- chol2inv(R)
    XtX.inv <- matrix(0.0, nrow = p, ncol = p)
    XtX.inv[qrX$pivot, qrX$pivot] <- XtX.inv.pivoted

    # Factor the heteroskedastic sandwich as a crossproduct of coefficient
    # influences, preserving the existing residuals and finite-sample factor.
    scores <- (X %*% XtX.inv) * (response - train.fit)
    vcov <- crossprod(scores) * (nrow(X) / (nrow(X) - p - as.integer(zdim)))
    stderr <- sqrt(diag(vcov))
  }

  list(
    coef = beta,
    vcov = vcov,
    se = stderr,
    train.fit = train.fit,
    qr = qrX
  )
}

.np_plreg_pad_eval_vector <- function(x, keep.eval) {
  out <- rep(NA_real_, length(keep.eval))
  out[keep.eval] <- as.double(x)
  out
}

.np_plot_plreg_local_fit <-
  function(bws,
           xdat,
           ydat,
           zdat,
           exdat,
           ezdat,
           .np.empty.report = NULL,
           se = TRUE) {
    se <- npValidateScalarLogical(se, "se")
    activity <- .np_plot_activity_begin("Computing partially linear plot fit")
    on.exit(.np_plot_activity_end(activity), add = TRUE)

    xdat <- toFrame(xdat)
    zdat <- toFrame(zdat)

    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- attr(na.omit(data.frame(xdat, ydat, zdat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    xdat <- xdat[keep.rows, , drop = FALSE]
    ydat <- ydat[keep.rows]
    zdat <- zdat[keep.rows, , drop = FALSE]

    no.exz <- missing(exdat)
    if (!no.exz) {
      exdat <- toFrame(exdat)
      ezdat <- toFrame(ezdat)

      keep.eval <- rep_len(TRUE, nrow(exdat))
      rows.omit <- attr(na.omit(data.frame(exdat, ezdat)), "na.action")
      if (length(rows.omit) > 0L)
        keep.eval[as.integer(rows.omit)] <- FALSE

      if (!any(keep.eval))
        stop("Evaluation data has no rows without NAs")

      exdat <- exdat[keep.eval, , drop = FALSE]
      ezdat <- ezdat[keep.eval, , drop = FALSE]
    }

    tmp.ty <- .np_plreg_numeric_response(ydat, bws$bw$yzbw$ydati)

    local.direct <- isTRUE(identical(bws$type, "generalized_nn"))
    cat.profile.cache <- new.env(parent = emptyenv())

    reg_mean_cat_cached <- function(regbw, ytrain, zeval = NULL) {
      if (!npUseCategoricalCompress(ncon = regbw$ncon,
                                    ncat = regbw$nuno + regbw$nord))
        return(NULL)
      if (!identical(regbw$type, "fixed"))
        return(NULL)
      if (!isTRUE(regbw$ncon == 0L) || (regbw$nuno + regbw$nord) < 1L)
        return(NULL)

      regtype <- if (is.null(regbw$regtype)) "lc" else as.character(regbw$regtype)
      if (!(identical(regtype, "lc") || identical(regtype, "lp")))
        return(NULL)
      if (identical(regtype, "lp") && length(regbw$degree) && any(regbw$degree > 0L))
        return(NULL)

      if (ncol(zdat) != length(regbw$bw) || length(ytrain) != nrow(zdat))
        return(NULL)
      if (!is.null(zeval) && ncol(zeval) != ncol(zdat))
        return(NULL)
      if (!is.null(zeval) && !(zdat %~% zeval))
        return(NULL)
      if (!all(vapply(zdat, function(z) is.factor(z) || is.ordered(z), logical(1))))
        return(NULL)
      if (!is.null(zeval) &&
          !all(vapply(zeval, function(z) is.factor(z) || is.ordered(z), logical(1))))
        return(NULL)

      if (is.null(cat.profile.cache$xdati)) {
        tx.adj <- adjustLevels(zdat, regbw$xdati)
        eval.plot <- if (no.exz) tx.adj else adjustLevels(ezdat, regbw$xdati, allowNewCells = TRUE)

        if (!all(vapply(seq_along(tx.adj), function(j) {
          identical(is.ordered(tx.adj[[j]]), is.ordered(eval.plot[[j]])) &&
            identical(levels(tx.adj[[j]]), levels(eval.plot[[j]]))
        }, logical(1))))
          return(NULL)

        train.codes <- .np_cat_profile_code_matrix(tx.adj)
        if (anyNA(train.codes))
          return(NULL)
        train.keys <- .np_cat_profile_keys(train.codes)
        train.profile.keys <- unique(train.keys)
        train.id <- match(train.keys, train.profile.keys)
        train.rep <- match(train.profile.keys, train.keys)
        G <- length(train.profile.keys)

        eval_cache <- function(dat) {
          codes <- .np_cat_profile_code_matrix(dat)
          if (anyNA(codes))
            return(NULL)
          keys <- .np_cat_profile_keys(codes)
          profile.keys <- unique(keys)
          rep.idx <- match(profile.keys, keys)
          list(
            id = match(keys, profile.keys),
            profile.codes = codes[rep.idx, , drop = FALSE],
            profile.keys = profile.keys
          )
        }

        eval.train.cache <- eval_cache(tx.adj)
        eval.plot.cache <- eval_cache(eval.plot)
        if (is.null(eval.train.cache) || is.null(eval.plot.cache))
          return(NULL)

        cat.profile.cache$xdati <- regbw$xdati
        cat.profile.cache$train.id <- train.id
        cat.profile.cache$train.profile.codes <- train.codes[train.rep, , drop = FALSE]
        cat.profile.cache$train.profile.dat <- tx.adj[train.rep, , drop = FALSE]
        cat.profile.cache$G <- G
        cat.profile.cache$counts <- as.double(tabulate(train.id, nbins = G))
        cat.profile.cache$eval.train <- eval.train.cache
        cat.profile.cache$eval.plot <- eval.plot.cache
      } else if (!identical(cat.profile.cache$xdati, regbw$xdati)) {
        return(NULL)
      }

      if (is.factor(ytrain)) {
        ytrain <- adjustLevels(data.frame(ytrain), regbw$ydati)[, 1L]
        ytrain <- (regbw$ydati$all.dlev[[1L]])[as.integer(ytrain)]
      } else {
        ytrain <- as.double(ytrain)
      }
      if (anyNA(ytrain))
        return(NULL)

      eval.cache <- if (is.null(zeval)) {
        cat.profile.cache$eval.train
      } else {
        cat.profile.cache$eval.plot
      }
      L.eval <- tryCatch(
        .np_regression_cat_profile_kernel_matrix(
          eval.codes = eval.cache$profile.codes,
          train.codes = cat.profile.cache$train.profile.codes,
          xdat = cat.profile.cache$train.profile.dat,
          bws = regbw
        ),
        error = function(e) NULL
      )
      if (is.null(L.eval))
        return(NULL)

      sums <- as.vector(.np_cat_profile_rowsum(ytrain,
                                               cat.profile.cache$train.id,
                                               cat.profile.cache$G))
      den <- as.vector(L.eval %*% cat.profile.cache$counts)
      if (any(!is.finite(den)) || any(!(abs(den) > .Machine$double.xmin)))
        return(NULL)

      profile.mean <- as.vector(L.eval %*% sums / den)
      profile.mean[eval.cache$id]
    }

    empty.state <- new.env(hash = FALSE, parent = emptyenv())
    empty.state$rows <- NULL
    reg_mean <- function(regbw, ytrain, zeval = NULL) {
      out <- reg_mean_cat_cached(regbw = regbw, ytrain = ytrain, zeval = zeval)
      if (is.null(out)) {
        out <- .np_regression_cat_profile_mean(
          bws = regbw,
          txdat = zdat,
          tydat = ytrain,
          exdat = zeval
        )
      }
      if (!is.null(out))
        return(as.vector(out))

      args <- list(
        bws = regbw,
        txdat = zdat,
        tydat = ytrain,
        local.mode = local.direct
      )
      if (!is.null(zeval))
        args$exdat <- zeval
      args$allow.empty.rows <- !is.null(zeval)
      fit <- do.call(.np_regression_direct, args)
      empty.state$rows <- .npreg_merge_empty_rows(empty.state$rows,
        attr(fit, ".np.empty.rows", exact = TRUE))
      as.vector(fit$mean)
    }

    yhat.train <- reg_mean(regbw = bws$bw$yzbw, ytrain = ydat)
    resy <- tmp.ty - yhat.train

    if (!no.exz)
      yhat.eval <- reg_mean(regbw = bws$bw$yzbw, ytrain = ydat, zeval = ezdat)

    ntrain <- nrow(xdat)
    neval <- if (no.exz) ntrain else nrow(exdat)
    p <- ncol(xdat)
    resx <- matrix(0.0, nrow = ntrain, ncol = p)
    resx.eval <- matrix(0.0, nrow = neval, ncol = p)
    formation.error <- numeric(p)

    for (j in seq_len(p)) {
      xhat.train <- reg_mean(regbw = bws$bw[[j + 1L]], ytrain = xdat[, j])

      if (is.factor(xdat[1L, j])) {
        tmp.dat <- adjustLevels(xdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati)
        x.num.train <- (bws$bw[[j + 1L]]$ydati$all.dlev[[1L]])[as.integer(tmp.dat[, 1L])]
      } else {
        x.num.train <- as.double(xdat[, j])
      }
      resx[, j] <- x.num.train - xhat.train
      formation.error[j] <- .np_plreg_residual_formation_error(x.num.train,
                                                               xhat.train)

      if (!no.exz) {
        xhat.eval <- reg_mean(regbw = bws$bw[[j + 1L]], ytrain = xdat[, j], zeval = ezdat)
        if (is.factor(xdat[1L, j])) {
          tmp.dat <- adjustLevels(exdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati, allowNewCells = TRUE)
          x.num.eval <- (bws$bw[[j + 1L]]$ydati$all.dlev[[1L]])[as.integer(tmp.dat[, 1L])]
        } else {
          x.num.eval <- as.double(exdat[, j])
        }
        resx.eval[, j] <- x.num.eval - xhat.eval
      }
    }

    solved <- .np_plreg_linear_solve(
      resy = resy,
      resx = resx,
      response = tmp.ty,
      yhat.train = yhat.train,
      zdim = ncol(zdat),
      where = ".np_plot_plreg_local_fit",
      se = se,
      formation.error = formation.error
    )
    B <- solved$coef
    train.fit <- solved$train.fit
    Bvcov <- solved$vcov
    Berr <- solved$se

    RSQ <- RSQfunc(tmp.ty, train.fit)
    MSE <- MSEfunc(tmp.ty, train.fit)
    MAE <- MAEfunc(tmp.ty, train.fit)
    MAPE <- MAPEfunc(tmp.ty, train.fit)
    CORR <- CORRfunc(tmp.ty, train.fit)
    SIGN <- SIGNfunc(tmp.ty, train.fit)

    ply <- if (no.exz) {
      train.fit
    } else {
      as.vector(yhat.eval + resx.eval %*% B)
    }

    fit <- do.call(plregression, list(
      bws = bws,
      xcoef = B,
      xcoeferr = Berr,
      xcoefvcov = Bvcov,
      se = se,
      evalx = if (no.exz) xdat else exdat,
      evalz = if (no.exz) zdat else ezdat,
      mean = ply,
      ntrain = ntrain,
      trainiseval = no.exz,
      residuals = FALSE,
      xtra = c(RSQ, MSE, MAE, MAPE, CORR, SIGN)
    ))
    .npreg_publish_plot_rows(fit, empty.state$rows, report = .np.empty.report,
      omitted = if(no.exz) integer(0) else which(!keep.eval),
      row.labels = if(no.exz) NULL else row.names(ezdat))
  }


npplreg.plbandwidth <- 
  function(bws,
           txdat = stop("training data txdat missing"),
           tydat = stop("training data tydat missing"),
           tzdat = stop("training data tzdat missing"),
           exdat, eydat, ezdat, residuals = FALSE, ..., se = FALSE){

    fit.start <- proc.time()[3]
    residuals <- npValidateScalarLogical(residuals, "residuals")
    se <- npValidateScalarLogical(se, "se")
    dots <- list(...)
    fit.progress.handoff <- isTRUE(dots$.np_fit_progress_handoff)

    txdat = toFrame(txdat)
    tzdat = toFrame(tzdat)
    .np_require_paired_rows(txdat, tydat, "txdat", "tydat")
    .np_require_paired_rows(txdat, tzdat, "txdat", "tzdat")
    bws <- .np_bws_retain_native_training(bws, xdat = txdat, ydat = tydat, zdat = tzdat)
    
    ## catch and destroy NA's, part 1
    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(data.frame(txdat,tydat,tzdat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows]
    tzdat <- tzdat[keep.rows,,drop = FALSE]

    native.newdata <- dots[["newdata", exact = TRUE]]
    if (missing(exdat) && missing(ezdat) && !is.null(native.newdata)) {
      native.eval <- .np_native_newdata_parts(
        native.newdata, list(exdat = bws$xnames, ezdat = bws$znames), "npplreg")
      exdat <- native.eval$exdat
      ezdat <- native.eval$ezdat
    }

    no.exz = missing(exdat)
    no.ey = missing(eydat)

    if (!no.exz){
      exdat = toFrame(exdat)
      ezdat = toFrame(ezdat)
      .np_require_paired_rows(exdat, ezdat, "exdat", "ezdat")
      if (!no.ey)
        .np_require_paired_rows(exdat, eydat, "exdat", "eydat")
      exdat.full <- exdat
      ezdat.full <- ezdat

      ## c& d NA's, part 2

      keep.eval <- rep_len(TRUE, nrow(exdat))
      eval.df <- data.frame(exdat, ezdat)
      if (!no.ey)
        eval.df <- data.frame(eval.df, eydat)
      rows.omit <- attr(na.omit(eval.df), "na.action")
      if (length(rows.omit) > 0L)
        keep.eval[as.integer(rows.omit)] <- FALSE

      if (!any(keep.eval))
        stop("Evaluation data has no rows without NAs")

      exdat <- exdat[keep.eval,,drop = FALSE]
      if (!no.ey)
        eydat <- eydat[keep.eval]
      ezdat <- ezdat[keep.eval,,drop = FALSE]
    }

    ## tmp.ty and tmp.ey are the numeric representations of tydat and eydat
    tmp.ty <- .np_plreg_numeric_response(tydat, bws$bw$yzbw$ydati)

    if (!no.ey)
      tmp.ey <- .np_plreg_numeric_response(eydat, bws$bw$yzbw$ydati)
    
    empty.state <- new.env(hash = FALSE, parent = emptyenv())
    empty.state$rows <- NULL
    reg_mean <- function(regbw, ytrain, zeval = NULL) {
      out <- .np_regression_cat_profile_mean(
        bws = regbw,
        txdat = tzdat,
        tydat = ytrain,
        exdat = zeval
      )
      if (!is.null(out))
        return(as.vector(out))

      args <- list(txdat = tzdat, tydat = ytrain, bws = regbw)
      if (!is.null(zeval))
        args$exdat <- zeval
      args$.np.defer.empty.rows <- TRUE
      fit <- do.call(npreg, args)
      empty.state$rows <- .npreg_merge_empty_rows(
        empty.state$rows, attr(fit, ".np.empty.rows", exact = TRUE))
      as.vector(fitted(fit))
    }

    ## y on z
    yhat.train <- reg_mean(regbw = bws$bw$yzbw, ytrain = tydat)

    resy <- tmp.ty - yhat.train

    if (!no.exz)
      yhat.eval <- reg_mean(regbw = bws$bw$yzbw, ytrain = tydat, zeval = ezdat)

    
    ## x on z
    nrow = nrow(txdat)
    nrow.eval = (if (no.exz) 0 else nrow(exdat))
    ncol = ncol(txdat)
    B = double(ncol)
    resx = matrix(data = 0, nrow = nrow, ncol = ncol)
    resx.eval = matrix(data = 0, nrow = nrow.eval, ncol = ncol)
    formation.error <- numeric(ncol)
    fit.progress.targets <- .np_plreg_fit_progress_targets(names(txdat))
    fit.progress <- .np_plreg_fit_progress_begin(
      xnames = names(txdat),
      handoff = fit.progress.handoff
    )
    fit.progress.active <- TRUE
    on.exit({
      if (isTRUE(fit.progress.active))
        .np_progress_abort(fit.progress)
    }, add = TRUE)
    fit.progress <- .np_progress_step(
      fit.progress,
      done = 1L,
      detail = fit.progress.targets[1L]
    )

    for (i in seq_len(ncol)) {
      xhat.train <- reg_mean(regbw = bws$bw[[i+1]], ytrain = txdat[, i])

      if (is.factor(txdat[1,i])){
        tmp.dat <- adjustLevels(txdat[,i, drop=FALSE], bws$bw[[i+1]]$ydati)
        x.num.train <- (bws$bw[[i+1]]$ydati$all.dlev[[1]])[as.integer(tmp.dat[,1])]
        resx[,i] <- x.num.train - xhat.train
      } else {
        x.num.train <- txdat[,i]
        resx[,i] <- x.num.train - xhat.train
      }
      formation.error[i] <- .np_plreg_residual_formation_error(x.num.train,
                                                               xhat.train)

      if(!no.exz) {
        xhat.eval <- reg_mean(regbw = bws$bw[[i+1]], ytrain = txdat[, i], zeval = ezdat)

        if (is.factor(txdat[1,i])){
          tmp.dat <- adjustLevels(exdat[,i, drop=FALSE], bws$bw[[i+1]]$ydati)
          resx.eval[,i] <- (bws$bw[[i+1]]$ydati$all.dlev[[1]])[as.integer(tmp.dat[,1])] - xhat.eval
        } else {
          resx.eval[,i] <- exdat[,i] - xhat.eval
        }
      }

      fit.progress <- .np_progress_step(
        fit.progress,
        done = i + 1L,
        detail = fit.progress.targets[i + 1L]
      )
    }

    fit.progress <- .np_progress_step(
      fit.progress,
      done = ncol + 2L,
      detail = fit.progress.targets[ncol + 2L]
    )

    solved <- .np_plreg_linear_solve(
      resy = resy,
      resx = resx,
      response = tmp.ty,
      yhat.train = yhat.train,
      zdim = dim(tzdat)[2],
      where = "npplreg",
      se = se,
      formation.error = formation.error
    )
    B <- solved$coef
    Bvcov <- solved$vcov
    Berr <- solved$se

    train.ply <- solved$train.fit
    ply.complete <- if (no.exz) train.ply else as.vector(yhat.eval + resx.eval %*% B)
    ply <- if (no.exz) train.ply else .np_plreg_pad_eval_vector(ply.complete, keep.eval)

    if (!no.ey) {
      RSQ = RSQfunc(tmp.ey, ply.complete)
      MSE = MSEfunc(tmp.ey, ply.complete)
      MAE = MAEfunc(tmp.ey, ply.complete)
      MAPE = MAPEfunc(tmp.ey, ply.complete)
      CORR = CORRfunc(tmp.ey, ply.complete)
      SIGN = SIGNfunc(tmp.ey, ply.complete)

    } else if (!no.exz) {
      RSQ = MSE = MAE = MAPE = CORR = SIGN = NA_real_
    } else {
      RSQ = RSQfunc(tmp.ty, train.ply)
      MSE = MSEfunc(tmp.ty, train.ply)
      MAE = MAEfunc(tmp.ty, train.ply)
      MAPE = MAPEfunc(tmp.ty, train.ply)
      CORR = CORRfunc(tmp.ty, train.ply)
      SIGN = SIGNfunc(tmp.ty, train.ply)
    }

    ev.args <- list(
      bws = bws,
      xcoef = B,
      xcoeferr = Berr,
      xcoefvcov = Bvcov,
      se = se,
      evalx = if (no.exz) txdat else exdat.full,
      evalz = if (no.exz) tzdat else ezdat.full,
      mean = ply,
      ntrain = nrow,
      trainiseval = no.exz,
      residuals = residuals,
      xtra = c(RSQ, MSE, MAE, MAPE, CORR, SIGN)
    )
    if (residuals)
      ev.args$resid <- tmp.ty - train.ply
    ev <- do.call(plregression, ev.args)
    if (!no.exz) {
      ev$eval.keep <- keep.eval
      ev$eval.rows.omit <- which(!keep.eval)
    }

    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)
    ev$timing <- bws$timing
    ev$total.time <- total.time
    ev$optim.time <- optim.time
    ev$fit.time <- fit.elapsed
    ev$nomad.time <- if (!is.null(bws$nomad.time) && is.finite(bws$nomad.time)) as.double(bws$nomad.time) else NA_real_
    ev$powell.time <- if (!is.null(bws$powell.time) && is.finite(bws$powell.time)) as.double(bws$powell.time) else NA_real_

    
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    fit.progress <- .np_progress_end(fit.progress)
    fit.progress.active <- FALSE
    return(.npreg_finish_empty_rows(ev, empty.state$rows,
      omitted = if(no.exz) integer(0) else which(!keep.eval),
      defer = isTRUE(list(...)[[".np.defer.empty.rows", exact = TRUE]]),
      owner = "npplreg", row.labels = if(no.exz) NULL else row.names(exdat)))
  }


npplreg.default <- function(bws, txdat, tydat, tzdat, nomad = FALSE, ..., se = FALSE) {
  sc <- .np_formula_default_call(sys.call(), sys.function(), parent.frame())
  sc.names <- names(sc)
  nomad <- npValidateNomadControl(nomad, "nomad")
  se <- npValidateScalarLogical(se, "se")

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
  has.explicit.bws <- (!no.bws) && isa(bws, "plbandwidth")

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tzdat.named)
    tzdat <- toFrame(tzdat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npplregbw)

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
  sc.bw <- .np_public_dots_filter_call(sc.bw, "npplregbw")
  formula.input <- .np_formula_dispatch_args(
    NULL, substitute(list(...))[-1L], environment())[["formula", exact = TRUE]]
  frame.state <- if (!has.explicit.bws &&
      (bws.formula || inherits(formula.input, "formula") ||
       (!no.txdat && inherits(txdat, "formula"))) &&
      no.tydat && no.tzdat &&
      (no.txdat || inherits(txdat, "formula")))
    new.env(parent = emptyenv()) else NULL
  if (!is.null(frame.state)) {
    sc.bw$.np.formula.state <- frame.state
    on.exit(rm(list = ls(frame.state, all.names = TRUE), envir = frame.state), add = TRUE)
  }
    
  use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
    sc.bw,
    caller_env = parent.frame()
  )

  tbw <- if (!has.explicit.bws) {
    if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting partially linear regression bandwidth",
        .np_eval_bw_call(sc.bw, caller_env = parent.frame(),
                        native.map = c(bws = "bws", xdat = "txdat", ydat = "tydat", zdat = "tzdat"), native.frame = environment(),
                        formula.value = .np_formula_value(formula.input, bws, txdat))
      )
    } else {
        .np_eval_bw_call(sc.bw, caller_env = parent.frame(),
                        native.map = c(bws = "bws", xdat = "txdat", ydat = "tydat", zdat = "tzdat"), native.frame = environment(),
                        formula.value = .np_formula_value(formula.input, bws, txdat))
    }
  } else {
        .np_eval_bw_call(sc.bw, caller_env = parent.frame(),
                        native.map = c(bws = "bws", xdat = "txdat", ydat = "tydat", zdat = "tzdat"), native.frame = environment(),
                        formula.value = .np_formula_value(formula.input, bws, txdat))
  }
  
  call.args <- list(bws = tbw, se = se)
  if (!is.null(frame.state)) {
    call.args$.np.formula.state <- frame.state
  } else if (no.bws) {
    call.args$txdat <- txdat
    call.args$tydat <- tydat
    call.args$tzdat <- tzdat
  } else {
    if (txdat.named) call.args$txdat <- txdat
    if (tydat.named) call.args$tydat <- tydat
    if (tzdat.named) call.args$tzdat <- tzdat
    if ((!bws.named) && (!txdat.named) && (!no.tzdat) && (!tzdat.named)) {
      call.args <- c(call.args, list(tzdat))
    }
  }
  if (!has.explicit.bws)
    call.args$.np_fit_progress_handoff <- TRUE
  fit.dots <- if (!is.null(tbw[["formula", exact = TRUE]]))
    .np_formula_dispatch_args(NULL, substitute(list(...))[-1L], environment())
  else list(...)
  do.call(npplreg, c(call.args, fit.dots))
}
