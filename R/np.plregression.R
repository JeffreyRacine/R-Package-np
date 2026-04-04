npplreg <-
  function(bws, ...){
    args <- list(...)

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
  function(bws, data = NULL, newdata = NULL, y.eval = FALSE, ...){
    
    tt <- terms(bws)
    tt.xf <- bws$xterms
    
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf.xf <- tmf <- bws$call[c(1,m)]
    
    tmf[[1]] <- as.name("model.frame")
    tmf.xf[[1]] <- as.name("model.frame")

    tmf.xf[["formula"]] <- tt.xf
    tmf[["formula"]] <- tt

    mf.args <- as.list(tmf)[-1L]
    mf.xf.args <- as.list(tmf.xf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))
    tmf.xf <- do.call(stats::model.frame, mf.xf.args, envir = environment(tt.xf))
    
    response.name <- attr(tmf, "names")[attr(attr(tmf, "terms"), "response")]
    tydat <- model.response(tmf)
    txdat <- tmf.xf
    tzdat <- tmf[, bws$chromoly[[3]], drop = FALSE]

    has.eval <- !is.null(newdata)
    if (has.eval) {
      if (!y.eval){
        tt <- delete.response(tt)
        
        bronze <- lapply(bws$chromoly, paste, collapse = " + ")
        formula.xz <- terms(as.formula(paste(" ~ ",bronze[[2]], " + ",bronze[[3]]),
                                       env = environment(bws$formula)))

        orig.ts <- .np_terms_ts_mask(terms_obj = formula.xz, data = newdata)

        arguments.mfx <- bws$chromoly[[2]]
        arguments.mf <- bws$chromoly[[3]]

        if(all(orig.ts)){
          arguments <- (as.list(attr(formula.xz, "variables"))[-1])
          attr(tt, "predvars") <- bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mf,arguments)),drop = FALSE])
          attr(tt.xf, "predvars") <- bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mfx,arguments)),drop = FALSE])
        }else if(any(orig.ts)){
          arguments <- (as.list(attr(formula.xz, "variables"))[-1])
          arguments.normal <- arguments[which(!orig.ts)]
          arguments.timeseries <- arguments[which(orig.ts)]

          ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
          attr(tt, "predvars") <- bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mf,arguments)),drop = FALSE])
          attr(tt.xf, "predvars") <- bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mfx,arguments)),drop = FALSE])
        }else{
          attr(tt, "predvars") <- attr(tt, "variables")
          attr(tt.xf, "predvars") <- attr(tt.xf, "variables")
        }
          
      }
      
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      emf <- umf
      emf.xf.args <- list(formula = tt.xf, data = newdata)
      emf.xf <- do.call(stats::model.frame, emf.xf.args, envir = parent.frame())
      
      if (y.eval)
        eydat <- model.response(emf)

      exdat <- emf.xf
      ezdat <- emf[, bws$chromoly[[3]], drop = FALSE]
    }

    pl.args <- list(txdat = txdat, tydat = tydat, tzdat = tzdat)
    if (has.eval) {
      pl.args$exdat <- exdat
      pl.args$ezdat <- ezdat
      if (y.eval)
        pl.args$eydat <- eydat
    }
    pl.args$bws <- bws
    ev <- do.call(npplreg, c(pl.args, list(...)))

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

npplreg.call <-
  function(bws, ...) {
    npplreg(txdat = .np_eval_bws_call_arg(bws, "xdat"),
            tydat = .np_eval_bws_call_arg(bws, "ydat"),
            tzdat = .np_eval_bws_call_arg(bws, "zdat"),
            bws = bws, ...)
  }

.np_plreg_fit_progress_targets <- function(xnames) {
  c("y~z", sprintf("%s~z", xnames))
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

.np_plot_plreg_local_fit <-
  function(bws,
           xdat,
           ydat,
           zdat,
           exdat,
           ezdat) {
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

    if (is.factor(ydat)) {
      tmp.ty <- adjustLevels(data.frame(ydat), bws$bw$yzbw$ydati)
      tmp.ty <- (bws$bw$yzbw$ydati$all.dlev[[1L]])[as.integer(tmp.ty)]
    } else {
      tmp.ty <- as.double(ydat)
    }

    local.direct <- isTRUE(identical(bws$type, "generalized_nn"))

    reg_mean <- function(regbw, ytrain, zeval = NULL) {
      args <- list(
        bws = regbw,
        txdat = zdat,
        tydat = ytrain,
        local.mode = local.direct
      )
      if (!is.null(zeval))
        args$exdat <- zeval
      as.vector(.npRmpi_with_local_regression(do.call(.np_regression_direct, args))$mean)
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

    for (j in seq_len(p)) {
      xhat.train <- reg_mean(regbw = bws$bw[[j + 1L]], ytrain = xdat[, j])

      if (is.factor(xdat[1L, j])) {
        tmp.dat <- adjustLevels(xdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati)
        x.num.train <- (bws$bw[[j + 1L]]$ydati$all.dlev[[1L]])[as.integer(tmp.dat[, 1L])]
      } else {
        x.num.train <- as.double(xdat[, j])
      }
      resx[, j] <- x.num.train - xhat.train

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

    model <- lm(resy ~ resx - 1)
    B <- coef(model)

    train.fit <- as.vector(yhat.train + resx %*% B)
    Bvcov <- sum((tmp.ty - train.fit)^2) /
      (ntrain - ncol(xdat) - ncol(zdat)) *
      chol2inv(chol(t(model.matrix(model)) %*% model.matrix(model)))
    Berr <- sqrt(diag(Bvcov))

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

    do.call(plregression, list(
      bws = bws,
      xcoef = B,
      xcoeferr = Berr,
      xcoefvcov = Bvcov,
      evalx = if (no.exz) xdat else exdat,
      evalz = if (no.exz) zdat else ezdat,
      mean = ply,
      ntrain = ntrain,
      trainiseval = no.exz,
      residuals = FALSE,
      xtra = c(RSQ, MSE, MAE, MAPE, CORR, SIGN)
    ))
  }


npplreg.plbandwidth <- 
  function(bws,
           txdat = stop("training data txdat missing"),
           tydat = stop("training data tydat missing"),
           tzdat = stop("training data tzdat missing"),
           exdat, eydat, ezdat, residuals = FALSE, ...){

    fit.start <- proc.time()[3]
    residuals <- npValidateScalarLogical(residuals, "residuals")
    dots <- list(...)
    fit.progress.handoff <- isTRUE(dots$.np_fit_progress_handoff)
    .npRmpi_require_active_slave_pool(where = "npplreg()")
    use.master.fit.progress <- .npRmpi_autodispatch_active() &&
      is.null(bws$degree.search) &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      isTRUE(.np_progress_enabled(domain = "bandwidth"))
    if (!use.master.fit.progress &&
        .npRmpi_autodispatch_active() &&
        is.null(bws$degree.search) &&
        !isTRUE(.npRmpi_autodispatch_called_from_bcast())) {
      result <- .npRmpi_autodispatch_call(match.call(), parent.frame())
      return(.npRmpi_restore_nomad_fit_bws_metadata(result, bws))
    }

    txdat = toFrame(txdat)
    tzdat = toFrame(tzdat)
    
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

    no.exz = missing(exdat)
    no.ey = missing(eydat)

    if (!no.exz){
      exdat = toFrame(exdat)
      ezdat = toFrame(ezdat)
      if (!no.ey)
        eydat = as.double(eydat)

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
    if (is.factor(tydat)){
      tmp.ty <- adjustLevels(data.frame(tydat), bws$bw$yzbw$ydati)
      tmp.ty <- (bws$bw$yzbw$ydati$all.dlev[[1]])[as.integer(tmp.ty)]
    } else {
      tmp.ty <- as.double(tydat)
    }

    if (!no.ey){
      if (is.factor(tydat)){
        tmp.ey <- adjustLevels(data.frame(eydat), bws$bw$yzbw$ydati)
        tmp.ey <- (bws$bw$yzbw$ydati$all.dlev[[1]])[as.integer(tmp.ey)]
      } else {
        tmp.ey <- as.double(eydat)
      }
    }
    
    ## y on z
    mmy = npreg(txdat = tzdat, tydat = tydat, bws = bws$bw$yzbw)

    resy <- tmp.ty - mmy$mean

    if (!no.exz)
      mmy.eval = npreg(txdat = tzdat, tydat = tydat, exdat = ezdat, bws = bws$bw$yzbw)

    
    ## x on z
    nrow = nrow(txdat)
    nrow.eval = (if (no.exz) 0 else nrow(exdat))
    ncol = ncol(txdat)
    B = double(ncol)
    resx = matrix(data = 0, nrow = nrow, ncol = ncol)
    resx.eval = matrix(data = 0, nrow = nrow.eval, ncol = ncol)
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
      mm = npreg(txdat=tzdat, tydat=txdat[,i], bws = bws$bw[[i+1]])

      if (is.factor(txdat[1,i])){
        tmp.dat <- adjustLevels(txdat[,i, drop=FALSE], bws$bw[[i+1]]$ydati)
        resx[,i] <- (bws$bw[[i+1]]$ydati$all.dlev[[1]])[as.integer(tmp.dat[,1])] - mm$mean
      } else {
        resx[,i] <- txdat[,i] - mm$mean
      }

      if(!no.exz) {
        mm = npreg(txdat=tzdat, tydat=txdat[,i], exdat=ezdat, bws = bws$bw[[i+1]])

        if (is.factor(txdat[1,i])){
          tmp.dat <- adjustLevels(exdat[,i, drop=FALSE], bws$bw[[i+1]]$ydati)
          resx.eval[,i] <- (bws$bw[[i+1]]$ydati$all.dlev[[1]])[as.integer(tmp.dat[,1])] - mm$mean
        } else {
          resx.eval[,i] <- exdat[,i] - mm$mean
        }
      }

      fit.progress <- .np_progress_step(
        fit.progress,
        done = i + 1L,
        detail = fit.progress.targets[i + 1L]
      )
    }

    B = coef((model = lm(resy ~ resx - 1)))

    ## computes the standard errors of B using the model matrix
    ## and the MSE of the training data predictions

    Bvcov = sum((tmp.ty-(mmy$mean + resx  %*% B))^2)/
      (dim(txdat)[1]-dim(txdat)[2]-dim(tzdat)[2])*
      chol2inv(chol(t(model.matrix(model))%*%model.matrix(model)))

    Berr = sqrt(diag(Bvcov))

    train.ply =  mmy$mean + resx %*% B
    ply = if (no.exz) train.ply else mmy.eval$mean + resx.eval %*% B

    if (!no.ey) {
      RSQ = RSQfunc(tmp.ey, ply)
      MSE = MSEfunc(tmp.ey, ply)
      MAE = MAEfunc(tmp.ey, ply)
      MAPE = MAPEfunc(tmp.ey, ply)
      CORR = CORRfunc(tmp.ey, ply)
      SIGN = SIGNfunc(tmp.ey, ply)

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
      evalx = if (no.exz) txdat else exdat,
      evalz = if (no.exz) tzdat else ezdat,
      mean = ply,
      ntrain = nrow,
      trainiseval = no.exz,
      residuals = residuals,
      xtra = c(RSQ, MSE, MAE, MAPE, CORR, SIGN)
    )
    if (residuals)
      ev.args$resid <- tmp.ty - train.ply
    ev <- do.call(plregression, ev.args)

    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)
    ev$timing <- bws$timing
    ev$total.time <- total.time
    ev$optim.time <- optim.time
    ev$fit.time <- fit.elapsed

    
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    fit.progress <- .np_progress_end(fit.progress)
    fit.progress.active <- FALSE
    return(ev)
  }


npplreg.default <- function(bws, txdat, tydat, tzdat, nomad = FALSE, ...) {
  .npRmpi_require_active_slave_pool(where = "npplreg()")
  explicit.plbandwidth <- (!missing(bws)) && inherits(bws, "plbandwidth")
  nomad <- npValidateScalarLogical(nomad, "nomad")
  degree.select.value <- if (isTRUE(nomad)) {
    "coordinate"
  } else if ("degree.select" %in% names(list(...))) {
    match.arg(list(...)$degree.select, c("manual", "coordinate", "exhaustive"))
  } else {
    "manual"
  }
  if (.npRmpi_autodispatch_active() &&
      (explicit.plbandwidth || identical(degree.select.value, "manual")) &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()))
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
  has.explicit.bws <- (!no.bws) && isa(bws, "plbandwidth")

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  if(tydat.named)
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
    
  use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
    sc.bw,
    caller_env = parent.frame()
  )

  tbw <- if (!has.explicit.bws) {
    if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting partially linear regression bandwidth",
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
  do.call(npplreg, c(call.args, list(...)))
}
