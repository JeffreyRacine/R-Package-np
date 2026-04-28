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

.npqreg.fit.control.names <- c("data", "newdata", "tau", "gradients", "tol", "small", "itmax")
.npqreg.removed.solver.controls <- c("ftol",
                                     "lbc.dir", "dfc.dir", "cfac.dir", "initc.dir",
                                     "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir")

.npqreg_fit_dots <- function(dots, allow.bandwidth.controls = FALSE) {
  dot.names <- names(dots)
  if (is.null(dot.names))
    return(dots)

  stale <- intersect(dot.names[nzchar(dot.names)], .npqreg.removed.solver.controls)
  if (length(stale) && !allow.bandwidth.controls) {
    stop(sprintf(
      "'%s' %s no longer accepted by npqreg; the canonical one-dimensional quantile extractor is controlled by 'tol', 'small', and 'itmax'",
      paste(stale, collapse = "', '"),
      if (length(stale) == 1L) "is" else "are"
    ))
  }

  keep <- (!nzchar(dot.names)) | (dot.names %in% .npqreg.fit.control.names)
  dots[keep]
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

    tbw$quantile <- napredict(tbw$omit, tbw$quantile)
    tbw$quanterr <- napredict(tbw$omit, tbw$quanterr)

    if(tbw$gradients){
        tbw$quantgrad <- napredict(tbw$omit, tbw$quantgrad)
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

.npRmpi_npqreg_should_localize <- function(bws) {
  isa(bws, "condbandwidth")
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
    .npRmpi_require_active_slave_pool(where = "npqreg()")
    if (.npRmpi_npqreg_should_localize(bws) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
      return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    no.ex = missing(exdat)

    txdat = toFrame(txdat)
    tydat = toFrame(tydat)

    if (!is.numeric(tau) || length(tau) != 1 || is.na(tau) || tau <= 0 || tau >= 1)
      stop("'tau' must be a single numeric value in (0,1)")

    if (dim(tydat)[2] != 1)
      stop("'tydat' has more than one column")

    if (!no.ex){
      exdat = toFrame(exdat)
      
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
    }

    if(gradients)
      stop("gradients not currently supported for this object")
    
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

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.
    
    tydat = toMatrix(tydat)

    txdat = toMatrix(txdat)

    txuno = txdat[, bws$ixuno, drop = FALSE]
    txcon = txdat[, bws$ixcon, drop = FALSE]
    txord = txdat[, bws$ixord, drop = FALSE]

    if (!no.ex){
      exdat = toMatrix(exdat)

      exuno = exdat[, bws$ixuno, drop = FALSE]
      excon = exdat[, bws$ixcon, drop = FALSE]
      exord = exdat[, bws$ixord, drop = FALSE]
    } else {
      exuno = data.frame()
      excon = data.frame()
      exord = data.frame()
    }

    myopti = list(
      num_obs_train = tnrow,
      num_obs_eval = enrow,
      int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
      BANDWIDTH_den_extern = switch(bws$type,
        fixed = BW_FIXED,
        generalized_nn = BW_GEN_NN,
        adaptive_nn = BW_ADAP_NN),
      int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
      xkerneval = switch(bws$cxkertype,
        gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
        epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
        uniform = CKER_UNI,
        "truncated gaussian" = CKER_TGAUSS),
      ykerneval = switch(bws$cykertype,
        gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
        epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
        uniform = CKER_UNI,
        "truncated gaussian" = CKER_TGAUSS),
      uxkerneval = switch(bws$uxkertype,
        aitchisonaitken = UKER_AIT,
        liracine = UKER_LR),
      uykerneval = switch(bws$uykertype,
        aitchisonaitken = UKER_AIT,
        liracine = UKER_LR),
      oxkerneval = switch(bws$oxkertype,
        wangvanryzin = OKER_WANG,
        liracine = OKER_LR,
        "racineliyan" = OKER_RLY),
      oykerneval = switch(bws$oykertype,
        wangvanryzin = OKER_WANG,
        liracine = OKER_NLR,
        "racineliyan" = OKER_RLY),
      num_yuno = bws$ynuno,
      num_yord = bws$ynord,
      num_ycon = bws$yncon,
      num_xuno = bws$xnuno,
      num_xord = bws$xnord,
      num_xcon = bws$xncon,
      no.ex = no.ex,
      gradients = gradients,
      itmax = itmax,
      xmcv.numRow = attr(bws$xmcv, "num.row"),
      nmulti = itmax,
      qreg.unused = 0L)

    myoptd = c(
      qreg.unused = 0.0,
      tol = tol,
      small = small,
      rep(0.0, 7L))
    
    myout <-
      .Call("C_np_quantile_conditional",
            as.double(tydat),
            as.double(txuno), as.double(txord), as.double(txcon),
            as.double(exuno), as.double(exord), as.double(excon),
            as.double(tau),
            as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                        bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                        bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
            as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
            as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
            as.integer(myopti),
            as.double(myoptd),
            as.integer(enrow),
            as.integer(bws$xndim),
            as.logical(gradients),
            PACKAGE="npRmpi")[c("yq", "yqerr", "yqgrad")]

    if (all(!is.finite(myout$yqerr) | myout$yqerr <= 0.0)) {
      dens.bw <- tryCatch(
        conbandwidth(
          xbw = bws$xbw,
          ybw = bws$ybw,
          bwmethod = "manual",
          bwscaling = bws$scaling,
          bwtype = bws$type,
          cxkertype = bws$cxkertype,
          cxkerorder = bws$cxkerorder,
          cxkerbound = bws$cxkerbound,
          cxkerlb = bws$cxkerlb,
          cxkerub = bws$cxkerub,
          uxkertype = bws$uxkertype,
          oxkertype = bws$oxkertype,
          cykertype = bws$cykertype,
          cykerorder = bws$cykerorder,
          cykerbound = bws$cykerbound,
          cykerlb = bws$cykerlb,
          cykerub = bws$cykerub,
          uykertype = bws$uykertype,
          oykertype = bws$oykertype,
          nobs = nrow(txdat.df),
          xdati = bws$xdati,
          ydati = bws$ydati,
          xnames = bws$xnames,
          ynames = bws$ynames,
          sfactor = bws$sfactor,
          bandwidth = bws$bw,
          bandwidth.compute = FALSE
        ),
        error = function(e) NULL
      )
      dens.obj <- tryCatch(
        npcdens(
          txdat = txdat.df,
          tydat = tydat.df,
          exdat = if (no.ex) txdat.df else exdat.df,
          eydat = data.frame(y = as.double(myout$yq)),
          bws = dens.bw
        ),
        error = function(e) NULL
      )
      if (!is.null(dens.obj)) {
        dens.q <- as.double(dens.obj$condens)
        qvar <- tau * (1.0 - tau) / (tnrow * NZD(dens.q)^2)
        myout$yqerr <- sqrt(pmax(qvar, 0.0))
        myout$yqerr[!is.finite(myout$yqerr)] <- NA_real_
      }
    }

    ##need to untangle yqgrad

    if(gradients){
      myout$yqgrad = matrix(data=myout$yqgrad, nrow = enrow, ncol = bws$xndim, byrow = FALSE) 
      rorder = numeric(bws$xndim)
      xidx <- seq_len(bws$xndim)
      rorder[c(xidx[bws$ixcon], xidx[bws$ixuno], xidx[bws$ixord])] <- xidx
      myout$yqgrad = myout$yqgrad[, rorder, drop = FALSE]

    } else {
      myout$yqgrad = NA
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
                ntrain = tnrow,
                trainiseval = no.ex,
                gradients = gradients,
                timing = bws$timing, total.time = total.time,
                optim.time = optim.time, fit.time = fit.elapsed)
  }


npqreg.default <- function(bws, txdat, tydat, ...){
  .npRmpi_require_active_slave_pool(where = "npqreg()")
  if (!missing(bws) &&
      .npRmpi_npqreg_should_localize(bws) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
  if (.npRmpi_autodispatch_active())
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

  if (!missing(bws) && inherits(bws, "formula")) {
    dots <- list(...)
    tbw <- do.call(npcdistbw, c(list(formula = bws), dots))
    return(npqreg(bws = tbw, ...))
  }

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

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npcdistbw)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat')
  nstxy <- c('xdat','ydat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }
    
  tbw <- if (!has.explicit.bws) {
    .np_progress_select_bandwidth(
      "Selecting conditional distribution bandwidth",
      .np_eval_bw_call(sc.bw, caller_env = parent.frame())
    )
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
