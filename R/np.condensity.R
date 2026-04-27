npcdens <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npcdens",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npcdens",bws$call)
        else if (!is.call(bws))
          UseMethod("npcdens",bws)
        else
          UseMethod("npcdens",NULL)
      } else {
        UseMethod("npcdens", NULL)
      }
    } else {
      UseMethod("npcdens", NULL)
    }
  }

npcdens.formula <-
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
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      emf <- umf

      eydat <- emf[, bws$variableNames[["response"]], drop = FALSE]
      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }

    cd.args <- list(txdat = txdat, tydat = tydat)
    if (has.eval) {
      cd.args$exdat <- exdat
      cd.args$eydat <- eydat
    }
    cd.args$bws <- bws
    ev <- do.call(npcdens, c(cd.args, list(...)))

    ev$omit <- attr(umf,"na.action")
    ev$rows.omit <- as.vector(ev$omit)
    ev$nobs.omit <- length(ev$rows.omit)

    ev$condens <- napredict(ev$omit, ev$condens)
    ev$conderr <- napredict(ev$omit, ev$conderr)
    if (!is.null(ev$condens.raw))
      ev$condens.raw <- napredict(ev$omit, ev$condens.raw)

    if(ev$gradients){
        ev$congrad <- napredict(ev$omit, ev$congrad)
        ev$congerr <- napredict(ev$omit, ev$congerr)
    }

    return(ev)
  }

npcdens.call <-
  function(bws, ...) {
    npcdens(txdat = .np_eval_bws_call_arg(bws, "xdat"),
            tydat = .np_eval_bws_call_arg(bws, "ydat"),
            bws = bws, ...)
  }


npcdens.conbandwidth <- function(bws,
                                 txdat = stop("invoked without training data 'txdat'"),
                                 tydat = stop("invoked without training data 'tydat'"),
                                 exdat, eydat, gradients = FALSE,
                                 proper = FALSE,
                                 proper.method = c("project"),
                                 proper.control = list(),
                                 ...){

  dots <- list(...)
  fit.start <- proc.time()[3]
  fit.progress.handoff <- isTRUE(dots$.np_fit_progress_handoff)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  proper.args <- .np_condens_validate_proper_args(
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control
  )
  .npRmpi_require_active_slave_pool(where = "npcdens()")
  .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = "npcdens()")
  keep_local_shadow_nn <- identical(bws$regtype.engine, "lp") &&
    identical(bws$type %in% c("generalized_nn", "adaptive_nn"), TRUE)
  if (.npRmpi_autodispatch_active() &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
      identical(.npRmpi_safe_int(mpi.comm.size(0)), 1L)) {
    return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
  }
  if (.npRmpi_autodispatch_active() && !keep_local_shadow_nn) {
    out <- .npRmpi_autodispatch_call(match.call(), parent.frame())
    out <- .npRmpi_restore_nomad_fit_bws_metadata(out, bws)
    if (inherits(out, "condensity") &&
        !is.null(out$proper.requested) &&
        !is.null(out$proper.applied) &&
        !is.null(out$proper.info))
      return(out)
    return(.np_condens_finalize_proper_object(
      object = out,
      proper = proper.args$proper.requested,
      proper.method = proper.args$proper.method,
      proper.control = proper.args$proper.control,
      where = "npcdens()"
    ))
  }

  if (xor(missing(exdat),missing(eydat)))
    stop("evaluation data must be supplied for both 'exdat' and 'eydat'")

  no.exy = missing(exdat)

  txdat = toFrame(txdat)
  tydat = toFrame(tydat)

  if (!no.exy){
    exdat = toFrame(exdat)
    eydat = toFrame(eydat)

    if (! txdat %~% exdat )
      stop("'txdat' and 'exdat' are not similar data frames!")

    if (! tydat %~% eydat )
      stop("'tydat' and 'eydat' are not similar data frames!")

  }

  if (length(bws$xbw) != length(txdat))
    stop("length of bandwidth vector does not match number of columns of 'txdat'")

  if (length(bws$ybw) != length(tydat))
    stop("length of bandwidth vector does not match number of columns of 'tydat'")

  if ((any(bws$ixcon) &&
       !all(vapply(txdat[, bws$ixcon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
      (any(bws$ixord) &&
       !all(vapply(txdat[, bws$ixord, drop = FALSE], inherits, logical(1), "ordered"))) ||
      (any(bws$ixuno) &&
       !all(vapply(txdat[, bws$ixuno, drop = FALSE], inherits, logical(1), "factor"))))
    stop("supplied bandwidths do not match 'txdat' in type")

  if ((any(bws$iycon) &&
       !all(vapply(tydat[, bws$iycon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
      (any(bws$iyord) &&
       !all(vapply(tydat[, bws$iyord, drop = FALSE], inherits, logical(1), "ordered"))) ||
      (any(bws$iyuno) &&
       !all(vapply(tydat[, bws$iyuno, drop = FALSE], inherits, logical(1), "factor"))))
    stop("supplied bandwidths do not match 'tydat' in type")
  
  ## catch and destroy NA's
  keep.rows <- rep_len(TRUE, nrow(txdat))
  rows.omit <- attr(na.omit(data.frame(txdat, tydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE

  if (!any(keep.rows))
    stop("Data has no rows without NAs")

  txdat <- txdat[keep.rows,,drop = FALSE]
  tydat <- tydat[keep.rows,,drop = FALSE]

  if (!no.exy){
    keep.eval <- rep_len(TRUE, nrow(exdat))
    rows.omit <- attr(na.omit(data.frame(exdat, eydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.eval[as.integer(rows.omit)] <- FALSE

    if (!any(keep.eval))
      stop("Data has no rows without NAs")

    exdat <- exdat[keep.eval,,drop = FALSE]
    eydat <- eydat[keep.eval,,drop = FALSE]
  }


  tnrow = nrow(txdat)
  enrow = (if (no.exy) tnrow else nrow(exdat))

  ## re-assign levels in training and evaluation data to ensure correct
  ## conversion to numeric type.
  
  txdat <- adjustLevels(txdat, bws$xdati)
  tydat <- adjustLevels(tydat, bws$ydati)
  
  if (!no.exy){
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    eydat <- adjustLevels(eydat, bws$ydati, allowNewCells = TRUE)
    npKernelBoundsCheckEval(exdat, bws$ixcon, bws$cxkerlb, bws$cxkerub, argprefix = "cxker")
    npKernelBoundsCheckEval(eydat, bws$iycon, bws$cykerlb, bws$cykerub, argprefix = "cyker")
  }

  proper.slice.context <- list(
    txdat = txdat,
    tydat = tydat,
    exdat = if (no.exy) NULL else exdat,
    eydat = if (no.exy) NULL else eydat
  )

  ## grab the evaluation data before it is converted to numeric
  if(no.exy){
    txeval <- txdat
    tyeval <- tydat
  } else {
    txeval <- exdat
    tyeval <- eydat
  }

  ## at this stage, data to be sent to the c routines must be converted to
  ## numeric type.
  
  tydat = toMatrix(tydat)

  tyuno = tydat[, bws$iyuno, drop = FALSE]
  tycon = tydat[, bws$iycon, drop = FALSE]
  tyord = tydat[, bws$iyord, drop = FALSE]


  txdat = toMatrix(txdat)

  txuno = txdat[, bws$ixuno, drop = FALSE]
  txcon = txdat[, bws$ixcon, drop = FALSE]
  txord = txdat[, bws$ixord, drop = FALSE]

  if (!no.exy){
    eydat = toMatrix(eydat)

    eyuno = eydat[, bws$iyuno, drop = FALSE]
    eycon = eydat[, bws$iycon, drop = FALSE]
    eyord = eydat[, bws$iyord, drop = FALSE]


    exdat = toMatrix(exdat)

    exuno = exdat[, bws$ixuno, drop = FALSE]
    excon = exdat[, bws$ixcon, drop = FALSE]
    exord = exdat[, bws$ixord, drop = FALSE]
  } else {
    eyuno = data.frame()
    eycon = data.frame()
    eyord = data.frame()

    exuno = data.frame()
    excon = data.frame()
    exord = data.frame()
  }

  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  basis.engine <- if (is.null(bws$basis.engine)) {
    if (is.null(bws$basis)) "glp" else bws$basis
  } else {
    bws$basis.engine
  }
  degree.engine <- if (is.null(bws$degree.engine)) {
    if (bws$xncon > 0L) {
      if (identical(reg.engine, "lc")) rep.int(0L, bws$xncon) else npValidateGlpDegree(
        regtype = "lp",
        degree = bws$degree,
        ncon = bws$xncon
      )
    } else {
      integer(0)
    }
  } else {
    as.integer(bws$degree.engine)
  }
  bernstein.engine <- if (is.null(bws$bernstein.basis.engine)) {
    isTRUE(bws$bernstein.basis)
  } else {
    isTRUE(bws$bernstein.basis.engine)
  }

  reg.c <- npRegtypeToC(
    regtype = if (identical(reg.engine, "lp")) "lp" else "lc",
    degree = degree.engine,
    ncon = bws$xncon,
    context = "npcdens"
  )
  degree.c <- if (bws$xncon > 0L) {
    as.integer(if (is.null(reg.c$degree)) rep.int(0L, bws$xncon) else reg.c$degree)
  } else {
    integer(0)
  }
  basis.code <- as.integer(npLpBasisCode(basis.engine))

  myopti <- list(
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
          liracine = OKER_NLR,
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
      no.exy = no.exy,
      gradients = gradients,
      ymcv.numRow = attr(bws$ymcv, "num.row"),
      xmcv.numRow = attr(bws$xmcv, "num.row"),
      densOrDist = NP_DO_DENS,
      int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO)

  cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
  cyker.bounds.c <- npKernelBoundsMarshal(bws$cykerlb[bws$iycon], bws$cykerub[bws$iycon])

  myout <- .np_with_compiled_fit_progress(
    label = "Fitting conditional density",
    total = .np_condensdist_fit_total(bws = bws, tnrow = tnrow, enrow = enrow),
    handoff = fit.progress.handoff,
    handoff.detail = if (fit.progress.handoff) "starting" else NULL,
    if (keep_local_shadow_nn) {
      .npRmpi_with_local_regression(
        .Call("C_np_density_conditional",
              as.double(tyuno), as.double(tyord), as.double(tycon),
              as.double(txuno), as.double(txord), as.double(txcon),
              as.double(eyuno), as.double(eyord), as.double(eycon),
              as.double(exuno), as.double(exord), as.double(excon),
              as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                          bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                          bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
              as.double(bws$ymcv), as.double(attr(bws$ymcv, "pad.num")),
              as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
              as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
              as.integer(myopti),
              as.integer(enrow),
              as.integer(bws$xndim),
              as.double(cxker.bounds.c$lb),
              as.double(cxker.bounds.c$ub),
              as.double(cyker.bounds.c$lb),
              as.double(cyker.bounds.c$ub),
              as.integer(reg.c$code),
              as.integer(degree.c),
              as.integer(bernstein.engine),
              basis.code,
              PACKAGE = "npRmpi")
      )
    } else {
      .Call("C_np_density_conditional",
            as.double(tyuno), as.double(tyord), as.double(tycon),
            as.double(txuno), as.double(txord), as.double(txcon),
            as.double(eyuno), as.double(eyord), as.double(eycon),
            as.double(exuno), as.double(exord), as.double(excon),
            as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                        bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                        bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
            as.double(bws$ymcv), as.double(attr(bws$ymcv, "pad.num")),
            as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
            as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
            as.integer(myopti),
            as.integer(enrow),
            as.integer(bws$xndim),
            as.double(cxker.bounds.c$lb),
            as.double(cxker.bounds.c$ub),
            as.double(cyker.bounds.c$lb),
            as.double(cyker.bounds.c$ub),
            as.integer(reg.c$code),
            as.integer(degree.c),
            as.integer(bernstein.engine),
            basis.code,
            PACKAGE = "npRmpi")
    }
  )

  if(gradients){
    myout$congrad = matrix(data=myout$congrad, nrow = enrow, ncol = bws$xndim, byrow = FALSE) 
    rorder = numeric(bws$xndim)
    xidx <- seq_len(bws$xndim)
    rorder[c(xidx[bws$ixcon], xidx[bws$ixuno], xidx[bws$ixord])] <- xidx
    myout$congrad = myout$congrad[, rorder, drop = FALSE]

    myout$congerr = matrix(data=myout$congerr, nrow = enrow, ncol = bws$xndim, byrow = FALSE)
    myout$congerr = myout$congerr[, rorder, drop = FALSE]
  } else {
    myout$congrad = NA
    myout$congerr = NA
  }


  fit.elapsed <- proc.time()[3] - fit.start
  optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
  total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)

  out <- condensity(bws = bws,
                    xeval = txeval,
                    yeval = tyeval,
                    condens = myout$condens, conderr = myout$conderr,
                    congrad = myout$congrad, congerr = myout$congerr,
                    ll = myout$log_likelihood,
                    ntrain = tnrow, trainiseval = no.exy, gradients = gradients,
                    rows.omit = rows.omit,
                    timing = bws$timing, total.time = total.time,
                    optim.time = optim.time, fit.time = fit.elapsed)

  .np_condens_finalize_proper_object(
    object = out,
    proper = proper.args$proper.requested,
    proper.method = proper.args$proper.method,
    proper.control = proper.args$proper.control,
    slice.context = proper.slice.context,
    where = "npcdens()"
  )

}


npcdens.default <- function(bws, txdat, tydat, nomad = FALSE, ...){
  .npRmpi_require_active_slave_pool(where = "npcdens()")
  .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = "npcdens()")
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
  has.explicit.bws <- (!no.bws) && isa(bws, "conbandwidth")
  bws.formula <- (!no.bws) && inherits(bws, "formula")
  regtype.request <- if (has.explicit.bws) {
    if (is.null(bws$regtype.engine)) {
      if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
    } else {
      as.character(bws$regtype.engine)
    }
  } else if ("regtype" %in% sc.names) {
    as.character(eval(sc$regtype, parent.frame()))
  } else {
    "lc"
  }
  bwtype.request <- if (has.explicit.bws) {
    if (is.null(bws$type)) "fixed" else as.character(bws$type)
  } else if ("bwtype" %in% sc.names) {
    as.character(eval(sc$bwtype, parent.frame()))
  } else {
    "fixed"
  }
  bwmethod.request <- if (has.explicit.bws) {
    if (is.null(bws$method)) "cv.ls" else as.character(bws$method)
  } else if ("bwmethod" %in% sc.names) {
    as.character(eval(sc$bwmethod, parent.frame()))
  } else {
    "cv.ls"
  }
  degree.request <- if ("degree" %in% sc.names) eval(sc$degree, parent.frame()) else NULL
  bernstein.request <- if ("bernstein.basis" %in% sc.names) {
    isTRUE(eval(sc$bernstein.basis, parent.frame()))
  } else {
    FALSE
  }
  if (has.explicit.bws &&
      .npRmpi_autodispatch_active() &&
      !isTRUE(nomad) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
      identical(.npRmpi_safe_int(mpi.comm.size(0)), 1L)) {
    return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
  }
  keep_local_shadow_nn <- (identical(regtype.request[1L], "lp") ||
    identical(regtype.request[1L], "ll")) &&
    identical(bwtype.request[1L] %in% c("generalized_nn", "adaptive_nn"), TRUE)
  keep_local_raw_degree1_cvls <- !has.explicit.bws &&
    identical(bwmethod.request[1L], "cv.ls") &&
    identical(bwtype.request[1L], "fixed") &&
    npIsRawDegreeOneConditionalRequest(
      regtype = regtype.request[1L],
      degree = degree.request,
      bernstein.basis = bernstein.request
    )
  if (.npRmpi_autodispatch_active() &&
      !isTRUE(nomad) &&
      !keep_local_shadow_nn &&
      !keep_local_raw_degree1_cvls)
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

  ## autodispatch normalizes calls via match.call(), which can turn an
  ## originally unnamed formula first argument into named bws=... .
  ## Preserve legacy formula behavior by rewriting npcdensbw() call shape.
  if (bws.named && no.txdat && no.tydat && inherits(bws, "formula")) {
    sc$`bws` <- NULL
    sc$formula <- bws
    sc.bw <- sc
    sc.bw[[1]] <- quote(npcdensbw)
    bws.named <- FALSE
  } else {
    sc.bw <- sc
    sc.bw[[1]] <- quote(npcdensbw)
  }

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  if(bws.named && !bws.formula){
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
        "Selecting conditional density bandwidth",
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
  if (!has.explicit.bws)
    call.args$.np_fit_progress_handoff <- TRUE
  do.call(npcdens, c(call.args, list(...)))
}
