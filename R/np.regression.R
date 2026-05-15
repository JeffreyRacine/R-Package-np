npreg <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (inherits(bws, "formula") && is.null(args$txdat))
        UseMethod("npreg", bws)
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npreg",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npreg",bws$call)
        else if (!is.call(bws))
          UseMethod("npreg",bws)
        else
          UseMethod("npreg",NULL)
      } else {
        UseMethod("npreg", NULL)
      }
    } else {
      UseMethod("npreg", NULL)
    }
  }

npreg.formula <-
  function(bws, data = NULL, newdata = NULL, y.eval = FALSE, ...){

    tt <- terms(bws)
    tmf <- if (!is.null(bws$call)) {
      m <- match(c("formula", "data", "subset", "na.action"),
                 names(bws$call), nomatch = 0)
      bws$call[c(1, m)]
    } else {
      mc <- match.call(expand.dots = FALSE)
      m <- match(c("bws", "data", "subset", "na.action"),
                 names(mc), nomatch = 0)
      tmf <- mc[c(1, m)]
      if ("bws" %in% names(tmf))
        names(tmf)[names(tmf) == "bws"] <- "formula"
      tmf
    }
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    response.name <- attr(tmf, "names")[attr(attr(tmf, "terms"), "response")]
    tydat <- model.response(tmf)
    txdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
    has.eval <- !is.null(newdata)
    if (has.eval) {
      if (!y.eval){
        npValidateNewdataFormula(newdata, tt, include.response = FALSE)
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
      
      if (y.eval)
        npValidateNewdataFormula(newdata, tt, include.response = TRUE)
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      emf <- umf

      if (y.eval)
        eydat <- model.response(emf)
      
      exdat <- emf[, attr(attr(emf, "terms"),"term.labels"), drop = FALSE]
    }

    reg.args <- list(txdat = txdat, tydat = tydat, bws = bws)
    if (has.eval) {
      reg.args$exdat <- exdat
      if (y.eval)
        reg.args$eydat <- eydat
    }
    ev <- do.call(npreg, c(reg.args, list(...)))
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()

    if (length(response.name) == 1L && !is.na(response.name) && nzchar(response.name)) {
      if (!is.null(ev$bws))
        ev$bws$ynames <- response.name
    }

    ev$omit <- attr(umf,"na.action")
    ev$rows.omit <- as.vector(ev$omit)
    ev$nobs.omit <- length(ev$rows.omit)

    ev$mean <- napredict(ev$omit, ev$mean)
    ev$merr <- napredict(ev$omit, ev$merr)

    if(ev$gradients){
        ev$grad <- napredict(ev$omit, ev$grad)
        ev$gerr <- napredict(ev$omit, ev$gerr)
    }

    if(ev$residuals){
        ev$resid <- naresid(ev$omit, ev$resid)
    }    
    return(ev)
  }

npreg.call <-
  function(bws, ...) {
    ev <- npreg(txdat = .np_eval_bws_call_arg(bws, "xdat"),
                tydat = .np_eval_bws_call_arg(bws, "ydat"),
                bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

.np_reg_fit_total <- function(bws, tnrow, enrow) {
  if (identical(as.character(bws$type)[1L], "adaptive_nn")) {
    as.integer(tnrow)
  } else {
    as.integer(enrow)
  }
}

npreg.rbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat, eydat, gradient.order = 1L, gradients = FALSE,
           residuals = FALSE,
           ...){
    fit.start <- proc.time()[3]

    no.ex = missing(exdat)
    no.ey = missing(eydat)
    dots <- list(...)
    fit.progress.handoff <- isTRUE(dots$.np_fit_progress_handoff)
    if ("remin" %in% names(dots)) {
      warning("npreg: bandwidth-selection argument 'remin' is ignored when a bandwidth object is supplied",
              call. = FALSE)
      dots$remin <- NULL
    }
    npRejectLegacyLpArgs(names(dots), where = "npreg")
    warn.glp.gradient <- if (is.null(dots$warn.glp.gradient)) TRUE else isTRUE(dots$warn.glp.gradient)

    gradients <- npValidateScalarLogical(gradients, "gradients")
    residuals <- npValidateScalarLogical(residuals, "residuals")

    txdat = toFrame(txdat)

    if (!(is.vector(tydat) || is.factor(tydat)))
      stop("'tydat' must be a vector or a factor")


    ## if no.ex then if !no.ey then ey and tx must match, to get oos errors
    ## alternatively if no.ey you get is errors
    ## if !no.ex then if !no.ey then ey and ex must match, to get oos errors
    ## alternatively if no.ey you get NO errors since we don't evaluate on the training
    ## data
    
     
    if (!no.ex){
      exdat = toFrame(exdat)

      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")

      if (!no.ey){
        if (!(is.vector(eydat) || is.factor(eydat)))
          stop("'eydat' must be a vector or a factor")
        if (dim(exdat)[1] != length(eydat))
          stop("number of evaluation data 'exdat' and dependent data 'eydat' do not match")
        if (!identical(coarseclass(eydat),coarseclass(tydat)))
          stop("type of evaluation data 'eydat' does not match that of 'tydat'")
      }
      
    } else if(!no.ey) {
      if (dim(txdat)[1] != length(eydat))
        stop("number of training data 'txdat' and dependent data 'eydat' do not match")
    }

    if (length(bws$bw) != length(txdat))
      stop("length of bandwidth vector does not match number of columns of 'txdat'")

    npValidateRegressionNnLowerBound(bws, where = "npreg")

    bws$basis <- npValidateLpBasis(regtype = bws$regtype,
                                   basis = bws$basis)
    bws$degree <- npValidateGlpDegree(regtype = bws$regtype,
                                      degree = bws$degree,
                                      ncon = bws$ncon)
    bws$bernstein.basis <- npValidateGlpBernstein(regtype = bws$regtype,
                                                  bernstein.basis = bws$bernstein.basis)
    reg.spec <- npCanonicalConditionalRegSpec(
      regtype = bws$regtype,
      basis = bws$basis,
      degree = bws$degree,
      bernstein.basis = bws$bernstein.basis,
      ncon = bws$ncon,
      where = "npreg"
    )
    glp.gradient.order <- if (identical(reg.spec$regtype.engine, "lp")) {
      if (identical(bws$regtype, "lp")) {
        npValidateGlpGradientOrder(regtype = bws$regtype,
                                   gradient.order = gradient.order,
                                   ncon = bws$ncon)
      } else if (bws$ncon > 0L) {
        rep.int(1L, bws$ncon)
      } else {
        integer(0)
      }
    } else {
      NULL
    }
    if (isTRUE(gradients) &&
        identical(reg.spec$regtype.engine, "lp") &&
        (bws$ncon > 0L) &&
        all(reg.spec$degree.engine == 0L)) {
      stop("regtype='lp' with degree=0 does not support derivatives; use gradients=FALSE for fitted/predicted values")
    }

    reg.c <- npRegtypeToC(regtype = reg.spec$regtype.engine,
                          degree = reg.spec$degree.engine,
                          ncon = bws$ncon,
                          context = "npreg")
    degree.c <- if (bws$ncon > 0) {
      as.integer(if (is.null(reg.c$degree)) rep.int(0L, bws$ncon) else reg.c$degree)
    } else {
      integer(1)
    }

    if ((any(bws$icon) &&
         !all(vapply(txdat[, bws$icon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iord) &&
         !all(vapply(txdat[, bws$iord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$iuno) &&
         !all(vapply(txdat[, bws$iuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'txdat' in type")

    if (dim(txdat)[1] != length(tydat))
      stop("number of explanatory data 'txdat' and dependent data 'tydat' do not match")

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(data.frame(txdat,tydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows]

    ## no.ex = missing(exdat)
    ## no.ey = missing(eydat)

    if (!no.ex){
      keep.eval <- rep_len(TRUE, nrow(exdat))
      eval.df <- data.frame(exdat)
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
    }

    if (identical(bws$regtype, "lp") &&
        isTRUE(bws$bernstein.basis) &&
        !no.ex &&
        any(bws$icon)) {
      out.of.support <- character(0)
      for (ii in which(bws$icon)) {
        tr <- range(as.numeric(txdat[[ii]]))
        ex <- as.numeric(exdat[[ii]])
        if (any(ex < tr[1] | ex > tr[2], na.rm = TRUE))
          out.of.support <- c(out.of.support, colnames(txdat)[ii])
      }
      if (length(out.of.support)) {
        .np_warning(
          "bernstein.basis=TRUE: evaluation continuous predictor(s) outside training support (",
          paste(unique(out.of.support), collapse = ", "),
          "); proceeding, but interpret extrapolated results with care",
          call. = FALSE
        )
      }
    }

    ## evaluate residuals before data conversion ...

    if (residuals){
      resid <- tydat - npreg(txdat = txdat, tydat = tydat, bws = bws)$mean
    }


    tnrow = dim(txdat)[1]
    enrow = (if (no.ex) tnrow else dim(exdat)[1])
    ncol = dim(txdat)[2]

    ## convert tydat, eydat to numeric, from a factor with levels from the y-data
    ## used during bandwidth selection.
    
    if (is.factor(tydat)){
      tydat <- adjustLevels(data.frame(tydat), bws$ydati)[,1]
      tydat <- (bws$ydati$all.dlev[[1]])[as.integer(tydat)]
    }
    else
      tydat <- as.double(tydat)


    if (no.ey)
      eydat <- double()
    else {
      if (is.factor(eydat)){
        eydat <- adjustLevels(data.frame(eydat), bws$ydati, allowNewCells = TRUE)
        eydat <- toMatrix(eydat)[,1]
      }
      else
        eydat <- as.double(eydat)
    }

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati)
      
    if (!no.ex)
      exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)

    if (!no.ex)
      npKernelBoundsCheckEval(exdat, bws$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")

    ## grab the evaluation data before it is converted to numeric
    if(no.ex)
      teval <- txdat
    else
      teval <- exdat

    ## put the unordered, ordered, and continuous data in their own objects
    ## data that is not a factor is continuous.

    txdat = toMatrix(txdat)

    tuno = txdat[, bws$iuno, drop = FALSE]
    tcon = txdat[, bws$icon, drop = FALSE]
    tord = txdat[, bws$iord, drop = FALSE]

    npCheckRegressionDesignCondition(reg.code = reg.c$code,
                                     xcon = tcon,
                                     basis = reg.spec$basis.engine,
                                     degree = reg.spec$degree.engine,
                                     bernstein.basis = reg.spec$bernstein.basis.engine,
                                     where = "npreg")

    if (!no.ex){
      exdat = toMatrix(exdat)

      euno = exdat[, bws$iuno, drop = FALSE]
      econ = exdat[, bws$icon, drop = FALSE]
      eord = exdat[, bws$iord, drop = FALSE]

    } else {
      euno = data.frame()
      eord = data.frame()
      econ = data.frame()
    }

    myopti = list(
      num_obs_train = tnrow,
      num_obs_eval = enrow,
      num_uno = bws$nuno, num_ord = bws$nord,
      num_con = bws$ncon,
      int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
      BANDWIDTH_reg_extern = switch(bws$type,
        fixed = BW_FIXED,
        generalized_nn = BW_GEN_NN,
        adaptive_nn = BW_ADAP_NN),
      int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE, 
      kerneval = switch(bws$ckertype,
        gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
        epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
        uniform = CKER_UNI,
        "truncated gaussian" = CKER_TGAUSS),
      ukerneval = switch(bws$ukertype,
        aitchisonaitken = UKER_AIT,
        liracine = UKER_LR),
      okerneval = switch(bws$okertype,
        wangvanryzin = OKER_WANG,
        liracine = OKER_LR,
        "racineliyan" = OKER_RLY),
      ey_is_ty = no.ey,
      do_grad = gradients,
      regtype = reg.c$code,
      no.ex = no.ex,
      mcv.numRow = attr(bws$xmcv, "num.row"),
      int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
      old.reg = FALSE)

    cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])
    
   asDouble <- function(data){
	   if (is.null(data)){
	 	   result <- as.double(0.0)
	   }
	   else {
		   result <- as.double(data)
	   }
	   return(result)
   }


    glp.gradient.order.c <- if (bws$ncon > 0L) {
      as.integer(if (is.null(glp.gradient.order)) rep.int(1L, bws$ncon) else glp.gradient.order)
    } else {
      integer(1)
    }

    myout <- .np_with_compiled_fit_progress(
      label = "Fitting regression",
      total = .np_reg_fit_total(bws = bws, tnrow = tnrow, enrow = enrow),
      handoff = fit.progress.handoff,
      handoff.detail = if (fit.progress.handoff) "starting" else NULL,
      .Call("C_np_regression",
            asDouble(tuno), asDouble(tord), asDouble(tcon), asDouble(tydat),
            asDouble(euno), asDouble(eord), asDouble(econ), asDouble(eydat),
            asDouble(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
            asDouble(bws$xmcv), asDouble(attr(bws$xmcv, "pad.num")),
            asDouble(bws$nconfac), asDouble(bws$ncatfac), asDouble(bws$sdev),
            as.integer(myopti),
            as.integer(degree.c),
            as.integer(glp.gradient.order.c),
            as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
            as.integer(npLpBasisCode(reg.spec$basis.engine)),
            as.integer(enrow),
            as.integer(ncol),
            as.logical(gradients),
            as.double(cker.bounds.c$lb),
            as.double(cker.bounds.c$ub),
            PACKAGE = "np")
    )

    if (gradients){
      myout$g = matrix(data=myout$g, nrow = enrow, ncol = ncol, byrow = FALSE) 
      rorder = numeric(ncol)
      ord_idx <- seq_len(ncol)
      rorder[c(ord_idx[bws$icon], ord_idx[bws$iuno], ord_idx[bws$iord])] <- ord_idx
      myout$g = as.matrix(myout$g[,rorder])

      myout$gerr = matrix(data=myout$gerr, nrow = enrow, ncol = ncol, byrow = FALSE) 
      myout$gerr = as.matrix(myout$gerr[,rorder])

      if (identical(bws$regtype, "lp")) {
        cont.idx <- which(bws$icon)
        if (length(cont.idx)) {
          invalid.order <- glp.gradient.order > bws$degree
          if (any(invalid.order)) {
            bad.idx <- cont.idx[invalid.order]
            myout$g[, bad.idx] <- NA_real_
            myout$gerr[, bad.idx] <- NA_real_
            if (warn.glp.gradient)
              .np_warning("some requested glp derivatives exceed polynomial degree; returning NA for those components")
          }
        }
      }
    }


    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)

    ev.args <- list(
      bws = bws,
      eval = teval,
      mean = myout$mean,
      merr = myout$merr,
      ntrain = tnrow,
      trainiseval = no.ex,
      gradients = gradients,
      residuals = residuals,
      xtra = myout$xtra,
      rows.omit = rows.omit,
      timing = bws$timing,
      total.time = total.time,
      optim.time = optim.time,
      fit.time = fit.elapsed
    )
    if (gradients) {
      ev.args$grad <- myout$g
      ev.args$gerr <- myout$gerr
    }
    if (residuals)
      ev.args$resid <- resid
    if (identical(bws$regtype, "lp"))
      ev.args$gradient.order <- glp.gradient.order
    ev <- do.call(npregression, ev.args)
    ev$nomad.time <- if (!is.null(bws$nomad.time) && is.finite(bws$nomad.time)) as.double(bws$nomad.time) else NA_real_
    ev$powell.time <- if (!is.null(bws$powell.time) && is.finite(bws$powell.time)) as.double(bws$powell.time) else NA_real_


    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npreg.default <- function(bws, txdat, tydat, nomad = FALSE, ...){
  sc <- sys.call()
  sc.names <- names(sc)
  nomad <- npValidateScalarLogical(nomad, "nomad")

  if (!missing(bws) &&
      !isa(bws, "rbandwidth") &&
      (inherits(bws, "formula") || is.call(bws))) {
    dots <- list(...)
    dots$nomad <- nomad
    bw.args <- if (missing(txdat) && missing(tydat)) {
      list(formula = bws)
    } else {
      list(xdat = txdat, ydat = tydat)
    }
    tbw <- do.call(npregbw, c(bw.args, dots))
    reg.args <- list(bws = tbw)
    if (!missing(txdat))
      reg.args$txdat <- txdat
    if (!missing(tydat))
      reg.args$tydat <- tydat
    fit.dots <- dots
    fit.dots$remin <- NULL
    fit.dots$.np_fit_progress_handoff <- TRUE
    return(do.call(npreg, c(reg.args, fit.dots)))
  }

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  has.explicit.bws <- (!no.bws) && isa(bws, "rbandwidth")

  ## match.call() can normalize an originally unnamed formula first
  ## argument into named bws=... . Rewrite the bw call so formula-driven
  ## fits continue down the formula route instead of disabling search.
  if (bws.named && no.txdat && no.tydat && inherits(bws, "formula")) {
    sc$`bws` <- NULL
    sc$formula <- bws
    sc.bw <- sc
    sc.bw[[1]] <- quote(npregbw)
    bws.named <- FALSE
  } else {
    sc.bw <- sc
    sc.bw[[1]] <- quote(npregbw)
  }

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

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
        "Selecting regression bandwidth",
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
  ev <- do.call(npreg, c(call.args, list(...)))

  ev$call <- match.call(expand.dots = FALSE)
  environment(ev$call) <- parent.frame()
  return(ev)
}
