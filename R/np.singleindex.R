# We implement Ichimura's single index model and Klein and Spady's
# single index model using npksum() and the nlm() minimization
# routine in R. These semiparametric models are used to reduce
# dimensionality to a one-dimensional nonparametric estimator, though
# at the potential cost of misspecification.

# Note also that we will use the so-called scale normalization, i.e.,
# that beta_1=1 (no need to estimate) which reduces search by 1
# parameter (this is obviously restricted search subject to beta_1=1).

# Define the index function model... it is a simple local constant
# estimator of y on a linear index X\beta where beta_1 is presumed to
# be 1 by restriction though, at this stage, the user may feed in any
# value they so desire.

npindex <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (length(args) > 0L &&
          inherits(args[[1L]], "formula") &&
          is.null(args$txdat)) {
        formula <- args[[1L]]
        args <- args[-1L]
        args$.np_index_explicit_bws <- bws
        return(do.call(npindex.formula,
                       c(list(bws = formula), args),
                       envir = parent.frame()))
      }
      if (inherits(bws, "formula") && is.null(args$txdat))
        UseMethod("npindex", bws)
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npindex",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npindex",bws$call)
        else if (!is.call(bws))
          UseMethod("npindex",bws)
        else
          UseMethod("npindex",NULL)
      } else {
        UseMethod("npindex", NULL)
      }
    } else {
      UseMethod("npindex", NULL)
    }
  }

npindex.formula <-
    function(bws, data = NULL, newdata = NULL, y.eval = FALSE, ...){

        mc <- match.call(expand.dots = FALSE)
        tt <- terms(bws)
        tmf <- if (!is.null(bws$call)) {
          m <- match(c("formula", "data", "subset", "na.action"),
                     names(bws$call), nomatch = 0)
          bws$call[c(1, m)]
        } else {
          m <- match(c("bws", "data", "subset", "na.action"),
                     names(mc), nomatch = 0)
          tmf <- mc[c(1, m)]
          if ("bws" %in% names(tmf))
            names(tmf)[names(tmf) == "bws"] <- "formula"
          tmf
        }
        tmf[[1]] <- as.name("model.frame")
        tmf[["formula"]] <- tt
        if (!missing(data) && !is.null(data))
          tmf[["data"]] <- substitute(data)
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

        dots <- list(...)
        si.bws <- if (!is.null(dots$.np_index_explicit_bws)) {
            out <- dots$.np_index_explicit_bws
            dots$.np_index_explicit_bws <- NULL
            out
        } else if (!is.null(dots$bws)) {
            out <- dots$bws
            dots$bws <- NULL
            out
        } else {
            bws
        }

        si.args <- list(txdat = txdat, tydat = tydat)
        if (has.eval) {
          si.args$exdat <- exdat
          if (y.eval)
            si.args$eydat <- eydat
        }
        si.args$bws <- si.bws
        ev <- do.call(npindex, c(si.args, dots))
        ev$call <- mc
        environment(ev$call) <- parent.frame()

        if (length(response.name) == 1L && !is.na(response.name) && nzchar(response.name)) {
            if (!is.null(ev$bws))
                ev$bws$ynames <- response.name
        }
        if (!is.null(ev$bws) && inherits(bws, "formula")) {
            ev$bws$formula <- bws
            ev$bws$terms <- attr(tmf, "terms")
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

npindex.call <-
  function(bws, ...) {
    npindex(txdat = .np_eval_bws_call_arg(bws, "xdat"),
            tydat = .np_eval_bws_call_arg(bws, "ydat"),
            bws = bws, ...)
  }

.np_index_regression_bw_state <- function(index.df, ydat, bws) {
  .np_semihat_make_regbw_state(
    source = bws,
    xdat = index.df,
    ydat = ydat,
    bw = bws$bw
  )
}

.np_index_formula_reentry_xdat <- function(mf) {
  terms.obj <- attr(mf, "terms")
  mf[, attr(terms.obj, "term.labels"), drop = FALSE]
}

.np_index_formula_reentry_rhs_terms <- function(formula, xdat) {
  delete.response(terms(formula, data = toFrame(xdat)))
}

.np_index_formula_reentry_response_name <- function(formula) {
  response.vars <- all.vars(formula[[2L]])
  if (length(response.vars) != 1L)
    stop("direct formula 'bws' with explicit 'txdat'/'tydat' requires a single response variable",
         call. = FALSE)
  response.vars
}

.np_index_formula_reentry_model_frame <- function(formula, txdat, tydat, call, caller_env) {
  data <- toFrame(txdat)
  if (is.data.frame(tydat) && ncol(tydat) == 1L)
    tydat <- tydat[[1L]]
  data[[.np_index_formula_reentry_response_name(formula)]] <- tydat
  data[[".np_index_formula_reentry_response"]] <- tydat

  rhs.formula <- formula(.np_index_formula_reentry_rhs_terms(formula, txdat))
  mf.formula <- as.formula(
    as.call(list(as.name("~"),
                 as.name(".np_index_formula_reentry_response"),
                 rhs.formula[[2L]])),
    env = environment(formula)
  )

  mf.call <- as.call(list(quote(stats::model.frame),
                          formula = mf.formula,
                          data = as.name(".np_index_formula_reentry_data")))
  call.names <- names(call)
  for (arg in c("subset", "na.action")) {
    pos <- match(arg, call.names, nomatch = 0L)
    if (pos > 0L)
      mf.call[[arg]] <- call[[pos]]
  }

  eval.env <- new.env(parent = caller_env)
  eval.env$.np_index_formula_reentry_data <- data
  eval(mf.call, envir = eval.env)
}

.np_index_formula_reentry_eval_xdat <- function(formula, xdat, caller_env) {
  tt <- .np_index_formula_reentry_rhs_terms(formula, xdat)
  mf.call <- as.call(list(quote(stats::model.frame),
                          formula = tt,
                          data = as.name(".np_index_formula_reentry_data")))
  eval.env <- new.env(parent = caller_env)
  eval.env$.np_index_formula_reentry_data <- toFrame(xdat)
  mf <- eval(mf.call, envir = eval.env)
  .np_index_formula_reentry_xdat(mf)
}

npindex.default <- function(bws, txdat, tydat, nomad = FALSE, ...){
  .npRmpi_require_active_slave_pool(where = "npindex()")
  sc <- sys.call()
  sc.names <- names(sc)

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  explicit.sibandwidth <- (!no.bws) && inherits(bws, "sibandwidth")
  bws.formula <- (!no.bws) && inherits(bws, "formula")
  bws.call <- (!no.bws) && is.call(bws)
  nomad <- npValidateNomadControl(nomad, "nomad")
  degree.select.value <- if (npNomadControlRequested(nomad)) {
    "coordinate"
  } else if ("degree.select" %in% names(list(...))) {
    match.arg(list(...)$degree.select, c("manual", "coordinate", "exhaustive"))
  } else {
    "manual"
  }
  if (.npRmpi_autodispatch_active() &&
      !bws.formula &&
      !bws.call &&
      (explicit.sibandwidth || identical(degree.select.value, "manual"))) {
    dispatch.call <- match.call()
    if (explicit.sibandwidth && isTRUE(.np_progress_enabled(domain = "bandwidth")))
      dispatch.call$.np_lc_fixed_progress_route <- TRUE
    return(.npRmpi_autodispatch_call(dispatch.call, parent.frame()))
  }

  if (!explicit.sibandwidth && (bws.formula || bws.call) && !no.txdat && !no.tydat) {
    txdat <- toFrame(txdat)
    if (bws.formula) {
      mf <- .np_index_formula_reentry_model_frame(
        formula = bws,
        txdat = txdat,
        tydat = tydat,
        call = sc,
        caller_env = parent.frame()
      )
      txdat <- .np_index_formula_reentry_xdat(mf)
      tydat <- model.response(mf)
    }

    dots <- list(...)
    bw.dots <- dots
    bw.dots$.np_fit_progress_handoff <- NULL
    if (bws.formula)
      bw.dots[c("subset", "na.action")] <- NULL
    bw.args <- c(list(xdat = txdat, ydat = tydat, nomad = nomad), bw.dots)
    tbw <- do.call(npindexbw, bw.args)

    fit.args <- list(bws = tbw, txdat = txdat, tydat = tydat)
    fit.dots <- dots
    fit.dots$.np_fit_progress_handoff <- TRUE
    if (bws.formula) {
      fit.dots[c("subset", "na.action")] <- NULL
      if (!is.null(fit.dots$exdat))
        fit.dots$exdat <- .np_index_formula_reentry_eval_xdat(
          formula = bws,
          xdat = fit.dots$exdat,
          caller_env = parent.frame()
        )
    }
    return(do.call(npindex, c(fit.args, fit.dots)))
  }

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")

  has.explicit.bws <- (!no.bws) && isa(bws, "sibandwidth")

  ## if bws was passed in explicitly, do not compute bandwidths

  if(txdat.named)
    txdat <- toFrame(txdat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npindexbw)

  if (bws.formula) {
    ib <- match("bws", names(sc.bw), nomatch = 0L)
    if (ib > 0L)
      names(sc.bw)[ib] <- "formula"
    drop.xy <- names(sc.bw) %in% c("txdat", "tydat")
    if (any(drop.xy))
      sc.bw <- sc.bw[!drop.xy]
  }

  if(bws.named && !bws.formula){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat')
  nstxy <- c('xdat','ydat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(!bws.formula && any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }
    
  use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
    sc.bw,
    caller_env = parent.frame()
  )

  tbw <- if (!has.explicit.bws) {
    if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting single-index bandwidth",
        .np_eval_bw_call(sc.bw, caller_env = parent.frame())
      )
    } else {
      .np_eval_bw_call(sc.bw, caller_env = parent.frame())
    }
  } else {
    bws
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
  if (no.bws || bws.formula || bws.call)
    call.args$.np_fit_progress_handoff <- TRUE
  do.call(npindex, c(call.args, list(...)))
}

.npindex_local_regression_fit <- function(source,
                                          idx.train,
                                          ydat,
                                          idx.eval = NULL,
                                          gradients = FALSE,
                                          gradient.order = 1L) {
  rbw <- .np_semihat_make_regbw_state(
    source = source,
    xdat = idx.train,
    ydat = ydat,
    bw = source$bw
  )

  direct.args <- list(
    bws = rbw,
    txdat = idx.train,
    tydat = ydat,
    exdat = idx.eval,
    gradients = gradients,
    gradient.order = gradient.order
  )

  if (isTRUE(gradients) && identical(source$type, "generalized_nn")) {
    return(.npRmpi_with_local_regression(do.call(.np_regression_direct, direct.args)))
  }

  do.call(.np_regression_direct, direct.args)
}

npindex.sibandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           eydat,
           boot.num = 399,
           errors = FALSE,
           gradients = FALSE,
           residuals = FALSE, ...) {

    fit.start <- proc.time()[3]
    dots <- list(...)
    fit.progress.handoff <- isTRUE(dots$.np_fit_progress_handoff)
    lc.fixed.progress.intent <- isTRUE(dots$.np_lc_fixed_progress_route)
    dots$.np_lc_fixed_progress_route <- NULL
    fit.progress.allow <- isTRUE(.np_progress_enabled(domain = "bandwidth"))
    gradients <- npValidateScalarLogical(gradients, "gradients")
    residuals <- npValidateScalarLogical(residuals, "residuals")
    errors <- npValidateScalarLogical(errors, "errors")
    if (!is.numeric(boot.num) || length(boot.num) != 1L || is.na(boot.num) ||
        !is.finite(boot.num) || boot.num < 1 || boot.num != floor(boot.num))
      stop("'boot.num' must be a positive integer")
    boot.num <- as.integer(boot.num)
    .npRmpi_require_active_slave_pool(where = "npindex()")
    spec <- .npindex_resolve_spec(bws, where = "npindex")
    no.ex = missing(exdat)
    no.ey = missing(eydat)

    txdat <- toFrame(txdat)
    if (!(is.vector(tydat) || is.factor(tydat)))
      stop("'tydat' must be a vector or a factor")
    tydat <-
      if (is.factor(tydat))
        as.numeric(levels(tydat))[as.integer(tydat)]
      else
        as.double(tydat)
    if (!no.ex)
      exdat <- toFrame(exdat)

    use.local.exact.exec <- .npRmpi_autodispatch_active() &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      identical(spec$regtype.engine, "lc") &&
      !identical(bws$type, "fixed") &&
      !isTRUE(gradients) &&
      !isTRUE(errors) &&
      !isTRUE(residuals) &&
      missing(eydat)
    .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = "npindex()")
    if (use.local.exact.exec) {
      .npRmpi_bcast_cmd_expr(quote(invisible(NULL)), comm = 1L, caller.execute = FALSE)

      old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
      options(npRmpi.autodispatch.disable = TRUE)
      on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)

      local.args <- c(
        list(
          bws = .npRmpi_autodispatch_untag(bws),
          txdat = txdat,
          tydat = tydat,
          boot.num = boot.num,
          errors = errors,
          gradients = gradients,
          residuals = residuals
        ),
        if (!missing(exdat)) list(exdat = exdat) else NULL,
        dots
      )
      result <- do.call(npindex, local.args)
      result <- .npRmpi_autodispatch_tag_result(result, mode = "auto")
      return(.npRmpi_restore_nomad_fit_bws_metadata(result, bws))
    }
    if (.npRmpi_autodispatch_active() &&
        !(isTRUE(gradients) && identical(bws$type, "generalized_nn"))) {
      dispatch.call <- match.call()
      if (identical(bws$method, "ichimura") &&
          identical(spec$regtype.engine, "lc") &&
          identical(bws$type, "fixed") &&
          !gradients &&
          !errors &&
          !residuals &&
          (no.ex || (!no.ex && no.ey)) &&
          (fit.progress.allow || fit.progress.handoff)) {
        dispatch.call$.np_lc_fixed_progress_route <- TRUE
      }
      result <- .npRmpi_autodispatch_call(dispatch.call, parent.frame())
      return(.npRmpi_restore_nomad_fit_bws_metadata(result, bws))
    }

    ## if no.ex then if !no.ey then ey and tx must match, to get
    ## oos errors alternatively if no.ey you get is errors if
    ## !no.ex then if !no.ey then ey and ex must match, to get
    ## oos errors alternatively if no.ey you get NO errors since we
    ## don't evaluate on the training data

    if (!no.ex){
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")

      if (!no.ey){
        if (dim(exdat)[1] != length(eydat))
          stop("number of evaluation data 'exdat' and dependent data 'eydat' do not match")
      }

    } else if(!no.ey) {
      if (dim(txdat)[1] != length(eydat))
        stop("number of training data 'txdat' and dependent data 'eydat' do not match")
    }

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(data.frame(txdat,tydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows]

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
        eydat <- adjustLevels(data.frame(eydat), bws$ydati)[,1]
        eydat <- (bws$ydati$all.dlev[[1]])[as.integer(eydat)]
      }
      else
        eydat <- as.double(eydat)
    }

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.

    txdat <- adjustLevels(txdat, bws$xdati)

    if (!no.ex)
      exdat <- adjustLevels(exdat, bws$xdati)

    ## grab the evaluation data before it is converted to numeric
    if(no.ex)
      teval <- txdat
    else
      teval <- exdat

    ## put the unordered, ordered, and continuous data in their own objects
    ## data that is not a factor is continuous.

    txdat = toMatrix(txdat)

    if (!no.ex){
      exdat = toMatrix(exdat)
    }

    ## from this point on txdat and exdat have been recast as matrices

    ## First, create the scalar index (n x 1 vector)

    index <- txdat %*% bws$beta

    if(no.ex) {
      index.eval <- index
      exdat <- txdat
      eydat <- tydat
    } else {
      index.eval <- exdat %*% bws$beta
    }
    index.df <- data.frame(index = as.vector(index))
    index.eval.df <- data.frame(index = as.vector(index.eval))

    regtype <- spec$regtype.engine
    lc.fixed.progress.route <- identical(bws$method, "ichimura") &&
      identical(regtype, "lc") &&
      identical(bws$type, "fixed") &&
      !gradients &&
      !errors &&
      !residuals &&
      (no.ex || (!no.ex && no.ey)) &&
      (fit.progress.allow || fit.progress.handoff || lc.fixed.progress.intent)
    npreg.idx.args <- list(
      txdat = if (gradients) index.df else index,
      tydat = tydat,
      bws = bws$bw,
      bwtype = bws$type,
      ckertype = bws$ckertype,
      ckerorder = bws$ckerorder,
      regtype = regtype,
      warn.glp.gradient = FALSE
    )
    if (identical(regtype, "lp")) {
      npreg.idx.args$basis <- spec$basis.engine
      npreg.idx.args$degree <- spec$degree.engine
      npreg.idx.args$bernstein.basis <- spec$bernstein.basis.engine
    }
    npreg.idx.bw <- if (identical(regtype, "lp") || lc.fixed.progress.route) {
      .np_index_regression_bw_state(
        index.df = index.df,
        ydat = tydat,
        bws = bws
      )
    } else {
      NULL
    }
    next_npreg_fit_args <- function(exdat = NULL, gradients = FALSE) {
      args <- if (identical(regtype, "lp") || lc.fixed.progress.route) {
        c(
          list(
            bws = npreg.idx.bw,
            txdat = index.df,
            tydat = tydat,
            gradients = gradients,
            warn.glp.gradient = FALSE
          ),
          if (!is.null(exdat)) list(exdat = exdat) else list()
        )
      } else {
        c(
          npreg.idx.args,
          list(gradients = gradients),
          if (!is.null(exdat)) list(exdat = exdat) else list()
        )
      }
      if (fit.progress.handoff) {
        args$.np_fit_progress_handoff <- TRUE
        fit.progress.handoff <<- FALSE
      }
      args
    }

    run_npreg_fit <- function(args) {
      if (isTRUE(lc.fixed.progress.route) ||
          isTRUE(.npRmpi_autodispatch_called_from_bcast()) ||
          isTRUE(.npRmpi_autodispatch_in_context()))
        .npRmpi_with_local_regression(do.call(npreg, args))
      else
        do.call(npreg, args)
    }

    fast.largeh <- FALSE
    fast.largeh.eval.mean <- NULL
    fast.largeh.train.mean <- NULL
    if (identical(regtype, "lc") && !gradients && !errors &&
        identical(bws$type, "fixed") && !lc.fixed.progress.route) {
      gate.index <- if (no.ex) index else c(index, index.eval)
      fast.largeh <- .npindexbw_fast_eligible(
        h = as.double(bws$bw),
        bws = bws,
        eval.index = gate.index
      )
      if (fast.largeh) {
        fast.largeh.eval.mean <- {
          tww.fast <- npksum(
            txdat = as.matrix(txdat) %*% as.matrix(bws$beta),
            tydat = as.matrix(data.frame(tydat, 1)),
            weights = as.matrix(data.frame(tydat, 1)),
            exdat = matrix(index.eval[1L], nrow = 1L),
            bws = bws$bw,
            bwtype = bws$type,
            ckertype = bws$ckertype,
            ckerorder = bws$ckerorder
          )$ksum
          as.double(tww.fast[1, 2, 1L] / NZD(tww.fast[2, 2, 1L]))
        }

        if (!no.ex && (no.ey || residuals)) {
          fast.largeh.train.mean <- {
            tww.fast <- npksum(
              txdat = as.matrix(txdat) %*% as.matrix(bws$beta),
              tydat = as.matrix(data.frame(tydat, 1)),
              weights = as.matrix(data.frame(tydat, 1)),
              exdat = matrix(index[1L], nrow = 1L),
              bws = bws$bw,
              bwtype = bws$type,
              ckertype = bws$ckertype,
              ckerorder = bws$ckerorder
            )$ksum
            as.double(tww.fast[1, 2, 1L] / NZD(tww.fast[2, 2, 1L]))
          }
        }
      }
    }

    ## Next, if no gradients are requested, use (faster) npksum

    if(gradients==FALSE) {
      if (identical(regtype, "lc") && !lc.fixed.progress.route) {
        if (fast.largeh) {
          index.mean <- rep.int(fast.largeh.eval.mean, length(index.eval))
        } else if (identical(bws$type, "fixed")) {
          tww <- npksum(txdat=as.matrix(txdat) %*% as.matrix(bws$beta),
                        tydat=as.matrix(data.frame(tydat,1)),
                        weights=as.matrix(data.frame(tydat,1)),
                        exdat=as.matrix(exdat) %*% as.matrix(bws$beta),
                        bws=bws$bw,
                        bwtype = bws$type,
                        ckertype = bws$ckertype,
                        ckerorder = bws$ckerorder)$ksum

          index.mean <- tww[1,2,]/NZD(tww[2,2,])
        } else {
          index.mean <- .np_indexhat_exact(
            bws = bws,
            idx.train = index.df,
            idx.eval = index.eval.df,
            y = tydat,
            output = "apply",
            s = 0L
          )
        }

        if (!no.ex && (no.ey || residuals)) {

          ## want to evaluate on training data for in sample errors even
          ## if evaluation x's are different from training but no y's
          ## are specified

          if (fast.largeh) {
            index.tmean <- rep.int(fast.largeh.train.mean, length(tydat))
          } else if (identical(bws$type, "fixed")) {
            tww <- npksum(txdat=as.matrix(txdat) %*% as.matrix(bws$beta),
                          tydat=as.matrix(data.frame(tydat,1)),
                          weights=as.matrix(data.frame(tydat,1)),
                          bws=bws$bw,
                          bwtype = bws$type,
                          ckertype = bws$ckertype,
                          ckerorder = bws$ckerorder)$ksum

            index.tmean <- tww[1,2,]/NZD(tww[2,2,])
          } else {
            index.tmean <- .np_indexhat_exact(
              bws = bws,
              idx.train = index.df,
              idx.eval = index.df,
              y = tydat,
              output = "apply",
              s = 0L
            )
          }

        }
      } else {
        model <- run_npreg_fit(next_npreg_fit_args(
          exdat = index.eval.df,
          gradients = FALSE
        ))
        index.mean <- model$mean

        if (!no.ex && (no.ey || residuals)) {
          model <- run_npreg_fit(next_npreg_fit_args(
            gradients = FALSE
          ))
          index.tmean <- model$mean
        }

      }

    } else if(gradients==TRUE) {
      if (identical(bws$type, "generalized_nn")) {
        model <- .npindex_local_regression_fit(
          source = bws,
          idx.train = index.df,
          ydat = tydat,
          idx.eval = index.eval.df,
          gradients = TRUE,
          gradient.order = 1L
        )
      } else {
        model <- run_npreg_fit(next_npreg_fit_args(
          exdat = index.eval.df,
          gradients = TRUE
        ))
      }

      index.mean <- model$mean

      ## index.grad is a matrix, one column for each variable, each
      ## equal to its coefficient beta_i times the first derivative of
      ## the local-constant model

      index.grad <- as.matrix(model$grad)%*%t(as.vector(bws$beta))

      if (!no.ex) {

        ## Want to evaluate on training data for in sample errors even
        ## if evaluation x's are different from training but no y's
        ## are specified. Also, needed for variance-covariance matrix
        ## (uses on ly the training data)

        model <- if (identical(bws$type, "generalized_nn")) {
          .npindex_local_regression_fit(
            source = bws,
            idx.train = index.df,
            ydat = tydat,
            gradients = TRUE,
            gradient.order = 1L
          )
        } else {
          run_npreg_fit(next_npreg_fit_args(
            gradients = TRUE
          ))
        }

        index.tmean <- model$mean

        index.tgrad <- model$grad

      }

    }

    if (no.ex) {
      index.tmean <- index.mean
    }

    if (no.ex && gradients) {
      index.tgrad <- index.grad
    }

    ## 5/3/2010, jracine, added vcov methods... thanks to Juan Carlos
    ## Escanciano <jescanci@indiana.edu> for pushing me on this for
    ## the Klein and Spady estimator... use index.tmean, index.tgrad
    ## (training X) - need gradients == TRUE in order for this to
    ## work.

    if (bws$method == "ichimura" && gradients) {

      ## First row & column of covariance matrix `Bvcov' are zero due
      ## to identification condition that beta_1=0. Note the n n^{-1}
      ## n in V^{-1}\Sigma V^{-1} and the \sqrt{n} in the
      ## normalization of \hat\beta will cancel.

      q <- ncol(txdat)
      Bvcov <- matrix(0,q,q)
      dimnames(Bvcov) <- list(bws$xnames,bws$xnames)

      ## Use the weight matrix so we can compute all expectations with
      ## only one call to npksum (the kernel arguments x\beta do not
      ## change, only the j for X_{ij} in E(X_{ij}|X_i'\beta)

      W <- txdat[,-1,drop=FALSE]

      if (identical(bws$type, "generalized_nn")) {
        kbw <- .np_indexhat_kbw(bws = bws, idx.train = index.df)
        kw <- .np_kernel_weights_direct(
          bws = kbw,
          txdat = index.df,
          bandwidth.divide = FALSE,
          kernel.pow = 1.0
        )
        if (!is.matrix(kw))
          kw <- matrix(kw, nrow = nrow(index.df))
        tyindex <- structure(as.vector(t(W) %*% kw), dim = nrow(index.df))
        tindex <- colSums(kw)
      } else {
        tyindex <- npksum(txdat = index,
                          tydat = rep(1,length(tydat)),
                          weights = W,
                          bws = bws$bw,
                          bwtype = bws$type,
                          ckertype = bws$ckertype,
                          ckerorder = bws$ckerorder)$ksum

        tindex <- npksum(txdat = index,
                         bws = bws$bw,
                         bwtype = bws$type,
                         ckertype = bws$ckertype,
                         ckerorder = bws$ckerorder)$ksum
      }

      ## Need to trap case where k-1=1... ksum will return a 1 D
      ## array, need a 1 x n matrix

      if(length(dim(tyindex))==1) tyindex <- matrix(tyindex,nrow=1,ncol=dim(tyindex))

      ## xmex = X_i-\hat E(X_i|X_i'\beta), dimension k\times n.

      xmex <- sapply(seq_along(tydat),function(i){W[i,]-tyindex[,i]/tindex[i]})

      ## Need to trap case where k-1=1..., sapply will return a
      ## vector, need a 1 x n matrix

      ## Need to trap case where k-1=1..., sapply will return a
      ## vector, need a 1 x n matrix
      if(is.vector(xmex)) {
        dg.db.xmex <- matrix(index.tgrad[,1]*xmex,nrow=1,ncol=length(xmex))
      } else {
        dg.db.xmex <- index.tgrad[,1]*xmex
      }

      uhat <- tydat - index.tmean ## Training y and training mean

      Vinv <- chol2inv(chol(dg.db.xmex%*%t(dg.db.xmex)))

      Sigma <- (uhat*dg.db.xmex)%*%t(uhat*dg.db.xmex)

      Bvcov[-1,-1] <- Vinv %*% Sigma %*% Vinv

      dimnames(Bvcov) <- list(bws$xnames,bws$xnames)

      ## Now export this in an S3 method...

    } else if (bws$method == "kleinspady" && gradients) {

      ## We divide by P(1-P) so test for P=0 or 1...

      keep <- which(index.tmean < 1 & index.tmean > 0)
      dg.db <- txdat[,-1,drop=FALSE]*index.tgrad[,1]

      ## First row & column of covariance matrix are zero due to
      ## identification condition that beta_1=0. Note the n^{-1} in
      ## the E and the \sqrt{n} in the normalization of \hat\beta will
      ## cancel.

      q <- ncol(txdat)
      Bvcov <- matrix(0,q,q)
      Bvcov[-1,-1] <- chol2inv(chol(t(dg.db[keep,])%*%(dg.db[keep,]/(index.tmean[keep]*
        (1-index.tmean[keep])))))

      dimnames(Bvcov) <- list(bws$xnames,bws$xnames)

      ## Now export this in an S3 method...

    }

    if (gradients){
      boofun = function(data, indices){
        rindex <- txdat[indices,] %*% bws$beta
        boot.args <- list(
          txdat = rindex,
          tydat = tydat[indices],
          exdat = index.eval,
          bws = bws$bw,
          ckertype = bws$ckertype,
          ckerorder = bws$ckerorder,
          regtype = regtype,
          gradients = TRUE,
          warn.glp.gradient = FALSE
        )
        if (identical(regtype, "lp")) {
          boot.args$basis <- spec$basis.engine
          boot.args$degree <- spec$degree.engine
          boot.args$bernstein.basis <- spec$bernstein.basis.engine
        }
        model <- run_npreg_fit(boot.args)[c('mean','grad')]
        
        c(model$mean, model$grad, mean(model$grad))
      }

    } else {
      boofun = function(data, indices){
        rindex = txdat[indices,] %*% bws$beta
        if (identical(regtype, "lc")) {
          tww <- npksum(txdat = rindex,
                        tydat = cbind(tydat[indices],1),
                        weights = cbind(tydat[indices],1),
                        exdat = index.eval,
                        bws = bws$bw,
                        bwtype = bws$type,
                        ckertype = bws$ckertype,
                        ckerorder = bws$ckerorder)$ksum

          tww[1,2,]/NZD(tww[2,2,])
        } else {
          boot.args <- list(
            txdat = rindex,
            tydat = tydat[indices],
            exdat = index.eval,
            bws = bws$bw,
            ckertype = bws$ckertype,
            ckerorder = bws$ckerorder,
            regtype = regtype,
            gradients = FALSE,
            warn.glp.gradient = FALSE
          )
          if (identical(regtype, "lp")) {
            boot.args$basis <- spec$basis.engine
            boot.args$degree <- spec$degree.engine
            boot.args$bernstein.basis <- spec$bernstein.basis.engine
          }
          run_npreg_fit(boot.args)$mean
        }
        
      }
    }

    if (errors){

      boot.out = suppressWarnings(boot(data.frame(txdat,tydat), boofun, R = boot.num))

      index.merr = matrix(data = 0, ncol = 1, nrow = length(index.eval))
      index.merr[,] = .np_plot_bootstrap_col_sds(boot.out$t[, seq_len(length(index.eval)), drop = FALSE])

      if (gradients) {
        index.gerr = matrix(data = 0, ncol = ncol(txdat), nrow = length(index.eval))
        index.gerr[,] = .np_plot_bootstrap_col_sds(
          boot.out$t[, (length(index.eval) + 1):(2 * length(index.eval)), drop = FALSE]
        )

        for (i in seq_len(ncol(txdat)))
          index.gerr[,i] = abs(bws$beta[i])*index.gerr[,i]

        index.mgerr = sd(boot.out$t[,2*length(index.eval)+1])
        index.mgerr = abs(bws$beta)*index.mgerr
      }
    }
    ## goodness of fit

    if(bws$method == "ichimura") {
      if (!no.ey) {
        RSQ = RSQfunc(eydat,index.mean)
        MSE = MSEfunc(eydat,index.mean)
        MAE = MAEfunc(eydat,index.mean)
        MAPE = MAPEfunc(eydat,index.mean)
        CORR = if (fast.largeh) suppressWarnings(CORRfunc(eydat,index.mean)) else CORRfunc(eydat,index.mean)
        SIGN = SIGNfunc(eydat,index.mean)
      } else {
        RSQ = RSQfunc(tydat,index.tmean)
        MSE = MSEfunc(tydat,index.tmean)
        MAE = MAEfunc(tydat,index.tmean)
        MAPE = MAPEfunc(tydat,index.tmean)
        CORR = if (fast.largeh) suppressWarnings(CORRfunc(tydat,index.tmean)) else CORRfunc(tydat,index.tmean)
        SIGN = SIGNfunc(tydat,index.tmean)
      }
      strgof = "xtra=c(RSQ,MSE,MAE,MAPE,CORR,SIGN),"
      strres = (if (residuals) "resid = tydat - index.tmean," else "")
    } else if(bws$method == "kleinspady") {
      index.pred =
        if (!no.ey) round(index.mean)
        else round(index.tmean)

      confusion.matrix =
        table(if (!no.ey) eydat else tydat,
              index.pred, dnn=c("Actual", "Predicted"))

      CCR.overall <- sum(diag(confusion.matrix))/sum(confusion.matrix)
      CCR.byoutcome <- diag(confusion.matrix)/rowSums(confusion.matrix)

      fit.mcfadden <- confusion.matrix/sum(confusion.matrix)

      fit.mcfadden <- sum(diag(fit.mcfadden)) -
        (sum(fit.mcfadden^2)-sum(diag(fit.mcfadden)^2))

      strgof = "confusion.matrix = confusion.matrix, CCR.overall = CCR.overall,
           CCR.byoutcome =  CCR.byoutcome, fit.mcfadden = fit.mcfadden,"
      strres = ""
    }

    ev.args <- list(
      bws = bws,
      index = index.eval,
      mean = index.mean,
      ntrain = nrow(txdat),
      trainiseval = no.ex,
      residuals = residuals,
      gradients = gradients
    )
    if (errors)
      ev.args$merr <- index.merr
    if (gradients) {
      ev.args$grad <- index.grad
      ev.args$mean.grad <- colMeans(index.grad)
      ev.args$betavcov <- Bvcov
    }
    if (errors && gradients) {
      ev.args$gerr <- index.gerr
      ev.args$mean.gerr <- index.mgerr
    }
    if (bws$method == "ichimura") {
      if (residuals)
        ev.args$resid <- tydat - index.tmean
      ev.args$xtra <- c(RSQ, MSE, MAE, MAPE, CORR, SIGN)
    } else if (bws$method == "kleinspady") {
      ev.args$confusion.matrix <- confusion.matrix
      ev.args$CCR.overall <- CCR.overall
      ev.args$CCR.byoutcome <- CCR.byoutcome
      ev.args$fit.mcfadden <- fit.mcfadden
    }
    ev <- do.call(singleindex, ev.args)
    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)
    ev$timing <- bws$timing
    ev$total.time <- total.time
    ev$optim.time <- optim.time
    ev$fit.time <- fit.elapsed
    ev$nomad.time <- if (!is.null(bws$nomad.time) && is.finite(bws$nomad.time)) as.double(bws$nomad.time) else NA_real_
    ev$powell.time <- if (!is.null(bws$powell.time) && is.finite(bws$powell.time)) as.double(bws$powell.time) else NA_real_
    ev
  }
