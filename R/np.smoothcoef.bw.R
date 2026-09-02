npscoefbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npscoefbw")
    .np_validate_public_dots(mc[["..."]], "npscoefbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat", "zdat"),
                                     eval_env = parent.frame())
    UseMethod("npscoefbw", target)
  }

npscoefbw.formula <-
  function(formula, data, subset, na.action, call, ...){
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
    formula.all <- if (missing(data)) {
      terms(as.formula(paste(" ~ ", paste(bronze, collapse = " + ")),
                       env = environment(formula)))
    } else {
      terms(as.formula(paste(" ~ ", paste(bronze, collapse = " + ")),
                       env = environment(formula)), data = data)
    }

    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = formula.all,
                        data = environment(formula.all),
                        eval_env = environment(formula.all))
    else .np_terms_ts_mask(terms_obj = formula.all,
                           data = data,
                           eval_env = environment(formula.all))

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

.npscoef_nn_candidate_bandwidth <- function(param, bwtype, nobs, icon = NULL) {
  if (identical(bwtype, "fixed"))
    return(as.double(param))

  if (is.null(icon)) {
    icon <- rep_len(TRUE, length(param))
  } else {
    icon <- as.logical(icon)
    if (length(icon) != length(param) || anyNA(icon))
      stop("npscoefbw: invalid nearest-neighbor coordinate map", call. = FALSE)
  }

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  hard.upper <- .Machine$integer.max / 2
  candidate <- as.double(param)
  candidate[!is.finite(candidate)] <- NA_real_
  candidate[icon] <- vapply(candidate[icon], function(h) {
    if (!is.finite(h))
      return(NA_real_)
    k <- .np_round_half_to_even(h)
    if (k < lower)
      return(NA_real_)
    if (k > upper &&
        npExtendedNnEnabled() &&
        bwtype %in% c("generalized_nn", "adaptive_nn") &&
        k <= hard.upper) {
      return(as.double(k))
    }
    if (k > upper)
      return(NA_real_)
    as.double(k)
  }, numeric(1))
  candidate
}

.npscoefbw_raw_objective_valid <- function(value,
                                           invalid.objective = sqrt(.Machine$double.xmax)) {
  is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value < invalid.objective
}

.npscoef_nn_assert_training_radius <- function(scbw, eval.zdat, owner) {
  if (!(scbw$type %in% c("generalized_nn", "adaptive_nn")) ||
      !any(scbw$icon))
    return(invisible(TRUE))

  eval.zdat <- toFrame(eval.zdat)
  if (ncol(eval.zdat) != length(scbw$icon))
    stop("internal npscoefbw NN radius check has an invalid coordinate map",
         call. = FALSE)
  continuous <- which(scbw$icon)
  bandwidth <- .npscoef_nn_candidate_bandwidth(
    param = scbw$bw,
    bwtype = scbw$type,
    nobs = nrow(eval.zdat),
    icon = scbw$icon
  )[continuous]
  if (any(!is.finite(bandwidth)))
    return(invisible(TRUE))

  tie.maximum <- vapply(continuous, function(j) {
    value <- eval.zdat[[j]]
    max(tabulate(match(value, value), nbins = length(value)))
  }, integer(1L))
  if (any(bandwidth < tie.maximum)) {
    .np_nn_abort_candidate_invalid(
      sprintf("%s produced a zero literal nearest-neighbor radius", owner),
      owner = owner,
      point = as.double(scbw$bw)
    )
  }
  invisible(TRUE)
}

.npscoef_bw_scale_multiplier <- function(scbw) {
  out <- rep.int(1.0, length(scbw$bw))
  if (!isTRUE(scbw$scaling))
    return(out)

  if (is.null(scbw$ncatfac) || !is.finite(scbw$ncatfac) ||
      is.null(scbw$nconfac) || !is.finite(scbw$nconfac))
    stop("scaled smooth-coefficient bandwidth state is missing scale factors",
         call. = FALSE)

  out[] <- as.double(scbw$ncatfac)
  if (any(scbw$icon)) {
    if (is.null(scbw$sdev) || any(!is.finite(scbw$sdev)))
      stop("scaled smooth-coefficient bandwidth state is missing continuous scales",
           call. = FALSE)
    icon.cumsum <- cumsum(scbw$icon)
    out[scbw$icon] <- as.double(scbw$nconfac) *
      as.double(scbw$sdev)[icon.cumsum[scbw$icon]]
  }
  out
}

.npscoef_apply_bw_to_scbw <- function(scbw, param, nobs = scbw$nobs) {
  param <- .npscoef_nn_candidate_bandwidth(param = param,
                                           bwtype = scbw$type,
                                           nobs = nobs,
                                           icon = scbw$icon)
  if (length(param) != length(scbw$bw))
    stop("smooth-coefficient bandwidth candidate has wrong length",
         call. = FALSE)
  scbw$bw <- as.double(param)
  if (isTRUE(scbw$scaling))
    scbw$bandwidth[[1L]] <- scbw$bw * .npscoef_bw_scale_multiplier(scbw)
  else
    scbw$bandwidth[[1L]] <- scbw$bw
  scbw
}

.npscoefbw_start_controls <- function(scale.factor.init.lower = 0.1,
                                     scale.factor.init.upper = 2.0,
                                     scale.factor.init = 0.5,
                                     lbd.init = 0.5,
                                     hbd.init = 1.5,
                                     dfac.init = 1.0,
                                     scale.factor.search.lower = 0,
                                     where = "npscoefbw") {
  cont.start <- npContinuousSearchStartControls(
    scale.factor.init.lower,
    scale.factor.init.upper,
    scale.factor.init,
    scale.factor.search.lower,
    where = where
  )
  lbd.init <- npValidatePositiveFiniteNumeric(lbd.init, "lbd.init")
  hbd.init <- npValidatePositiveFiniteNumeric(hbd.init, "hbd.init")
  dfac.init <- npValidatePositiveFiniteNumeric(dfac.init, "dfac.init")

  if (hbd.init < lbd.init) {
    stop(sprintf("%s: 'hbd.init' must be greater than or equal to 'lbd.init'", where),
         call. = FALSE)
  }
  if (lbd.init > 2 || hbd.init > 2 || dfac.init > 2) {
    stop(sprintf("%s: categorical start factors must be less than or equal to 2", where),
         call. = FALSE)
  }

  list(
    scale.factor.init.lower = cont.start$scale.factor.init.lower,
    scale.factor.init.upper = cont.start$scale.factor.init.upper,
    scale.factor.init = cont.start$scale.factor.init,
    scale.factor.search.lower = as.double(scale.factor.search.lower),
    lbd.init = as.double(lbd.init),
    hbd.init = as.double(hbd.init),
    dfac.init = as.double(dfac.init)
  )
}

.npscoef_start_factor_vector <- function(param,
                                         icon = NULL,
                                         iord = NULL,
                                         iuno = NULL,
                                         continuous.factor,
                                         categorical.factor,
                                         where = "npscoefbw") {
  ndim <- length(param)
  if (is.null(icon) || is.null(iord) || is.null(iuno))
    return(rep.int(as.double(continuous.factor), ndim))

  icon <- as.logical(icon)
  iord <- as.logical(iord)
  iuno <- as.logical(iuno)

  if (length(icon) != ndim || length(iord) != ndim || length(iuno) != ndim) {
    stop(sprintf("%s: invalid fixed-start coordinate map", where), call. = FALSE)
  }

  cat.mask <- iord | iuno
  factor <- rep_len(NA_real_, ndim)
  factor[icon] <- as.double(continuous.factor)
  factor[cat.mask] <- as.double(categorical.factor)

  if (anyNA(factor)) {
    stop(sprintf("%s: unsupported fixed-start coordinate type", where), call. = FALSE)
  }

  factor
}

.npscoef_default_start_bandwidth <- function(param,
                                             bwtype,
                                             nobs,
                                             start.controls = .npscoefbw_start_controls(),
                                             icon = NULL,
                                             iord = NULL,
                                             iuno = NULL) {
  if (identical(bwtype, "fixed"))
    return(as.double(
      param * .npscoef_start_factor_vector(
        param = param,
        icon = icon,
        iord = iord,
        iuno = iuno,
        continuous.factor = start.controls$scale.factor.init,
        categorical.factor = start.controls$dfac.init
      )
    ))

  start <- if (npExtendedNnEnabled() &&
               bwtype %in% c("generalized_nn", "adaptive_nn") &&
               is.finite(start.controls$scale.factor.init) &&
               start.controls$scale.factor.init > max(1L, as.integer(nobs) - 1L)) {
    start.controls$scale.factor.init
  } else {
    sqrt(nobs)
  }
  candidate <- as.double(
    param * .npscoef_start_factor_vector(
      param = param,
      icon = icon,
      iord = iord,
      iuno = iuno,
      continuous.factor = 1.0,
      categorical.factor = start.controls$dfac.init
    )
  )
  continuous <- if (is.null(icon)) rep_len(TRUE, length(param)) else as.logical(icon)
  candidate[continuous] <- as.double(start)
  .npscoef_nn_candidate_bandwidth(
    param = candidate,
    bwtype = bwtype,
    nobs = nobs,
    icon = icon
  )
}

.npscoef_random_start_bandwidth <- function(param,
                                            bwtype,
                                            nobs,
                                            start.controls = .npscoefbw_start_controls(),
                                            icon = NULL,
                                            iord = NULL,
                                            iuno = NULL) {
  if (identical(bwtype, "fixed")) {
    draws <- .npscoef_start_factor_vector(
      param = param,
      icon = icon,
      iord = iord,
      iuno = iuno,
      continuous.factor = start.controls$scale.factor.init.lower,
      categorical.factor = start.controls$lbd.init
    )

    for (ii in seq_along(draws)) {
      if (!is.null(icon) && !is.null(iord) && !is.null(iuno) &&
          !isTRUE(as.logical(icon)[ii])) {
        draws[ii] <- runif(1L, min = start.controls$lbd.init, max = start.controls$hbd.init)
      } else {
        draws[ii] <- runif(1L, min = start.controls$scale.factor.init.lower, max = start.controls$scale.factor.init.upper)
      }
    }

    return(as.double(draws * param))
  }

  upper <- max(1L, as.integer(nobs) - 1L)
  if (is.null(icon) || all(as.logical(icon))) {
    return(.npscoef_nn_candidate_bandwidth(
      param = runif(length(param), min = 2, max = max(2L, upper)),
      bwtype = bwtype,
      nobs = nobs,
      icon = icon
    ))
  }

  .npscoef_start_factor_vector(
    param = param,
    icon = icon,
    iord = iord,
    iuno = iuno,
    continuous.factor = 1.0,
    categorical.factor = 1.0
  )
  icon <- as.logical(icon)
  candidate <- as.double(param)
  candidate[icon] <- runif(sum(icon), min = 2, max = max(2L, upper))
  candidate[!icon] <- candidate[!icon] *
    runif(sum(!icon), min = start.controls$lbd.init,
           max = start.controls$hbd.init)
  .npscoef_nn_candidate_bandwidth(
    param = candidate,
    bwtype = bwtype,
    nobs = nobs,
    icon = icon
  )
}

.npscoef_candidate_is_admissible <- function(param,
                                            bwtype,
                                            nobs,
                                            lower = NULL,
                                            upper = NULL,
                                            icon = NULL) {
  candidate <- .npscoef_nn_candidate_bandwidth(
    param = param,
    bwtype = bwtype,
    nobs = nobs,
    icon = icon
  )
  if (any(!is.finite(candidate)))
    return(FALSE)
  if (identical(bwtype, "fixed")) {
    if (!is.null(lower))
      return(all(candidate >= lower))
    return(all(candidate > 0))
  }
  if (!is.null(icon)) {
    icon <- as.logical(icon)
    categorical <- !icon
    if (any(candidate[categorical] < 0))
      return(FALSE)
    if (!is.null(upper)) {
      if (length(upper) != length(candidate) || anyNA(upper))
        stop("npscoefbw: invalid nearest-neighbor upper-bound map", call. = FALSE)
      if (any(candidate[categorical] > upper[categorical]))
        return(FALSE)
    }
  }
  TRUE
}

.npscoef_finalize_bandwidth <- function(param,
                                        bwtype,
                                        nobs,
                                        lower = NULL,
                                        icon = NULL,
                                        where = "npscoefbw") {
  candidate <- .npscoef_nn_candidate_bandwidth(
    param = param, bwtype = bwtype, nobs = nobs, icon = icon
  )
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
  if (identical(bwtype, "fixed") && !is.null(lower) && any(candidate < lower)) {
    stop(sprintf("%s: bandwidth is below the continuous scale-factor lower bound", where),
         call. = FALSE)
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
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           lbd.init = 0.5,
           hbd.init = 1.5,
           dfac.init = 1.0,
           scale.factor.search.lower = NULL,
           ...,
           nomad.opts = list()){

    nomad.opts <- .np_nomad_normalize_user_opts(nomad.opts, "npscoefbw")
    dots <- list(...)
    if (length(nomad.opts))
      dots$nomad.opts <- nomad.opts
    npRejectUnsupportedBwsolver(dots, "npscoefbw")
    
    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)
    on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)


    miss.z <- missing(zdat)
    
    xdat <- toFrame(xdat)

    if (!miss.z)
      zdat <- toFrame(zdat)
    
    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(if (miss.z) NCOL(xdat) else NCOL(zdat))
    }
    regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
    cv.iterate <- npValidateScalarLogical(cv.iterate, "cv.iterate")
    backfit.iterate <- npValidateScalarLogical(backfit.iterate, "backfit.iterate")
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)
    backfit.maxiter <- npValidatePositiveInteger(backfit.maxiter, "backfit.maxiter")
    backfit.tol <- npValidatePositiveFiniteNumeric(backfit.tol, "backfit.tol")
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    start.controls <- .npscoefbw_start_controls(
      scale.factor.init.lower = scale.factor.init.lower,
      scale.factor.init.upper = scale.factor.init.upper,
      scale.factor.init = scale.factor.init,
      lbd.init = lbd.init,
      hbd.init = hbd.init,
      dfac.init = dfac.init,
      scale.factor.search.lower = scale.factor.search.lower,
      where = "npscoefbw"
    )
    if (cv.iterate)
      cv.num.iterations <- npValidatePositiveInteger(cv.num.iterations, "cv.num.iterations")
    spec <- .npscoef_canonical_spec(
      source = bws,
      zdat = if (miss.z) xdat else zdat,
      where = "npscoefbw"
    )
    reg.engine <- if (npIsCanonicalLp0Spec(spec, ncon = bws$ncon))
      "lc"
    else
      spec$regtype.engine

    if (!identical(reg.engine, "lc") && cv.iterate)
      stop("cv.iterate currently supports regtype='lc' for npscoefbw")

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
    apply_bw_to_scbw <- function(scbw, param) {
      scbw$sdev <- mysd
      scbw$nconfac <- nconfac
      scbw$ncatfac <- ncatfac
      .npscoef_apply_bw_to_scbw(scbw = scbw, param = param, nobs = n)
    }

    fast_largeh_tol <- npLargehRelTol()
    fast_disc_tol <- npDiscUpperRelTol()

    cont_utol <- switch(
      bws$ckertype,
      gaussian = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
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
      .npscoefbw_fast_eligible(sbw = sbw, eval.zdat = zdat.df)
    }
    objective.cache.enabled <- npObjectiveCacheEnabled()
    r.nn.cache.surface <- identical(bws$type %in% c("generalized_nn", "adaptive_nn"), TRUE) &&
      isTRUE(bws$ncon > 0L) &&
      isTRUE((bws$nuno + bws$nord) == 0L)
    r.exact.cache.surface <- identical(bws$type, "fixed") &&
      isTRUE(bws$ncon > 0L) &&
      isTRUE((bws$nuno + bws$nord) == 0L)
    r.objective.cache.kind <- if (isTRUE(r.nn.cache.surface)) {
      "nn"
    } else if (isTRUE(r.exact.cache.surface)) {
      "exact"
    } else {
      "none"
    }
    r.objective.cache.surface <- !identical(r.objective.cache.kind, "none")
    r.objective.cache.eligible <- isTRUE(bandwidth.compute) &&
      objective.cache.enabled &&
      r.objective.cache.surface
    r.objective.cache.stats <- list()
    r.objective.cache.disabled <- NULL
    if (isTRUE(bandwidth.compute) &&
        r.objective.cache.surface &&
        !objective.cache.enabled) {
      r.objective.cache.disabled <- .np_r_nn_cache_new(FALSE)
    }
    r_objective_cache_new <- function() {
      if (!r.objective.cache.eligible && is.null(r.objective.cache.disabled))
        return(NULL)
      key.length <- if (identical(r.objective.cache.kind, "nn")) bws$ncon else bws$ndim
      .np_r_nn_cache_new(r.objective.cache.eligible, key.length = key.length)
    }
    r_objective_cache_record <- function(cache) {
      st <- .np_r_nn_cache_stats(cache)
      if (!is.null(st))
        r.objective.cache.stats[[length(r.objective.cache.stats) + 1L]] <<- st
      invisible(NULL)
    }
    r_exact_cache_key <- function(x) {
      paste(sprintf("%a", as.double(x)), collapse = "\r")
    }
    r_objective_cache_lookup <- function(cache, sbw) {
      if (!is.environment(cache) || !isTRUE(cache$enabled))
        return(list(hit = FALSE, token = NULL, value = NULL))
      if (identical(r.objective.cache.kind, "nn"))
        return(.np_r_nn_cache_get(cache, as.integer(sbw$bw[sbw$icon])))
      if (identical(r.objective.cache.kind, "exact")) {
        token <- r_exact_cache_key(sbw$bw)
        return(.np_r_nn_cache_get_token(cache, token))
      }
      list(hit = FALSE, token = NULL, value = NULL)
    }
    r_objective_cache_store <- function(cache, token, value) {
      if (is.finite(value) && value < maxPenalty)
        .np_r_nn_cache_put(cache, token, value)
      invisible(NULL)
    }

    solve_cv_moment_system <- function(tyw, tww, W.eval.design, maxPenalty, Wz.eval = NULL) {
      neval <- ncol(tyw)
      ncoef <- nrow(tyw)
      pcoef <- ncol(W.eval.design)
      coef.out <- matrix(maxPenalty, nrow = pcoef, ncol = neval)
      theta.batch <- .npscoef_batch_zero_solve(tyw = tyw, tww = tww)
      if (!is.null(theta.batch)) {
        if (is.null(Wz.eval)) {
          coef.out[,] <- theta.batch
        } else {
          coef.out[,] <- .npscoef_batch_project(theta.batch, Wz.eval)
        }
        return(coef.out)
      }
      ridge.grid <- npRidgeSequenceAdditive(n.train = n, cap = 1.0)
      ridge <- rep.int(ridge.grid[1L], neval)
      ridge.idx <- rep.int(1L, neval)
      doridge <- rep.int(TRUE, neval)

      while(any(doridge)){
        iloo <- seq_len(neval)[doridge]
        for (ii in iloo) {
          doridge[ii] <- FALSE
          ridge.val <- npRidgeInterceptCorrection(
            ridge = ridge[ii], intercept = tyw[, ii][1L],
            pristine.anchor = tww[, , ii][1L, 1L]
          )
          if (is.finite(ridge.val)) {
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
          } else {
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

    cat.profile.cv <- NULL
    on.exit({ cat.profile.cv <- NULL }, add = TRUE)
    get_cat_profile_cv <- function() {
      if (is.null(cat.profile.cv)) {
        train.codes <- .np_cat_profile_code_matrix(zdat.df)
        train.keys <- .np_cat_profile_keys(train.codes)
        profile.keys <- unique(train.keys)
        train.id <- match(train.keys, profile.keys)
        train.rep <- match(profile.keys, train.keys)
        train.profile.codes <- train.codes[train.rep, , drop = FALSE]
        train.profile.dat <- zdat.df[train.rep, , drop = FALSE]
        yW.local <- cbind(ydat, W)
        p <- ncol(yW.local)
        cross.profile <- matrix(0.0, nrow = length(profile.keys), ncol = p * p)
        for (j in seq_len(p)) {
          for (k in seq_len(p)) {
            cross.profile[, (j - 1L) * p + k] <-
              .np_cat_profile_rowsum(yW.local[, j] * yW.local[, k],
                                     train.id, length(profile.keys))[, 1L]
          }
        }
        profile.index <- split(seq_along(train.id), train.id)
        profile.rows <- lapply(profile.index, function(idx) {
          Wg <- W[idx, , drop = FALSE]
          list(
            idx = idx,
            Wg = Wg,
            Wgy = Wg * ydat[idx],
            w1sq = Wg[, 1L] * Wg[, 1L]
          )
        })
        cat.profile.cv <<- list(
          train.id = train.id,
          profile.index = profile.index,
          profile.rows = profile.rows,
          train.profile.codes = train.profile.codes,
          train.profile.dat = train.profile.dat,
          G = length(profile.keys),
          yW = yW.local,
          p = p,
          cross.profile = cross.profile,
          ridge.grid = npRidgeSequenceAdditive(n.train = n, cap = 1.0)
        )
      }
      cat.profile.cv
    }

    lc_cat_profile_loo_mean <- function(sbw) {
      cp <- get_cat_profile_cv()
      train.id <- cp$train.id
      profile.index <- cp$profile.index
      profile.rows <- cp$profile.rows
      train.profile.codes <- cp$train.profile.codes
      train.profile.dat <- cp$train.profile.dat
      G <- cp$G
      yW.local <- cp$yW
      p <- cp$p
      cross.profile <- cp$cross.profile
      ridge.grid <- cp$ridge.grid

      L.profile <- .np_regression_cat_profile_kernel_matrix(
        eval.codes = train.profile.codes,
        train.codes = train.profile.codes,
        xdat = train.profile.dat,
        bws = sbw
      )

      flat.profile <- L.profile %*% cross.profile
      self.weight.profile <- L.profile[cbind(seq_len(G), seq_len(G))]
      mean.loo <- rep(maxPenalty, n)
      nc <- ncol(W)

      solve_one <- function(ii) {
        gg <- train.id[ii]
        mat.full <- matrix(flat.profile[gg, ], nrow = p, ncol = p)
        mat.full <- mat.full - self.weight.profile[gg] *
          tcrossprod(yW.local[ii, ])
        tyw.ii <- mat.full[-1L, 1L]
        tww.ii <- mat.full[-1L, -1L, drop = FALSE]

        ridge.idx <- 1L
        ridge <- ridge.grid[ridge.idx]
        repeat {
          ridge.val <- npRidgeInterceptCorrection(
            ridge = ridge, intercept = tyw.ii[1L],
            pristine.anchor = tww.ii[1L, 1L]
          )
          if (!is.finite(ridge.val))
            return(maxPenalty)
          beta.ii <- tryCatch(
            solve(tww.ii + diag(rep(ridge, nc)),
                  tyw.ii + c(ridge.val, rep(0, nc - 1L))),
            error = function(e) e
          )
          if (!inherits(beta.ii, "error")) {
            return(as.double(W[ii,, drop = FALSE] %*% beta.ii))
          }
          ridge.idx <- ridge.idx + 1L
          if (ridge.idx > length(ridge.grid))
            break
          ridge <- ridge.grid[ridge.idx]
        }
        maxPenalty
      }

      ridge <- ridge.grid[1L]
      for (gg in seq_len(G)) {
        row.state <- profile.rows[[gg]]
        idx <- row.state$idx
        mat.full <- matrix(flat.profile[gg, ], nrow = p, ncol = p)
        tyw.full <- mat.full[-1L, 1L]
        tww.full <- mat.full[-1L, -1L, drop = FALSE]
        inv.full <- tryCatch(
          solve(tww.full + diag(rep(ridge, nc))),
          error = function(e) e
        )
        if (inherits(inv.full, "error")) {
          mean.loo[idx] <- vapply(idx, solve_one, numeric(1))
          next
        }

        Wg <- row.state$Wg
        sg <- self.weight.profile[gg]
        rhs <- matrix(tyw.full, nrow = length(idx), ncol = nc, byrow = TRUE)
        rhs <- rhs - sg * row.state$Wgy
        tww11 <- tww.full[1L, 1L] - sg * row.state$w1sq
        ridge.correction <- npRidgeInterceptCorrection(
          ridge = ridge, intercept = rhs[, 1L],
          pristine.anchor = tww11
        )
        rhs[, 1L] <- rhs[, 1L] + ridge.correction

        u <- Wg %*% inv.full
        denom <- 1.0 - sg * rowSums(u * Wg)
        base <- rhs %*% inv.full
        alpha <- rowSums(u * rhs)
        beta <- base + (sg * alpha / denom) * u
        pred <- rowSums(Wg * beta)
        bad <- !is.finite(ridge.correction) | !is.finite(pred) |
          !is.finite(denom) |
          abs(denom) < sqrt(.Machine$double.eps)
        if (any(bad))
          pred[bad] <- vapply(idx[bad], solve_one, numeric(1))
        mean.loo[idx] <- pred
      }

      mean.loo
    }

    lc_cat_profile_partial_sums <- function(wj, partial.y) {
      cp <- get_cat_profile_cv()
      list(
        num.profile = .np_cat_profile_rowsum(partial.y * wj, cp$train.id, cp$G)[, 1L],
        den.profile = .np_cat_profile_rowsum(wj * wj, cp$train.id, cp$G)[, 1L]
      )
    }

    lc_cat_profile_partial_loo <- function(sbw, wj, partial.y, profile.sums = NULL) {
      cp <- get_cat_profile_cv()
      train.id <- cp$train.id
      train.profile.codes <- cp$train.profile.codes
      train.profile.dat <- cp$train.profile.dat
      G <- cp$G

      L.profile <- .np_regression_cat_profile_kernel_matrix(
        eval.codes = train.profile.codes,
        train.codes = train.profile.codes,
        xdat = train.profile.dat,
        bws = sbw
      )

      if (is.null(profile.sums))
        profile.sums <- lc_cat_profile_partial_sums(wj = wj, partial.y = partial.y)
      num.profile <- profile.sums$num.profile
      den.profile <- profile.sums$den.profile
      num <- as.vector(L.profile %*% num.profile)[train.id]
      den <- as.vector(L.profile %*% den.profile)[train.id]
      self.weight <- L.profile[cbind(seq_len(G), seq_len(G))][train.id]
      num <- num - self.weight * partial.y * wj
      den <- den - self.weight * wj * wj

      wj * num / NZD(den)
    }

    use_cat_profile_cv_lc <- function(sbw) {
      identical(reg.engine, "lc") &&
        identical(sbw$type, "fixed") &&
        npUseCategoricalCompress(ncon = sbw$ncon,
                                 ncat = sbw$nuno + sbw$nord) &&
        !miss.z &&
        isTRUE(sbw$ncon == 0L) &&
        isTRUE((sbw$nuno + sbw$nord) > 0L)
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
      main.ks <- do.call(.npscoef_npksum, ksum.args)$ksum
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
      main.ks <- do.call(.npscoef_npksum, ksum.args)$ksum
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

          partial_progress_detail <- function() {
            sprintf(
              "backfitting iteration %d of %d, partial residual %d of %d",
              cv_state$backfit_iteration,
              cv.num.iterations,
              cv_state$partial_index,
              ncol(W)
            )
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
              detail = partial_progress_detail()
            )
            invisible(NULL)
          }

          partial_progress_finish <- function(fv = NULL) {
            if (is.null(cv_state$partial_progress))
              return(invisible(NULL))

            cv_state$partial_progress$last_emit <- -Inf
            cv_state$partial_progress <- .np_progress_end(
              cv_state$partial_progress,
              detail = partial_progress_detail()
            )
            cv_state$partial_progress <- NULL
            invisible(NULL)
          }

          overall.cache <- NULL
          overall.cv.ls.raw <- function(param, account = TRUE) {
            if (account)
              cv_state$objective_fast <- FALSE
            sbw <- apply_bw_to_scbw(bws, param)
            if (!validateBandwidthTF(sbw) ||
                (!is.null(fixed.lower) && any(param < fixed.lower)) ||
                ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            .npscoef_nn_assert_training_radius(
              scbw = sbw,
              eval.zdat = zdat,
              owner = "npscoefbw overall CV objective"
            )
            cache.hit <- if (account) {
              r_objective_cache_lookup(overall.cache, sbw)
            } else {
              list(hit = FALSE, token = NULL)
            }
            if (isTRUE(cache.hit$hit)) {
              cv_progress_step()
              cv_state$fast_total <- cv_state$fast_total + 1L
              return(cache.hit$value)
            }
            objective.fast <- npscoef_fast_eligible(sbw) ||
              use_cat_profile_cv_lc(sbw)
            if (account)
              cv_state$objective_fast <- objective.fast

            if (identical(reg.engine, "lc")) {
              if (use_cat_profile_cv_lc(sbw)) {
                mean.loo <- lc_cat_profile_loo_mean(sbw)
              } else {
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
                    ridge.val <- npRidgeInterceptCorrection(
                      ridge = ridge[ii], intercept = tww[-1, 1, ii][1L],
                      pristine.anchor = tww[-1, -1, ii][1L, 1L]
                    )
                    if (!is.finite(ridge.val)) {
                      mean.loo[ii] <- NA_real_
                      next
                    }
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
                        next
                      }
                      mean.loo[ii] <- NA_real_
                      next
                    }
                    mean.loo[ii] <- W[ii,, drop = FALSE] %*% beta.ii
                  }
                }
              }
            } else {
              coef.loo <- lp_full_coef(sbw = sbw, leave.one.out.eval = TRUE)
              mean.loo <- rowSums(W * t(coef.loo))
            }

            if (!all(is.finite(mean.loo))) {
              if (account)
                cv_progress_step(ridging = TRUE)
              if (account && isTRUE(objective.fast))
                cv_state$fast_total <- cv_state$fast_total + 1L
              if (account)
                r_objective_cache_store(overall.cache, cache.hit$token, maxPenalty)
              return(maxPenalty)
            }

            if(!any(mean.loo == maxPenalty)){
              fv <- sum((ydat-mean.loo)^2)/n
              if (account)
                cv_progress_step()
            } else {
              if (account)
                cv_progress_step(ridging = TRUE)
              fv <- maxPenalty
            }

            if (account && isTRUE(objective.fast))
              cv_state$fast_total <- cv_state$fast_total + 1L
            if (account)
              r_objective_cache_store(overall.cache, cache.hit$token, fv)

            return((if (is.finite(fv)) fv else maxPenalty))

          }

          overall.cv.ls <- function(param) {
            tryCatch(
              overall.cv.ls.raw(param),
              np_nn_candidate_invalid = function(e) {
                cv_progress_step(ridging = TRUE)
                maxPenalty
              }
            )
          }

          scoef.loo.args <- list(
            bws = bws, txdat = xdat, tydat = ydat,
            leave.one.out = TRUE, iterate = TRUE,
            maxiter = backfit.maxiter, tol = backfit.tol,
            betas = TRUE, se = FALSE,
            .np_fit_progress_allow = FALSE
          )
          if (!miss.z)
            scoef.loo.args$tzdat <- zdat

          current.partial.profile <- NULL
          current.partial.cache <- NULL
          
          partial.cv.ls.raw <- function(param, partial.index) {
            cv_state$objective_fast <- FALSE
            sbw <- apply_bw_to_scbw(bws, param)

            if (!validateBandwidthTF(sbw) ||
                (!is.null(fixed.lower) && any(param < fixed.lower)) ||
                ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            .npscoef_nn_assert_training_radius(
              scbw = sbw,
              eval.zdat = zdat,
              owner = "npscoefbw partial CV objective"
            )
            cache.hit <- r_objective_cache_lookup(current.partial.cache, sbw)
            if (isTRUE(cache.hit$hit)) {
              cv_state$fast_total <- cv_state$fast_total + 1L
              partial_progress_step(fv = cache.hit$value)
              return(cache.hit$value)
            }
            cv_state$objective_fast <- npscoef_fast_eligible(sbw) ||
              use_cat_profile_cv_lc(sbw)
            
            if (backfit.iterate){
              local.bws <- bws
              if (is.null(local.bws$bw.fitted)) {
                local.bws$bw.fitted <- matrix(
                  local.bws$bw,
                  nrow = length(local.bws$bw),
                  ncol = n.part
                )
              }
              local.bws$bw.fitted[, partial.index] <- sbw$bw
              local.args <- scoef.loo.args
              local.args$bws <- local.bws
              scoef.loo <- do.call(npscoef, local.args)
              partial.loo <- W[,partial.index]*scoef.loo$beta[,partial.index]
            } else {
              wj <- W[,partial.index]
              if (identical(reg.engine, "lc")) {
                if (use_cat_profile_cv_lc(sbw)) {
                  partial.loo <- lc_cat_profile_partial_loo(
                    sbw = sbw,
                    wj = wj,
                    partial.y = partial.orig,
                    profile.sums = current.partial.profile
                  )
                } else {
                  tww <- npksum(txdat=zdat,
                                tydat=cbind(partial.orig * wj, wj * wj),
                                weights=cbind(partial.orig * wj, 1),
                                bws=sbw,
                                leave.one.out=TRUE)$ksum

                  partial.loo <- wj * tww[1,2,]/NZD(tww[2,2,])
                }
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
            r_objective_cache_store(current.partial.cache, cache.hit$token, fv)

            partial_progress_step(fv = fv)
            return((if (is.finite(fv)) fv else maxPenalty))
          }

          partial.cv.ls <- function(param, partial.index) {
            tryCatch(
              partial.cv.ls.raw(param, partial.index),
              np_nn_candidate_invalid = function(e) {
                partial_progress_step(fv = maxPenalty)
                maxPenalty
              }
            )
          }

          ## Now we implement multistarting

          fval.min <- .Machine$double.xmax
          have_best <- FALSE
          numimp <- 0
          value.overall <- numeric(nmulti)
          num.feval.overall <- 0
          overall.cache <- r_objective_cache_new()
          automatic.start <- all(bws$bw == 0)
          first.automatic.start <- NULL

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
          fixed.lower <- if (identical(bws$type, "fixed")) {
            out <- rep.int(0, length(x.scale))
            out[dati$icon] <- x.scale[dati$icon] * start.controls$scale.factor.search.lower
            out
          } else {
            NULL
          }
          candidate.upper <- rep.int(Inf, length(x.scale))
          candidate.upper[!dati$icon] <- 2.0*x.scale[!dati$icon]

          optim.control <- list(abstol = optim.abstol,
                                reltol = optim.reltol,
                                maxit = optim.maxit)
          configured.optim.control <- optim.control

          for (i in seq_len(nmulti)) {
            cv_state$multistart_index <- i
            cv_progress_begin()
            
            if (i == 1) {
              tbw <- .npscoef_default_start_bandwidth(
                param = x.scale,
                bwtype = bws$type,
                nobs = n,
                start.controls = start.controls,
                icon = dati$icon,
                iord = dati$iord,
                iuno = dati$iuno
              )
              if (all(bws$bw != 0) &&
                  .npscoef_candidate_is_admissible(param = bws$bw, bwtype = bws$type, nobs = n,
                                                   lower = fixed.lower,
                                                   upper = candidate.upper,
                                                   icon = dati$icon)) {
                tbw <- .npscoef_finalize_bandwidth(
                  param = bws$bw,
                  bwtype = bws$type,
                  nobs = n,
                  lower = fixed.lower,
                  icon = dati$icon,
                  where = "npscoefbw"
                )
              }
              if (automatic.start)
                first.automatic.start <- tbw
            } else {
              tbw <- .npscoef_random_start_bandwidth(
                param = x.scale,
                bwtype = bws$type,
                nobs = n,
                start.controls = start.controls,
                icon = dati$icon,
                iord = dati$iord,
                iuno = dati$iuno
              )
            }

            suppressWarnings(optim.return <- optim(tbw,
                                                   fn = overall.cv.ls,
                                                   method = optim.method,
                                                   control = optim.control))
            if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
              num.feval.overall <- num.feval.overall + optim.return$counts[1]
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              tbw <- .npscoef_random_start_bandwidth(
                param = x.scale,
                bwtype = bws$type,
                nobs = n,
                start.controls = start.controls,
                icon = dati$icon,
                iord = dati$iord,
                iuno = dati$iuno
              )
              optim.control <- lapply(optim.control, '*', 10.0)
              suppressWarnings(optim.return <- optim(tbw,
                                                     fn = overall.cv.ls,
                                                     method = optim.method,
                                                     control = optim.control))
              if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                num.feval.overall <- num.feval.overall + optim.return$counts[1]

            }

            cv_progress_finish()

            value.overall[i] <- optim.return$value

            if (.npscoefbw_raw_objective_valid(
                  optim.return$value,
                  invalid.objective = maxPenalty
                ) && .npscoef_candidate_is_admissible(
              param = optim.return$par,
              bwtype = bws$type,
              nobs = n,
              lower = fixed.lower,
              upper = candidate.upper,
              icon = dati$icon
            ) && (!have_best || optim.return$value < fval.min)) {
              param <- .npscoef_finalize_bandwidth(
                param = optim.return$par,
                bwtype = bws$type,
                nobs = n,
                lower = fixed.lower,
                icon = dati$icon,
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

          if (automatic.start &&
              bws$type %in% c("generalized_nn", "adaptive_nn") &&
              (!have_best || !is.finite(fval.min) || fval.min >= maxPenalty)) {
            cv_state$multistart_index <- if (have_best) best.overall else 1L
            cv_progress_begin()
            recovery.raw.eval <- function(point) {
              value <- overall.cv.ls.raw(point)
              if (!.npscoefbw_raw_objective_valid(
                    value,
                    invalid.objective = maxPenalty
                  )) {
                .np_nn_abort_candidate_invalid(
                  "npscoefbw overall CV objective returned an invalid raw objective",
                  owner = "npscoefbw overall CV objective",
                  point = point,
                  raw.objective = value
                )
              }
              as.double(value)
            }
            ordinary.caps <- rep.int(n - 1L, sum(dati$icon))
            ordinary.start <- first.automatic.start[dati$icon]
            recovery <- if (length(ordinary.start) &&
                            all(ordinary.caps >= 2L) &&
                            all(ordinary.start >= 2L) &&
                            all(ordinary.start <= ordinary.caps)) {
              .np_nn_find_raw_valid_start(
                point = first.automatic.start,
                nn.indices = which(dati$icon),
                caps = ordinary.caps,
                raw.eval = recovery.raw.eval
              )
            } else {
              list(found = FALSE, point = NULL, objective = NA_real_, evaluations = 0L)
            }
            num.feval.overall <- num.feval.overall + recovery$evaluations

            if (isTRUE(recovery$found)) {
              suppressWarnings(recovery.optim <- optim(
                recovery$point,
                fn = overall.cv.ls,
                method = optim.method,
                control = configured.optim.control
              ))
              if (!is.null(recovery.optim$counts) && length(recovery.optim$counts) > 0L)
                num.feval.overall <- num.feval.overall + recovery.optim$counts[1L]

              final.raw <- tryCatch(
                recovery.raw.eval(recovery.optim$par),
                np_nn_candidate_invalid = identity
              )
              num.feval.overall <- num.feval.overall + 1L
              if (!inherits(final.raw, "np_nn_candidate_invalid") &&
                  .npscoef_candidate_is_admissible(
                    param = recovery.optim$par,
                    bwtype = bws$type,
                    nobs = n,
                    upper = candidate.upper,
                    icon = dati$icon
                  )) {
                param <- .npscoef_finalize_bandwidth(
                  param = recovery.optim$par,
                  bwtype = bws$type,
                  nobs = n,
                  icon = dati$icon,
                  where = "npscoefbw"
                )
                min.overall <- as.double(final.raw)
                fval.min <- min.overall
                best.overall <- cv_state$multistart_index
                value.overall[best.overall] <- min.overall
                numimp.overall <- if (exists("numimp.overall", inherits = FALSE))
                  numimp.overall + 1L else 1L
                have_best <- TRUE
              }
            }
            cv_progress_finish()
          }
          r_objective_cache_record(overall.cache)

          if (!have_best) {
            if (identical(bws$type, "fixed")) {
              stop("npscoefbw: no feasible fixed bandwidths found", call. = FALSE)
            }
            stop("npscoefbw: no feasible bandwidths found", call. = FALSE)
          }
          if (!is.finite(fval.min) || fval.min >= maxPenalty)
            stop("npscoefbw: optimizer returned a bandwidth candidate with invalid objective",
                 call. = FALSE)

          if (bws$type %in% c("generalized_nn", "adaptive_nn") && !cv.iterate) {
            final.raw <- overall.cv.ls.raw(param, account = FALSE)
            if (!.npscoefbw_raw_objective_valid(final.raw, maxPenalty) ||
                !identical(as.numeric(final.raw), as.numeric(fval.min))) {
              stop("npscoefbw: optimizer endpoint failed raw-objective certification",
                   call. = FALSE)
            }
          }

          param.overall <- bws$bw <- .npscoef_finalize_bandwidth(
            param = param,
            bwtype = bws$type,
            nobs = n,
            icon = dati$icon,
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

            scoef.args <- list(
              bws = bws,
              txdat = xdat,
              tydat = ydat,
              iterate = FALSE,
              betas = TRUE,
              se = FALSE,
              .np_fit_progress_allow = FALSE
            )
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
                current.partial.profile <- if (identical(reg.engine, "lc") &&
                    npUseCategoricalCompress(ncon = bws$ncon,
                                             ncat = bws$nuno + bws$nord) &&
                    !miss.z &&
                    isTRUE(bws$ncon == 0L) &&
                    isTRUE((bws$nuno + bws$nord) > 0L)) {
                  lc_cat_profile_partial_sums(wj = W[, j], partial.y = partial.orig)
                } else {
                  NULL
                }
                partial_progress_begin(iteration = i, partial.index = j)
                current.partial.cache <- r_objective_cache_new()
                
                ## minimise
                partial.start <- bws$bw.fitted[, j]
                suppressWarnings(optim.return <-
                                 optim(partial.start, fn = partial.cv.ls,
                                       method = optim.method,
                                       control = optim.control,
                                       partial.index = j))
                if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                  num.feval.overall <- num.feval.overall + optim.return$counts[1]
                partial_progress_finish(fv = optim.return$value)
                r_objective_cache_record(current.partial.cache)
                current.partial.cache <- NULL
                current.partial.profile <- NULL
                
                ## grab parameter
                bws$bw.fitted[,j] <- optim.return$par

                if (backfit.iterate){
                  ## re-estimate all betas
                  scoef.args <- list(
                    bws = bws, txdat = xdat, tydat = ydat,
                    iterate = TRUE, maxiter = backfit.maxiter,
                    tol = backfit.tol, betas = TRUE, se = FALSE,
                    .np_fit_progress_allow = FALSE
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
              tol = backfit.tol, leave.one.out = TRUE,
              se = FALSE,
              .np_fit_progress_allow = FALSE
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
          if (length(r.objective.cache.stats)) {
            bws$nn.cache <- .np_r_nn_cache_combine_stats(r.objective.cache.stats)
          } else {
            bws$nn.cache <- .np_r_nn_cache_stats(r.objective.cache.disabled)
          }
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

    .np_seed_exit(seed.state, remove_if_absent = TRUE)
    nn.cache <- bws$nn.cache
    
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
                       bw.fitted = bws$bw.fitted,
                       fval.fitted = bws$fval.fitted,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       zdati = bws$zdati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       znames = bws$znames,
                       sfactor = bws$sfactor,
                       bandwidth = bws$bandwidth,
                       sdev = bws$sdev,
                       nconfac = bws$nconfac,
                       ncatfac = bws$ncatfac,
                       rows.omit = rows.omit,
                       bandwidth.compute = bandwidth.compute,
                       optim.method = optim.method,
                       total.time = total.time)
    bws$nn.cache <- nn.cache
    bws <- npSetScaleFactorSearchLower(bws, scale.factor.search.lower)

    bws
  }

.npscoefbw_build_scbandwidth <- function(xdat,
                                         ydat,
                                         zdat,
                                         bws,
                                         bandwidth.compute,
                                         reg.args) {
  miss.z <- is.null(zdat)
  zdati <- if (miss.z) NULL else untangle(zdat)
  znames <- if (miss.z) NULL else names(zdat)

  sbw.args <- c(
    list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      zdati = zdati,
      xnames = names(xdat),
      ynames = deparse(substitute(ydat)),
      znames = znames,
      bandwidth.compute = bandwidth.compute
    ),
    reg.args
  )

  out <- do.call(scbandwidth, sbw.args)
  if (!is.null(reg.args$scale.factor.search.lower))
    out$scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      reg.args$scale.factor.search.lower
    )
  out
}

.npscoefbw_run_fixed_degree <- function(xdat, ydat, zdat, bws, reg.args, opt.args) {
  tbw <- .npscoefbw_build_scbandwidth(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args
  )

  scbw.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
  if (!is.null(zdat))
    scbw.args$zdat <- zdat
  do.call(npscoefbw.scbandwidth, scbw.args)
}

.npscoefbw_nomad_controls <- function(search.engine) {
  .np_degree_search_engine_controls(search.engine)
}

.npscoefbw_fast_eligible <- function(sbw, eval.zdat) {
  if (!identical(sbw$type, "fixed"))
    return(FALSE)

  tdati <- if (is.null(sbw$zdati)) sbw$xdati else sbw$zdati
  eval.zdat <- toFrame(eval.zdat)

  ckerbound <- if (is.null(sbw$ckerbound) || !length(sbw$ckerbound)) {
    "none"
  } else {
    as.character(sbw$ckerbound)[1L]
  }
  if (!identical(ckerbound, "none"))
    return(FALSE)

  ckertype <- as.character(sbw$ckertype)[1L]
  ckerorder <- as.integer(sbw$ckerorder)[1L]
  if (!identical(ckertype, "uniform") &&
      (!is.finite(ckerorder) || ckerorder != 2L))
    return(FALSE)

  if (any(tdati$icon) && !npLogicalOption("np.largeh", TRUE))
    return(FALSE)
  if ((any(tdati$iuno) || any(tdati$iord)) &&
      !npLogicalOption("np.largelambda", TRUE))
    return(FALSE)

  fast_largeh_tol <- npLargehRelTol()
  fast_disc_tol <- npDiscUpperRelTol()

  cont_utol <- switch(
    ckertype,
    gaussian = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
    epanechnikov = sqrt(fast_largeh_tol),
    uniform = 1.0 - 32.0 * .Machine$double.eps,
    0.0
  )

  cont_hmin <- numeric(0)
  if (any(sbw$icon) && is.finite(cont_utol) && cont_utol > 0) {
    zcon <- eval.zdat[, sbw$icon, drop = FALSE]
    cont_hmin <- vapply(zcon, function(col) {
      vals <- as.double(col)
      if (!length(vals) || any(!is.finite(vals)))
        return(Inf)
      diff(range(vals)) / cont_utol
    }, numeric(1))
  }

  disc_upper_tol <- function(upper) {
    max(fast_disc_tol * abs(upper),
        16.0 * .Machine$double.eps * max(1.0, abs(upper)))
  }

  uno_upper <- numeric(0)
  if (any(sbw$iuno)) {
    uno_idx <- which(sbw$iuno)
    uno_upper <- vapply(uno_idx, function(i) {
      uMaxL(tdati$all.nlev[[i]], kertype = sbw$ukertype)
    }, numeric(1))
  }

  ord_upper <- numeric(0)
  if (any(sbw$iord)) {
    ord_idx <- which(sbw$iord)
    ord_upper <- vapply(ord_idx, function(i) {
      oMaxL(tdati$all.nlev[[i]], kertype = sbw$okertype)
    }, numeric(1))
  }

  bwv <- sbw$bandwidth[[1L]]
  if (!length(bwv) || length(bwv) != length(sbw$icon))
    return(FALSE)

  if (any(sbw$icon)) {
    bw_cont <- bwv[sbw$icon]
    if (any(!is.finite(bw_cont)) || any(bw_cont <= 0) ||
        any(bw_cont < cont_hmin))
      return(FALSE)
  }

  if (any(sbw$iuno)) {
    bw_uno <- bwv[sbw$iuno]
    ok_uno <- mapply(function(bw, upper) {
      is.finite(bw) && abs(bw - upper) <= disc_upper_tol(upper)
    }, bw = bw_uno, upper = uno_upper, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    if (!all(ok_uno))
      return(FALSE)
  }

  if (any(sbw$iord)) {
    bw_ord <- bwv[sbw$iord]
    ok_ord <- mapply(function(bw, upper) {
      is.finite(bw) && abs(bw - upper) <= disc_upper_tol(upper)
    }, bw = bw_ord, upper = ord_upper, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    if (!all(ok_ord))
      return(FALSE)
  }

  TRUE
}

.npscoefbw_nomad_context_prepare <- function(xdat, ydat, zdat = NULL) {
  miss.z <- is.null(zdat)
  xdat <- toFrame(xdat)

  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector or a factor")

  if (!miss.z)
    zdat <- toFrame(zdat)

  keep.rows <- rep_len(TRUE, nrow(xdat))
  train.df <- data.frame(xdat, ydat)
  if (!miss.z)
    train.df <- data.frame(train.df, zdat)
  rows.omit <- attr(na.omit(train.df), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE

  if (!any(keep.rows))
    stop("Data has no rows without NAs")

  xdat <- xdat[keep.rows, , drop = FALSE]
  ydat <- ydat[keep.rows]
  if (!miss.z)
    zdat <- zdat[keep.rows, , drop = FALSE]

  yvec <- if (is.factor(ydat)) {
    dlev(ydat)[as.integer(ydat)]
  } else {
    as.double(ydat)
  }

  xmat <- toMatrix(xdat)
  if (qr(xmat)$rank < ncol(xmat))
    stop("columns of the independent variable (xdat) are linearly dependent")

  zdf <- if (miss.z) xdat else if (is.data.frame(zdat)) zdat else as.data.frame(zdat)

  list(
    xdat = xdat,
    xmat = xmat,
    ydat = yvec,
    zdat.df = zdf,
    W = cbind(1.0, xmat),
    n = nrow(xmat)
  )
}

.npscoefbw_normalize_nomad_scbw <- function(scbw, eval.zdat, bw = scbw$bw) {
  rbw <- .npscoef_make_regbw(
    bws = scbw,
    zdat = eval.zdat,
    bw = bw
  )
  scbw$bw <- rbw$bw
  scbw$bandwidth[[1L]] <- rbw$bandwidth[[1L]]
  scbw$sfactor[[1L]] <- rbw$sfactor[[1L]]
  scbw$nconfac <- rbw$nconfac
  scbw$ncatfac <- rbw$ncatfac
  scbw$sdev <- rbw$sdev
  scbw
}

.npscoefbw_nomad_solve_cv_moment_system <- function(tyw,
                                                    tww,
                                                    W.eval.design,
                                                    maxPenalty,
                                                    n.train,
                                                    Wz.eval = NULL) {
  neval <- ncol(tyw)
  ncoef <- nrow(tyw)
  pcoef <- ncol(W.eval.design)
  coef.out <- matrix(maxPenalty, nrow = pcoef, ncol = neval)
  theta.batch <- .npscoef_batch_zero_solve(tyw = tyw, tww = tww)
  if (!is.null(theta.batch)) {
    if (is.null(Wz.eval)) {
      coef.out[,] <- theta.batch
    } else {
      coef.out[,] <- .npscoef_batch_project(theta.batch, Wz.eval)
    }
    return(coef.out)
  }
  ridge.grid <- npRidgeSequenceAdditive(n.train = n.train, cap = 1.0)
  ridge <- rep.int(ridge.grid[1L], neval)
  ridge.idx <- rep.int(1L, neval)
  doridge <- rep.int(TRUE, neval)

  while (any(doridge)) {
    iloo <- seq_len(neval)[doridge]
    for (ii in iloo) {
      doridge[ii] <- FALSE
      ridge.val <- npRidgeInterceptCorrection(
        ridge = ridge[ii], intercept = tyw[, ii][1L],
        pristine.anchor = tww[, , ii][1L, 1L]
      )
      if (is.finite(ridge.val)) {
        theta.ii <- tryCatch(
          solve(
            tww[, , ii] + diag(rep(ridge[ii], ncoef)),
            tyw[, ii] + c(ridge.val, rep(0, ncoef - 1L))
          ),
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
      } else {
        theta.ii <- rep(maxPenalty, ncoef)
      }

      if (is.null(Wz.eval)) {
        coef.out[, ii] <- theta.ii
      } else {
        coef.out[, ii] <- as.vector(crossprod(
          Wz.eval[ii, ],
          matrix(theta.ii, nrow = ncol(Wz.eval), ncol = pcoef)
        ))
      }
    }
  }

  coef.out
}

.npscoefbw_nomad_moment_state <- function(ctx,
                                         bws,
                                         leave.one.out,
                                         eval.idx = NULL) {
  if (is.null(eval.idx))
    eval.idx <- seq_len(ctx$n)
  eval.zdat <- ctx$zdat.df[as.integer(eval.idx), , drop = FALSE]
  spec <- .npscoef_canonical_spec(
    source = bws,
    zdat = ctx$zdat.df,
    where = "npscoefbw"
  )
  if (npIsCanonicalLp0Spec(spec, ncon = bws$ncon)) {
    z.train <- adjustLevels(ctx$zdat.df, bws$zdati)
    return(list(
      z.train = z.train,
      z.eval = if (leave.one.out) z.train else
        adjustLevels(eval.zdat, bws$zdati, allowNewCells = TRUE),
      kernel.bws = bws,
      ytensor = cbind(ctx$ydat, ctx$W),
      Wz.eval = NULL,
      bandwidth.divide = FALSE
    ))
  }

  lp.state <- .npscoef_lp_state(
    bws = bws,
    tzdat = ctx$zdat.df,
    ezdat = eval.zdat,
    leave.one.out = leave.one.out,
    where = "npscoefbw"
  )
  tensor.train <- .npscoef_row_tensor_design(ctx$W, lp.state$W.train)
  list(
    z.train = lp.state$z.train,
    z.eval = lp.state$z.eval,
    kernel.bws = lp.state$rbw,
    ytensor = cbind(ctx$ydat, tensor.train),
    Wz.eval = lp.state$W.eval,
    bandwidth.divide = TRUE
  )
}

.npscoefbw_nomad_invalid_objective <- function(bws,
                                               invalid.penalty = c("baseline", "large"),
                                               penalty.multiplier = 10,
                                               invalid.objective = NULL) {
  invalid.penalty <- match.arg(invalid.penalty)
  base.penalty <- switch(
    invalid.penalty,
    baseline = if (is.finite(bws$fval[1L])) as.numeric(bws$fval[1L]) else 1,
    large = 1
  )
  base.penalty <- max(abs(base.penalty), 1)
  if (is.null(invalid.objective)) {
    penalty.multiplier * base.penalty
  } else {
    invalid.objective <- as.numeric(invalid.objective)[1L]
    if (!is.finite(invalid.objective) || invalid.objective <= 0)
      stop("invalid NOMAD smooth-coefficient objective penalty")
    invalid.objective
  }
}

.npscoefbw_nomad_eval_direct <- function(ctx,
                                         bws,
                                         invalid.penalty = c("baseline", "large"),
                                         penalty.multiplier = 10,
                                         invalid.objective = NULL) {
  penalty <- .npscoefbw_nomad_invalid_objective(
    bws = bws,
    invalid.penalty = invalid.penalty,
    penalty.multiplier = penalty.multiplier,
    invalid.objective = invalid.objective
  )

  if (!is.list(ctx) || is.null(ctx$W) || is.null(ctx$ydat) || is.null(ctx$zdat.df))
    stop("invalid NOMAD smooth-coefficient context")

  if (!validateBandwidthTF(bws)) {
    return(list(
      objective = penalty,
      num.feval = 1L,
      num.feval.fast = 0L,
      raw.valid = FALSE
    ))
  }

  maxPenalty <- sqrt(.Machine$double.xmax)

  tryCatch({
    .npscoef_nn_assert_training_radius(
      scbw = bws,
      eval.zdat = ctx$zdat.df,
      owner = "npscoefbw NOMAD degree objective"
    )
    moment.state <- .npscoefbw_nomad_moment_state(
      ctx = ctx,
      bws = bws,
      leave.one.out = TRUE
    )
    main.ks <- .npscoef_npksum(
      txdat = moment.state$z.train,
      tydat = moment.state$ytensor,
      weights = moment.state$ytensor,
      bws = moment.state$kernel.bws,
      leave.one.out = TRUE,
      bandwidth.divide = moment.state$bandwidth.divide
    )$ksum
    tyw <- main.ks[-1L, 1L, , drop = FALSE]
    if (length(dim(tyw)) == 3L)
      dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
    tww <- main.ks[-1L, -1L, , drop = FALSE]
    coef.loo <- .npscoefbw_nomad_solve_cv_moment_system(
      tyw = tyw,
      tww = tww,
      W.eval.design = ctx$W,
      Wz.eval = moment.state$Wz.eval,
      maxPenalty = maxPenalty,
      n.train = ctx$n
    )
    mean.loo <- rowSums(ctx$W * t(coef.loo))

    if (any(!is.finite(mean.loo)) || any(mean.loo == maxPenalty)) {
      list(objective = penalty, num.feval = 1L, num.feval.fast = 0L,
           raw.valid = FALSE)
    } else {
      objective <- as.numeric(mean((ctx$ydat - mean.loo)^2))
      raw.valid <- .npscoefbw_raw_objective_valid(objective, maxPenalty)
      list(
        objective = if (raw.valid) objective else penalty,
        num.feval = 1L,
        num.feval.fast = if (.npscoefbw_fast_eligible(bws, eval.zdat = ctx$zdat.df)) 1L else 0L,
        raw.valid = raw.valid
      )
    }
  }, np_nn_candidate_invalid = function(e) {
    list(objective = penalty, num.feval = 1L, num.feval.fast = 0L,
         raw.valid = FALSE)
  })
}

.npscoefbw_nomad_search <- function(xdat,
                                    ydat,
                                    zdat,
                                    bws,
                                    reg.args,
                                    opt.args,
                                    degree.search,
                                    nomad.inner.nmulti = 0L,
                                    nomad.opts = list(),
                                    source = "explicit",
                                    reason = NULL,
                                    progress_label = NULL) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")
  if (is.null(opt.args$nomad.opts) && length(nomad.opts))
    opt.args$nomad.opts <- nomad.opts

  if (!identical(opt.args$bandwidth.compute, TRUE))
    stop("automatic degree search with search.engine='nomad' requires bandwidth.compute=TRUE")

  eval.zdat <- if (is.null(zdat)) xdat else zdat
  automatic.start <- is.numeric(bws) && length(bws) > 0L && all(bws == 0)

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template.type <- if (is.recursive(bws) && !is.null(bws$type) && length(bws$type)) {
    as.character(bws$type[1L])
  } else if (!is.null(template.reg.args$bwtype) && length(template.reg.args$bwtype)) {
    as.character(template.reg.args$bwtype[1L])
  } else {
    "fixed"
  }
  template.bws <- bws
  if (!identical(template.type, "fixed") &&
      is.numeric(template.bws) &&
      length(template.bws) > 0L &&
      all(template.bws == 0)) {
    nn.start <- max(2L, min(as.integer(NROW(eval.zdat)) - 1L, as.integer(round(sqrt(NROW(eval.zdat))))))
    template.bws <- rep.int(as.double(nn.start), length(template.bws))
  }

  template <- .npscoefbw_build_scbandwidth(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = template.bws,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args
  )

  if (!(template$type %in% c("fixed", "generalized_nn", "adaptive_nn")))
    stop("automatic degree search with search.engine='nomad' requires bwtype='fixed', 'generalized_nn', or 'adaptive_nn'")
  nn.search <- template$type %in% c("generalized_nn", "adaptive_nn")
  setup <- .npregbw_nomad_bw_setup(xdat = eval.zdat, template = template, allow.extended.nn = TRUE)
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)
  ndeg <- length(degree.search$start.degree)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(NCOL(eval.zdat)) else npValidateNmulti(opt.args$nmulti[1L])

  bw_bounds <- .npregbw_nomad_bw_bounds(template = template, setup = setup)
  opt.value <- function(name, default) {
    if (is.null(opt.args[[name]])) default else opt.args[[name]]
  }
  bw_start_bounds <- .np_nomad_bw_restart_start_bounds(
    bounds = bw_bounds,
    setup = setup,
    opt.value = opt.value,
    where = "npscoefbw"
  )

  x0 <- c(
    .npregbw_nomad_complete_bw_start_point(
      point = if (all(template$bw == 0)) NULL else .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup),
      bounds = bw_bounds,
      setup = setup,
      initial = bw_start_bounds$initial,
      where = "npscoefbw"
    ),
    as.integer(degree.search$start.degree)
  )
  lb <- c(bw_bounds$lower, degree.search$lower)
  ub <- c(bw_bounds$upper, degree.search$upper)
  bbin <- c(bw_bounds$bbin, rep.int(1L, ndeg))
  coordinate.roles <- .np_nomad_coordinate_roles(bw_bounds, degree.search)
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0
  ctx <- .npscoefbw_nomad_context_prepare(xdat = xdat, ydat = ydat, zdat = zdat)
  raw.valid.points <- new.env(hash = TRUE, parent = emptyenv())
  invalid.objective.override <- NULL
  evaluator.error <- NULL
  point.key <- function(point) {
    paste(sprintf("%.17g", as.numeric(point)), collapse = ",")
  }
  point.degree <- function(point) {
    .np_degree_clip_to_grid(
      as.integer(round(point[ncon + ncat + seq_len(ndeg)])),
      degree.search$candidates
    )
  }

  .np_nomad_baseline_note(degree.search$start.degree)

  evaluate.point <- function(point) {
    point <- as.numeric(point)
    degree <- point.degree(point)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)
    if (nn.search && !.npscoef_candidate_is_admissible(
          param = bw_vec,
          bwtype = template$type,
          nobs = NROW(eval.zdat),
          icon = template$icon
        )) {
      .np_nn_abort_candidate_invalid(
        "npscoefbw NOMAD candidate lies outside the admissible NN domain",
        owner = "npscoefbw NOMAD degree objective",
        point = point
      )
    }

    eval.reg.args <- reg.args
    eval.reg.args$regtype <- "lp"
    eval.reg.args$degree <- degree
    eval.reg.args$bernstein.basis <- degree.search$bernstein.basis
    eval.reg.args$bwtype <- template$type

    tbw <- .npscoefbw_build_scbandwidth(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = eval.reg.args
    )
    tbw <- .npscoefbw_normalize_nomad_scbw(
      scbw = tbw,
      eval.zdat = eval.zdat,
      bw = bw_vec
    )

    out <- .npscoefbw_nomad_eval_direct(
      ctx = ctx,
      bws = tbw,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier,
      invalid.objective = invalid.objective.override
    )
    if (nn.search)
      assign(point.key(point), isTRUE(out$raw.valid), envir = raw.valid.points)
    nomad.num.feval.total <<- nomad.num.feval.total + as.numeric(out$num.feval[1L])
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])

    list(
      objective = out$objective,
      degree = degree,
      num.feval = out$num.feval,
      raw.valid = out$raw.valid
    )
  }

  eval_fun <- function(point) {
    tryCatch(
      evaluate.point(point),
      error = function(e) {
        if (!inherits(e, "np_nn_candidate_invalid"))
          evaluator.error <<- e
        stop(e)
      }
    )
  }

  powell.payload.raw.valid <- FALSE
  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    degree <- as.integer(best_record$degree)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)
    powell.elapsed <- NA_real_

    build_direct_payload <- function() {
      final.reg.args <- reg.args
      final.reg.args$regtype <- "lp"
      final.reg.args$degree <- degree
      final.reg.args$bernstein.basis <- degree.search$bernstein.basis
      final.reg.args$bwtype <- template$type

      tbw <- .npscoefbw_build_scbandwidth(
        xdat = xdat,
        ydat = ydat,
        zdat = zdat,
        bws = bw_vec,
        bandwidth.compute = FALSE,
        reg.args = final.reg.args
      )
      tbw$fval <- as.numeric(best_record$objective)
      tbw$ifval <- as.numeric(best_record$objective)
      tbw$num.feval <- as.numeric(nomad.num.feval.total)
      tbw$num.feval.fast <- as.numeric(nomad.num.feval.fast.total)
      tbw$numimp <- 0
      tbw$fval.vector <- as.numeric(best_record$objective)
      tbw$total.time <- NA_real_
      if (!is.null(tbw$method) && length(tbw$method))
        tbw$pmethod <- bwmToPrint(as.character(tbw$method[1L]))
      tbw
    }

    direct.payload <- build_direct_payload()
    direct.objective <- as.numeric(best_record$objective)

    point.is.raw.valid <- !nn.search || (
      exists(point.key(point), envir = raw.valid.points, inherits = FALSE) &&
      isTRUE(get(point.key(point), envir = raw.valid.points, inherits = FALSE))
    )
    if (identical(degree.search$engine, "nomad+powell") &&
        point.is.raw.valid) {
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.reg.args$bwtype <- template$type
      hot.opt.args <- .np_nomad_powell_hotstart_opt_args(
        opt.args,
        strategy = "single_iteration",
        remin = isTRUE(opt.args$powell.remin)
      )
      hot.opt.args$bwsolver <- NULL
      powell.start <- proc.time()[3L]
      hot.payload <- .np_nomad_with_powell_progress(
        degree = degree,
        best_record = best_record,
        expr = local({
          .npscoefbw_run_fixed_degree(
            xdat = xdat,
            ydat = ydat,
            zdat = zdat,
            bws = bw_vec,
            reg.args = hot.reg.args,
            opt.args = hot.opt.args
          )
        })
      )
      powell.elapsed <- proc.time()[3L] - powell.start
      direct.payload$num.feval <- as.numeric(direct.payload$num.feval[1L]) + as.numeric(hot.payload$num.feval[1L])
      direct.payload$num.feval.fast <- as.numeric(direct.payload$num.feval.fast[1L]) + as.numeric(hot.payload$num.feval.fast[1L])
      hot.payload$num.feval <- direct.payload$num.feval
      hot.payload$num.feval.fast <- direct.payload$num.feval.fast
      if (!is.null(hot.payload$method) && length(hot.payload$method))
        hot.payload$pmethod <- bwmToPrint(as.character(hot.payload$method[1L]))
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min")) {
        powell.payload.raw.valid <<- TRUE
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  run.search <- function(start, nmulti, preserve.start.degree = FALSE) {
    start.degree <- if (isTRUE(preserve.start.degree)) {
      as.integer(start[ncon + ncat + seq_len(ndeg)])
    } else {
      degree.search$start.degree
    }
    result <- tryCatch(.np_nomad_search(
      engine = degree.search$engine,
      baseline_record = baseline.record,
      start_degree = start.degree,
      x0 = start,
      bbin = bbin,
      lb = lb,
      ub = ub,
      eval_fun = eval_fun,
      build_payload = build_payload,
      direction = "min",
      objective_name = "fval",
      nmulti = nmulti,
      nomad.inner.nmulti = nomad.inner.nmulti,
      random.seed = if (!is.null(opt.args$random.seed)) opt.args$random.seed else 42L,
      remin = isTRUE(opt.args$nomad.remin),
      nomad.opts = if (is.null(opt.args$nomad.opts)) list() else opt.args$nomad.opts,
      native.r.bridge = TRUE,
      source = source,
      reason = reason,
      progress_label = progress_label,
      handoff_before_build = nn.search,
      start.lower = c(bw_start_bounds$lower, degree.search$lower),
      start.upper = c(bw_start_bounds$upper, degree.search$upper),
      coordinate.roles = coordinate.roles,
      degree_spec = list(
        initial = start.degree,
        lower = degree.search$lower,
        upper = degree.search$upper,
        basis = degree.search$basis,
        nobs = degree.search$nobs,
        user_supplied = isTRUE(preserve.start.degree) || degree.search$start.user
      )
    ), error = function(e) {
      if (!is.null(evaluator.error))
        stop(evaluator.error)
      stop(e)
    })
    if (!is.null(evaluator.error))
      stop(evaluator.error)
    result
  }

  result.raw.valid <- function(result) {
    if (!is.null(result$best_point)) {
      key <- point.key(result$best_point)
      if (exists(key, envir = raw.valid.points, inherits = FALSE) &&
          isTRUE(get(key, envir = raw.valid.points, inherits = FALSE))) {
        return(TRUE)
      }
    }
    identical(degree.search$engine, "nomad+powell") &&
      isTRUE(powell.payload.raw.valid)
  }

  append.search.history <- function(incumbent, recovered) {
    offset <- length(incumbent$restart.results)
    if (length(recovered$restart.results)) {
      for (i in seq_along(recovered$restart.results))
        recovered$restart.results[[i]]$restart <- offset + i
    }
    recovered$baseline <- incumbent$baseline
    recovered$restart.start.info <- incumbent$restart.start.info
    if (!is.null(recovered$best.restart) && length(recovered$best.restart) &&
        !is.na(recovered$best.restart)) {
      recovered$best.restart <- offset + recovered$best.restart
    }
    if (!is.null(recovered$nomad.remin.index) &&
        length(recovered$nomad.remin.index) &&
        !is.na(recovered$nomad.remin.index)) {
      recovered$nomad.remin.index <- offset + recovered$nomad.remin.index
    }
    recovered$restart.starts <- c(
      incumbent$restart.starts,
      recovered$restart.starts
    )
    recovered$restart.degree.starts <- c(
      incumbent$restart.degree.starts,
      recovered$restart.degree.starts
    )
    recovered$restart.bandwidth.starts <- c(
      incumbent$restart.bandwidth.starts,
      recovered$restart.bandwidth.starts
    )
    recovered$restart.results <- c(
      incumbent$restart.results,
      recovered$restart.results
    )
    recovered$n.unique <- incumbent$n.unique + recovered$n.unique
    recovered$n.visits <- incumbent$n.visits + recovered$n.visits
    recovered$n.cached <- incumbent$n.cached + recovered$n.cached
    recovered$nomad.time <- incumbent$nomad.time + recovered$nomad.time
    recovered$optim.time <- incumbent$optim.time + recovered$optim.time
    trace <- rbind(incumbent$trace, recovered$trace)
    if (nrow(trace)) {
      trace$trace_id <- seq_len(nrow(trace))
      trace$eval_id <- seq_len(nrow(trace))
    }
    recovered$trace <- trace
    recovered
  }

  search.result <- run.search(x0, nomad.nmulti)
  if (!nn.search)
    return(search.result)

  if (!result.raw.valid(search.result) && automatic.start) {
    recovery.start <- as.numeric(search.result$restart.starts[[1L]])
    recovery.raw.eval <- function(point) {
      out <- eval_fun(point)
      if (!isTRUE(out$raw.valid)) {
        .np_nn_abort_candidate_invalid(
          "npscoefbw NOMAD degree objective returned an invalid raw objective",
          owner = "npscoefbw NOMAD degree objective",
          point = point,
          raw.objective = out$objective
        )
      }
      as.double(out$objective)
    }
    ordinary.caps <- rep.int(NROW(eval.zdat) - 1L, ncon)
    ordinary.start <- recovery.start[seq_len(ncon)]
    recovery <- if (length(ordinary.start) &&
                    all(ordinary.caps >= 2L) &&
                    all(ordinary.start >= 2L) &&
                    all(ordinary.start <= ordinary.caps)) {
      .np_nn_find_raw_valid_start(
        point = recovery.start,
        nn.indices = seq_len(ncon),
        caps = ordinary.caps,
        raw.eval = recovery.raw.eval
      )
    } else {
      list(found = FALSE, point = NULL, objective = NA_real_, evaluations = 0L)
    }
    if (isTRUE(recovery$found)) {
      incumbent <- search.result
      invalid.objective.override <- sqrt(.Machine$double.xmax)
      search.result <- append.search.history(
        incumbent = incumbent,
        recovered = run.search(recovery$point, 1L, preserve.start.degree = TRUE)
      )
    }
  }

  if (!result.raw.valid(search.result))
    stop("npscoefbw NOMAD degree search did not return a raw-valid solution",
         call. = FALSE)

  search.result$best_payload <- .npscoefbw_normalize_nomad_scbw(
    scbw = search.result$best_payload,
    eval.zdat = eval.zdat,
    bw = search.result$best_payload$bw
  )
  if (!isTRUE(powell.payload.raw.valid)) {
    final.raw <- .npscoefbw_nomad_eval_direct(
      ctx = ctx,
      bws = search.result$best_payload,
      invalid.penalty = "large",
      invalid.objective = sqrt(.Machine$double.xmax)
    )
    reported.objective <- as.numeric(search.result$best_payload$fval[1L])
    if (!isTRUE(final.raw$raw.valid) ||
        !identical(as.numeric(final.raw$objective), reported.objective)) {
      stop("npscoefbw NOMAD degree search failed final raw-objective certification",
           call. = FALSE)
    }
  }
  search.result
}

.npscoefbw_degree_search_controls <- function(regtype,
                                              regtype.named,
                                              cv.iterate,
                                              cv.iterate.named,
                                              bandwidth.compute,
                                              ncon,
                                              nobs,
                                              basis,
                                              degree.select,
                                              search.engine,
                                              degree.min,
                                              degree.max,
                                              degree.start,
                                              degree.restarts,
                                              degree.max.cycles,
                                              degree.verify,
                                              bernstein.basis,
                                              bernstein.named,
                                              nomad.source = "explicit",
                                              nomad.auto.filled = character()) {
  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  if (identical(degree.select, "manual"))
    return(NULL)
  resolved <- .np_degree_resolve_auto_engine(
    search.engine = search.engine,
    degree.select = degree.select,
    ncon = ncon,
    source = nomad.source,
    auto.filled = nomad.auto.filled
  )
  search.engine <- .npscoefbw_nomad_controls(resolved$search.engine)
  degree.select <- resolved$degree.select

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (!isTRUE(bandwidth.compute))
    stop("automatic degree search requires bandwidth.compute=TRUE")
  if (isTRUE(cv.iterate.named) && isTRUE(cv.iterate))
    stop("automatic degree search currently requires cv.iterate=FALSE")
  if (ncon < 1L)
    stop("automatic degree search requires at least one continuous smoothing predictor")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = ncon,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  baseline.degree <- rep.int(0L, ncon)
  default.start.degree <- if (identical(search.engine, "cell")) {
    baseline.degree
  } else {
    rep.int(1L, ncon)
  }
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, default.start.degree))
  } else {
    start.raw <- npValidateGlpDegree(regtype = "lp", degree = degree.start, ncon = ncon, argname = "degree.start")
    out.of.range <- vapply(seq_len(ncon), function(j) !(start.raw[j] %in% bounds$candidates[[j]]), logical(1))
    if (any(out.of.range))
      stop("degree.start must lie within the searched degree candidates for every continuous smoothing predictor")
    start.raw
  }

  list(
    method = if (identical(search.engine, "cell")) degree.select else search.engine,
    engine = search.engine,
    candidates = bounds$candidates,
    lower = bounds$lower,
    upper = bounds$upper,
    grid.size = bounds$grid.size,
    singleton = bounds$singleton,
    fixed.degree = bounds$fixed.degree,
    baseline.degree = baseline.degree,
    start.degree = start.degree,
    start.user = !is.null(degree.start),
    basis = if (missing(basis) || is.null(basis)) "glp" else as.character(basis[1L]),
    nobs = as.integer(nobs[1L]),
    restarts = npValidateNonNegativeInteger(degree.restarts, "degree.restarts"),
    max.cycles = npValidatePositiveInteger(degree.max.cycles, "degree.max.cycles"),
    verify = npValidateScalarLogical(degree.verify, "degree.verify"),
    bernstein.basis = bern.auto,
    source = resolved$source,
    reason = resolved$reason
  )
}

.npscoefbw_attach_degree_search <- function(bws, search_result) {
  metadata <- .np_degree_search_metadata(search_result, default_direction = "min")

  if (isTRUE(search_result$native) &&
      isTRUE(getOption("np.developer.native.nomad.diagnostics", FALSE)) &&
      !is.null(search_result$native.diagnostics)) {
    attr(bws, "native.nomad.diagnostics") <- search_result$native.diagnostics
  }

  if (!is.null(search_result$nomad.time))
    bws$nomad.time <- as.numeric(search_result$nomad.time[1L])
  if (!is.null(search_result$powell.time))
    bws$powell.time <- as.numeric(search_result$powell.time[1L])
  if (!is.null(search_result$optim.time) && is.finite(search_result$optim.time))
    bws$total.time <- as.numeric(search_result$optim.time[1L])
  bws <- .np_attach_nomad_restart_summary(bws, search_result)
  bws$degree.search <- metadata
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
           degree.select = c("manual", "coordinate", "exhaustive"),
           search.engine = c("nomad+powell", "cell", "nomad"),
           nomad = FALSE,
           nomad.nmulti = 0L,
           degree.min = NULL,
           degree.max = NULL,
           degree.start = NULL,
           degree.restarts = 0L,
           degree.max.cycles = 20L,
           degree.verify = FALSE,
           nmulti,
           nomad.remin = FALSE,
           powell.remin = TRUE,
           okertype,
           optim.abstol,
           optim.maxattempts,
           optim.maxit,
           optim.method,
           optim.reltol,
           random.seed,
           regtype,
           ukertype,
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           lbd.init = 0.5,
           hbd.init = 1.5,
           dfac.init = 1.0,
           scale.factor.search.lower = NULL,
           ...,
           nomad.opts = list()){

    nomad.opts <- .np_nomad_normalize_user_opts(nomad.opts, "npscoefbw")
    dots <- list(...)
    if (length(nomad.opts))
      dots$nomad.opts <- nomad.opts
    npRejectUnsupportedBwsolver(dots, "npscoefbw")

    if (!missing(bwmethod) && identical(match.arg(bwmethod, c("cv.ls", "manual")), "manual") &&
        missing(bws))
      stop("bwmethod='manual' requires argument 'bws'")

    miss.z <- missing(zdat)
    xdat <- toFrame(xdat)
    
    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector or a factor")

    if(!miss.z)
      zdat <- toFrame(zdat)

    mc <- match.call(expand.dots = FALSE)
    mc.names <- names(mc)
    nomad.shortcut <- .np_prepare_nomad_shortcut(
      nomad = nomad,
      call_names = mc.names,
      preset = list(
        regtype = "lp",
        search.engine = "nomad+powell",
        degree.select = "coordinate",
        bernstein.basis = TRUE,
        degree.min = 0L,
        degree.max = 10L,
        degree.verify = FALSE,
        bwtype = "fixed"
      ),
      values = list(
        regtype = if ("regtype" %in% mc.names) regtype else NULL,
        search.engine = if ("search.engine" %in% mc.names) search.engine else NULL,
        degree.select = if ("degree.select" %in% mc.names) degree.select else NULL,
        bernstein.basis = if ("bernstein.basis" %in% mc.names) bernstein.basis else NULL,
        degree.min = if ("degree.min" %in% mc.names) degree.min else NULL,
        degree.max = if ("degree.max" %in% mc.names) degree.max else NULL,
        degree.verify = if ("degree.verify" %in% mc.names) degree.verify else NULL,
        bwtype = if ("bwtype" %in% mc.names) bwtype else NULL,
        degree = if ("degree" %in% mc.names) degree else NULL
      ),
      where = "npscoefbw"
    )

    if (isTRUE(nomad.shortcut$enabled)) {
      if ("degree" %in% mc.names)
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$bwtype, c("fixed", "generalized_nn", "adaptive_nn")))[1L] %in%
              c("fixed", "generalized_nn", "adaptive_nn")))
        stop("nomad=TRUE requires bwtype='fixed', 'generalized_nn', or 'adaptive_nn'")
      if ("degree.select" %in% mc.names &&
          identical(as.character(match.arg(nomad.shortcut$values$degree.select, c("manual", "coordinate", "exhaustive")))[1L], "manual"))
        stop("nomad=TRUE requires automatic degree search; use degree.select='coordinate' or 'exhaustive'")
      if (!identical(nomad.shortcut$metadata$source, "auto") &&
          "search.engine" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$search.engine, c("nomad+powell", "cell", "nomad")))[1L] %in%
              c("nomad", "nomad+powell")))
        stop("nomad=TRUE requires search.engine='nomad' or 'nomad+powell'")
      if ("degree.verify" %in% mc.names &&
          isTRUE(npValidateScalarLogical(nomad.shortcut$values$degree.verify, "degree.verify")))
        stop("nomad=TRUE currently requires degree.verify=FALSE")
    }

    regtype.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "regtype")
    bernstein.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "bernstein.basis")
    cv.iterate.named <- any(mc.names == "cv.iterate")
    regtype.value <- if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else "lc"
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    degree.search <- .npscoefbw_degree_search_controls(
      regtype = regtype.value,
      regtype.named = regtype.named,
      cv.iterate = cv.iterate,
      cv.iterate.named = cv.iterate.named,
      bandwidth.compute = bandwidth.compute,
      ncon = sum(if (miss.z) untangle(xdat)$icon else untangle(zdat)$icon),
      nobs = NROW(xdat),
      basis = if ("basis" %in% mc.names) basis else "glp",
      degree.select = degree.select.value,
      search.engine = if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell",
      degree.min = nomad.shortcut$values$degree.min,
      degree.max = nomad.shortcut$values$degree.max,
      degree.start = if ("degree.start" %in% mc.names) degree.start else NULL,
      degree.restarts = if ("degree.restarts" %in% mc.names) degree.restarts else 0L,
      degree.max.cycles = if ("degree.max.cycles" %in% mc.names) degree.max.cycles else 20L,
      degree.verify = if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE,
      bernstein.basis = if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else bernstein.basis,
      bernstein.named = bernstein.named,
      nomad.source = nomad.shortcut$metadata$source,
      nomad.auto.filled = nomad.shortcut$metadata$auto.filled
    )
    nomad.inner.named <- "nomad.nmulti" %in% mc.names
    nomad.inner.nmulti <- if (nomad.inner.named) {
      npValidateNonNegativeInteger(nomad.nmulti, "nomad.nmulti")
    } else {
      0L
    }
    if (nomad.inner.named &&
        (is.null(degree.search) || !(degree.search$engine %in% c("nomad", "nomad+powell")))) {
      stop("nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad+powell'")
    }
    degree.setup <- npSetupGlpDegree(
      regtype = regtype.value,
      degree = if ("degree" %in% mc.names) degree else NULL,
      ncon = sum(if (miss.z) untangle(xdat)$icon else untangle(zdat)$icon),
      degree.select = degree.select.value
    )
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)

    ## first grab dummy args for scbandwidth() and perform 'bootstrap'
    ## bandwidth call

    margs <- c("regtype", "basis", "degree", "bernstein.basis",
               "bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype",
               "scale.factor.search.lower")

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
    if (!("degree" %in% names(sbw.args)) && !is.null(degree.setup))
      sbw.args$degree <- degree.setup
    reg.args <- sbw.args[setdiff(names(sbw.args), c("bw", "nobs", "xdati", "ydati", "zdati", "xnames", "ynames", "znames", "bandwidth.compute"))]
    if (!is.null(degree.search))
      reg.args$bernstein.basis <- degree.search$bernstein.basis
    tbw <- do.call(scbandwidth, sbw.args)
    tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)

    ## next grab dummies for actual bandwidth selection and perform call
    margs <- c("zdat",
               "nmulti",
               "powell.remin",
               "random.seed",
               "nomad.opts",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.factor.search.lower",
               "cv.iterate",
               "cv.num.iterations",
               "backfit.iterate",
               "backfit.maxiter",
               "backfit.tol",
               "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (any.m) {
      nms <- mc.names[m]
      opt.args <- mget(nms, envir = environment(), inherits = FALSE)
    } else {
      opt.args <- list()
    }
    opt.args <- c(list(bandwidth.compute = bandwidth.compute), opt.args)
    if ("nomad.opts" %in% names(dots))
      opt.args$nomad.opts <- dots$nomad.opts
    reg.args$scale.factor.search.lower <- scale.factor.search.lower
    opt.args$scale.factor.search.lower <- scale.factor.search.lower

    if (!is.null(degree.search)) {
      eval_fun <- function(degree.vec) {
        cell.reg.args <- reg.args
        cell.reg.args$regtype <- "lp"
        cell.reg.args$degree <- as.integer(degree.vec)
        cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
        cell.bws <- .npscoefbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          zdat = if (miss.z) NULL else zdat,
          bws = bws,
          reg.args = cell.reg.args,
          opt.args = opt.args
        )
        list(
          objective = as.numeric(cell.bws$fval[1L]),
          payload = cell.bws,
          num.feval = if (!is.null(cell.bws$num.feval)) as.numeric(cell.bws$num.feval[1L]) else NA_real_,
          nn.cache = cell.bws$nn.cache
        )
      }

      if (isTRUE(degree.search$singleton)) {
        search.result <- .np_degree_singleton_search_result(
          degree.search = degree.search,
          eval_result = eval_fun(degree.search$fixed.degree),
          direction = "min",
          objective_name = "fval"
        )
      } else if (identical(degree.search$engine, "cell")) {
        search.result <- .np_degree_search(
          method = degree.search$method,
          candidates = degree.search$candidates,
          baseline_degree = degree.search$baseline.degree,
          start_degree = degree.search$start.degree,
          restarts = degree.search$restarts,
          max_cycles = degree.search$max.cycles,
          verify = degree.search$verify,
          eval_fun = eval_fun,
          direction = "min",
          trace_level = "full",
          source = degree.search$source,
          reason = degree.search$reason,
          objective_name = "fval"
        )
      } else {
        search.result <- .npscoefbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          zdat = if (miss.z) NULL else zdat,
          bws = bws,
          reg.args = reg.args,
          opt.args = opt.args,
          degree.search = degree.search,
          nomad.inner.nmulti = nomad.inner.nmulti,
          nomad.opts = if (is.null(opt.args$nomad.opts)) list() else opt.args$nomad.opts,
          source = degree.search$source,
          reason = degree.search$reason,
          progress_label = .np_degree_search_label(degree.search$engine, degree.search$source)
        )
      }
      tbw <- .npscoefbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
    } else {
      scbw.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
      if (!miss.z)
        scbw.args$zdat <- zdat
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting smooth coefficient bandwidth",
        do.call(npscoefbw.scbandwidth, scbw.args)
      )
    }

    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)
    
  }
