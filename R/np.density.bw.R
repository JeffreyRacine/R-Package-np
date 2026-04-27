npudensbw <- function(...){
  mc <- match.call(expand.dots = FALSE)
  npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npudensbw")
  target <- .np_bw_dispatch_target(dots = mc$...,
                                   data_arg_names = "dat",
                                   eval_env = parent.frame())
  UseMethod("npudensbw", target)
}

npudensbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    formula.terms <- terms(formula)
    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = formula.terms,
                        data = environment(formula),
                        eval_env = environment(formula))
    else .np_terms_ts_mask(terms_obj = formula.terms,
                           data = data,
                           eval_env = environment(formula))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    if(all(orig.ts)){
      args <- (as.list(attr(formula.terms, "variables"))[-1])
      formula <- formula.terms
      attr(formula, "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
      mf[["formula"]] <- formula
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(formula.terms, "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      formula <- formula.terms
      attr(formula, "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
      mf[["formula"]] <- formula
    }

    mf[[1]] <- as.name("model.frame")
    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())

    if (attr(attr(mf, "terms"), "response") != 0)
      stop("invalid density formula")

    dat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    tbw <- npudensbw(dat = dat, ...)
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$terms <- attr(mf,"terms")
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw
  }


npudensbw.NULL <-
  function(dat = stop("invoked without input data 'dat'"),
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npudensbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    t.names <- NULL
    if(!is.data.frame(dat) && !is.matrix(dat))
      t.names <- paste(deparse(substitute(dat)), collapse = "")

    dat = toFrame(dat)

    if(!is.null(t.names))
      names(dat) <- t.names

    bws = double(dim(dat)[2])

    tbw <- npudensbw.default(dat = dat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw
  }

.npudensbw_assert_bounded_cvls_supported <- function(bws,
                                                     where = "npudensbw()") {
  method <- if (!is.null(bws$method) && length(bws$method)) {
    as.character(bws$method[1L])
  } else {
    "cv.ml"
  }

  if (!identical(method, "cv.ls"))
    return(invisible(TRUE))

  ckerlb <- if (is.null(bws$ckerlb)) numeric(0L) else bws$ckerlb[bws$icon]
  ckerub <- if (is.null(bws$ckerub)) numeric(0L) else bws$ckerub[bws$icon]
  bounded.x <- length(ckerlb) > 0L && any(is.finite(ckerlb) | is.finite(ckerub))

  if (!bounded.x)
    return(invisible(TRUE))

  if (bws$ncon < 1L || bws$ncon > 2L) {
    stop(
      sprintf(
        "%s bounded npudens cv.ls currently supports up to two continuous variables with optional ordered/unordered discrete components",
        where
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

npudensbw.bandwidth <-
  function(dat = stop("invoked without input data 'dat'"),
           bws,
           bandwidth.compute = TRUE,
           cfac.dir = 2.5*(3.0-sqrt(5)),
           scale.factor.init = 0.5,
           dfac.dir = 0.25*(3.0-sqrt(5)),
           dfac.init = 0.375,
           dfc.dir = 3,
           ftol = 1.490116e-07,
           scale.factor.init.upper = 2.0,
           hbd.dir = 1,
           hbd.init = 0.9,
           initc.dir = 1.0,
           initd.dir = 1.0,
           invalid.penalty = c("baseline","dbmax"),
           itmax = 10000,
           lbc.dir = 0.5,
           scale.factor.init.lower = 0.1,
           lbd.dir = 0.1,
           lbd.init = 0.1,
           nmulti,
           penalty.multiplier = 10,
           remin = TRUE,
           scale.init.categorical.sample = FALSE,
           scale.factor.search.lower = NULL,
           small = 1.490116e-05,
           tol = 1.490116e-04,
           transform.bounds = FALSE,
           ...){
    elapsed.start <- proc.time()[3]
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    remin <- npValidateScalarLogical(remin, "remin")
    scale.init.categorical.sample <-
      npValidateScalarLogical(scale.init.categorical.sample, "scale.init.categorical.sample")
    transform.bounds <- npValidateScalarLogical(transform.bounds, "transform.bounds")
    itmax <- npValidatePositiveInteger(itmax, "itmax")
    ftol <- npValidatePositiveFiniteNumeric(ftol, "ftol")
    tol <- npValidatePositiveFiniteNumeric(tol, "tol")
    small <- npValidatePositiveFiniteNumeric(small, "small")
    penalty.multiplier <- npValidatePositiveFiniteNumeric(penalty.multiplier, "penalty.multiplier")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    if (!missing(nmulti))
      nmulti <- npValidateNmulti(nmulti)
    .npRmpi_require_active_slave_pool(where = "npudensbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    dat = toFrame(dat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(dim(dat)[2])
    }
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)

    if (length(bws$bw) != dim(dat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of",
           "'dat'"))

    if ((any(bws$icon) &&
         !all(vapply(as.data.frame(dat[, bws$icon]), inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iord) &&
         !all(vapply(as.data.frame(dat[, bws$iord]), inherits, logical(1), "ordered"))) ||
        (any(bws$iuno) &&
         !all(vapply(as.data.frame(dat[, bws$iuno]), inherits, logical(1), "factor"))))
      stop(paste("supplied bandwidths do not match", "'dat'", "in type"))

    dat <- na.omit(dat)
    rows.omit <- unclass(na.action(dat))

    nrow = dim(dat)[1]
    ncol = dim(dat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    dat = toMatrix(dat)

    duno = dat[, bws$iuno, drop = FALSE]
    dcon = dat[, bws$icon, drop = FALSE]
    dord = dat[, bws$iord, drop = FALSE]

    tbw <- bws

    mysd <- EssDee(dcon)
    nconfac <- nrow^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- nrow^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

    if (bandwidth.compute){
      cont.start <- npContinuousSearchStartControls(
        scale.factor.init.lower,
        scale.factor.init.upper,
        scale.factor.init,
        scale.factor.search.lower,
        where = "npudensbw"
      )
      myopti = list(num_obs_train = dim(dat)[1],
        iMultistart = IMULTI_TRUE,
        iNum_Multistart = nmulti,
        int_use_starting_values = (if (all(bws$bw==0)) USE_START_NO else USE_START_YES),
        int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
        BANDWIDTH_den_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=(if (remin) RE_MIN_TRUE else RE_MIN_FALSE),
        int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
        bwmethod = switch(bws$method,
          cv.ml = BWM_CVML,
          cv.ls = BWM_CVLS),
        ckerneval = switch(bws$ckertype,
          gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
          uniform = CKER_UNI,
          "truncated gaussian" = CKER_TGAUSS),
        ukerneval = switch(bws$ukertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        okerneval = switch(bws$okertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_NLR,
        "racineliyan" = OKER_RLY),
        nuno = dim(duno)[2],
        nord = dim(dord)[2],
        ncon = dim(dcon)[2],
        old.dens = FALSE,
        int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
        scale.init.categorical.sample = scale.init.categorical.sample,
        dfc.dir = dfc.dir,
        transform.bounds = transform.bounds)


      myoptd = list(ftol=ftol, tol=tol, small=small,
        lbc.dir = lbc.dir, cfac.dir = cfac.dir, initc.dir = initc.dir,
        lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir,
        lbc.init = cont.start$scale.factor.init.lower,
        hbc.init = cont.start$scale.factor.init.upper,
        cfac.init = cont.start$scale.factor.init,
        lbd.init = lbd.init, hbd.init = hbd.init, dfac.init = dfac.init,
        nconfac = nconfac, ncatfac = ncatfac, memfac = 0,
        scale.factor.lower.bound = scale.factor.search.lower)
      cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

      .npudensbw_assert_bounded_cvls_supported(tbw, where = "npudensbw()")

      if (bws$method != "normal-reference"){
        myout <-
          .Call("C_np_density_bw",
                as.double(duno), as.double(dord), as.double(dcon),
                as.double(mysd),
                as.integer(myopti), as.double(myoptd),
                as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
                as.integer(nmulti),
                as.integer(penalty_mode),
                as.double(penalty.multiplier),
                as.double(cker.bounds.c$lb),
                as.double(cker.bounds.c$ub),
                PACKAGE="npRmpi")
        total.time <- proc.time()[3] - elapsed.start
      } else {
        nbw = double(ncol)
        if (bws$ncon > 0){
          con_idx <- seq_len(bws$ncon)
          nbw[con_idx] = 1.059224
          if(!bws$scaling)
            nbw[con_idx] = nbw[con_idx] * mysd * nconfac
        }
        myout= list( bw = nbw, fval = c(NA,NA) )
        total.time <- NA
      }

      rorder = numeric(ncol)
      ord_idx <- seq_len(ncol)
      rorder[c(ord_idx[bws$icon], ord_idx[bws$iuno], ord_idx[bws$iord])] <- ord_idx

      tbw$bw <- myout$bw[rorder]

      tbw$fval = myout$fval[1]
      tbw$ifval = myout$fval[2]
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
      tbw$fval.history <- myout$fval.history
      tbw$eval.history <- myout$eval.history
      tbw$invalid.history <- myout$invalid.history
      tbw$timing <- myout$timing
      tbw$total.time <- total.time
    }

    tbw$sfactor <- tbw$bandwidth <- tbw$bw

    if (tbw$nuno > 0){
      if(tbw$scaling){
        tbw$bandwidth[tbw$xdati$iuno] <- tbw$bandwidth[tbw$xdati$iuno]*ncatfac
      } else {
        tbw$sfactor[tbw$xdati$iuno] <- tbw$sfactor[tbw$xdati$iuno]/ncatfac
      }
    }

    if (tbw$nord > 0){
      if(tbw$scaling){
        tbw$bandwidth[tbw$xdati$iord] <- tbw$bandwidth[tbw$xdati$iord]*ncatfac
      } else {
        tbw$sfactor[tbw$xdati$iord] <- tbw$sfactor[tbw$xdati$iord]/ncatfac
      }
    }


    if (tbw$ncon > 0){
      dfactor <- mysd*nconfac

      if (tbw$scaling) {
        tbw$bandwidth[tbw$xdati$icon] <- tbw$bandwidth[tbw$xdati$icon]*dfactor
      } else {
        tbw$sfactor[tbw$xdati$icon] <- tbw$sfactor[tbw$xdati$icon]/dfactor
      }
    }

    tbw <- bandwidth(bw = tbw$bw,
                     bwmethod = tbw$method,
                     bwscaling = tbw$scaling,
                     bwtype = tbw$type,
                     ckertype = tbw$ckertype,
                     ckerorder = tbw$ckerorder,
                     ckerbound = tbw$ckerbound,
                     ckerlb = tbw$ckerlb,
                     ckerub = tbw$ckerub,
                     ukertype = tbw$ukertype,
                     okertype = tbw$okertype,
                     fval = tbw$fval,
                     ifval = tbw$ifval,
                     num.feval = tbw$num.feval,
                     num.feval.fast = tbw$num.feval.fast,
                     fval.history = tbw$fval.history,
                     eval.history = tbw$eval.history,
                     invalid.history = tbw$invalid.history,
                     nobs = tbw$nobs,
                     xdati = tbw$xdati,
                     xnames = tbw$xnames,
                     sfactor = tbw$sfactor,
                     bandwidth = tbw$bandwidth,
                     rows.omit = rows.omit,
                     nconfac = nconfac,
                     ncatfac = ncatfac,
                     sdev = mysd,
                     bandwidth.compute = bandwidth.compute,
                     timing = tbw$timing,
                     total.time = tbw$total.time)
    tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)

    tbw
  }

npudensbw.default <-
  function(dat = stop("invoked without input data 'dat'"),
           bws,
           bandwidth.compute = TRUE,
           bwmethod,
           bwscaling,
           bwtype,
           cfac.dir,
           scale.factor.init,
           ckerbound,
           ckerlb,
           ckerorder,
           ckertype,
           ckerub,
           dfac.dir,
           dfac.init,
           dfc.dir,
           ftol,
           scale.factor.init.upper,
           hbd.dir,
           hbd.init,
           initc.dir,
           initd.dir,
           invalid.penalty,
           itmax,
           lbc.dir,
           scale.factor.init.lower,
           lbd.dir,
           lbd.init,
           nmulti,
           okertype,
           penalty.multiplier,
           remin,
           scale.init.categorical.sample,
           scale.factor.search.lower = NULL,
           small,
           tol,
           transform.bounds,
           ukertype,
           ## dummy arguments for later passing into npudensbw.bandwidth
           ...){
    .npRmpi_require_active_slave_pool(where = "npudensbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    t.names <- NULL
    if(!is.data.frame(dat) && !is.matrix(dat))
      t.names <- paste(deparse(substitute(dat)), collapse = "")

    dat <- toFrame(dat)

    if(!is.null(t.names))
      names(dat) <- t.names

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    bw.args <- list(
      bw = bws,
      nobs = dim(dat)[1],
      xdati = untangle(dat),
      xnames = names(dat),
      bandwidth.compute = bandwidth.compute
    )
    if (!missing(bwmethod)) bw.args$bwmethod <- bwmethod
    if (!missing(bwscaling)) bw.args$bwscaling <- bwscaling
    if (!missing(bwtype)) bw.args$bwtype <- bwtype
    if (!missing(ckertype)) bw.args$ckertype <- ckertype
    if (!missing(ckerorder)) bw.args$ckerorder <- ckerorder
    if (!missing(ckerbound)) bw.args$ckerbound <- ckerbound
    if (!missing(ckerlb)) bw.args$ckerlb <- ckerlb
    if (!missing(ckerub)) bw.args$ckerub <- ckerub
    if (!missing(ukertype)) bw.args$ukertype <- ukertype
    if (!missing(okertype)) bw.args$okertype <- okertype
    tbw <- do.call(bandwidth, bw.args)


    ## next grab dummies for actual bandwidth selection and perform call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("bandwidth.compute", "nmulti", "remin", "itmax", "ftol", "tol",
               "small",
               "lbc.dir","dfc.dir","cfac.dir", "initc.dir",
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.init.categorical.sample",
               "scale.factor.search.lower",
               "invalid.penalty",
               "penalty.multiplier")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    bwsel.args <- list(dat = dat, bws = tbw)
    if (any.m) {
      nms <- mc.names[m]
      bwsel.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    tbw <- .np_progress_select_bandwidth_enhanced(
      "Selecting density bandwidth",
      do.call(npudensbw.bandwidth, bwsel.args)
    )

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }
