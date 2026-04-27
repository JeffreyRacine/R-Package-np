npudistbw <- function(...){
  mc <- match.call(expand.dots = FALSE)
  npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npudistbw")
  target <- .np_bw_dispatch_target(dots = mc$...,
                                   data_arg_names = c("dat", "gdat"),
                                   eval_env = parent.frame())
  UseMethod("npudistbw", target)
}

.npRmpi_npudistbw_bounded_adaptive_requested <- function(bwtype = NULL,
                                                         ckerbound = "none",
                                                         bandwidth.compute = TRUE) {
  bwtype.ch <- tolower(as.character(bwtype)[1L])
  ckerbound.ch <- tolower(as.character(ckerbound)[1L])
  isTRUE(bandwidth.compute) &&
    identical(bwtype.ch, "adaptive_nn") &&
    !identical(ckerbound.ch, "none")
}

npudistbw.formula <-
  function(formula, data, subset, na.action, call, gdata = NULL, ...){
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
      stop("invalid distribution formula")

    dat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    has.gval <- !is.null(gdata)
    if (has.gval) {
      gmf <- match.call(expand.dots = FALSE)
      gm <- match(c("formula", "gdata"), names(gmf), nomatch = 0)
      gmf <- gmf[c(1,gm)]

      gmf[[1]] <- as.name("model.frame")
      names(gmf)[3] <- "data"
      gmf.args <- as.list(gmf[-1L])
      gmf <- do.call(stats::model.frame, gmf.args, envir = parent.frame())

      gdat <- gmf[, attr(attr(gmf, "terms"),"term.labels"), drop = FALSE]

    }

    dots <- list(...)
    seed.args <- c(list(dat = dat), if (has.gval) list(gdat = gdat) else list(), dots)
    tbw <- do.call(npudistbw, seed.args)
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$terms <- attr(mf,"terms")
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw
  }


npudistbw.NULL <-
  function(dat = stop("invoked without input data 'dat'"),
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npudistbw()")
    dots <- list(...)
    if (.npRmpi_npudistbw_bounded_adaptive_requested(
      bwtype = dots$bwtype,
      ckerbound = if (is.null(dots$ckerbound)) "none" else dots$ckerbound,
      bandwidth.compute = if (is.null(dots$bandwidth.compute)) TRUE else dots$bandwidth.compute
    ))
      stop("bounded adaptive_nn remains unsupported for npudistbw() in npRmpi",
           call. = FALSE)
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    t.names <- NULL
    if(!is.data.frame(dat) && !is.matrix(dat))
      t.names <- paste(deparse(substitute(dat)), collapse = "")

    dat = toFrame(dat)

    if(!is.null(t.names))
      names(dat) <- t.names

    bws = double(dim(dat)[2])

    tbw <- npudistbw.default(dat = dat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw
  }

npudistbw.dbandwidth <-
  function(dat = stop("invoked without input data 'dat'"),
           bws,
           gdat = NULL,
           bandwidth.compute = TRUE,
           cfac.dir = 2.5*(3.0-sqrt(5)),
           scale.factor.init = 0.5,
           dfac.dir = 0.25*(3.0-sqrt(5)),
           dfac.init = 0.375,
           dfc.dir = 3,
           do.full.integral = FALSE,
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
           memfac = 500.0,
           ngrid = 100,
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
    do.full.integral <- npValidateScalarLogical(do.full.integral, "do.full.integral")
    scale.init.categorical.sample <-
      npValidateScalarLogical(scale.init.categorical.sample, "scale.init.categorical.sample")
    transform.bounds <- npValidateScalarLogical(transform.bounds, "transform.bounds")
    itmax <- npValidatePositiveInteger(itmax, "itmax")
    ngrid <- npValidatePositiveInteger(ngrid, "ngrid")
    ftol <- npValidatePositiveFiniteNumeric(ftol, "ftol")
    tol <- npValidatePositiveFiniteNumeric(tol, "tol")
    small <- npValidatePositiveFiniteNumeric(small, "small")
    memfac <- npValidatePositiveFiniteNumeric(memfac, "memfac")
    penalty.multiplier <- npValidatePositiveFiniteNumeric(penalty.multiplier, "penalty.multiplier")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    if (!missing(nmulti))
      nmulti <- npValidateNmulti(nmulti)
    .npRmpi_require_active_slave_pool(where = "npudistbw()")
    if (.npRmpi_npudistbw_bounded_adaptive_requested(
      bwtype = bws$type,
      ckerbound = bws$ckerbound,
      bandwidth.compute = bandwidth.compute
    ))
      stop("bounded adaptive_nn remains unsupported for npudistbw() in npRmpi",
           call. = FALSE)
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    dat = toFrame(dat)

    nofi <- missing(do.full.integral)
    nogi <- missing(ngrid)

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

    if(any(bws$iuno))
      stop("distribution bandwidth selection does not support unordered data types")

    dat <- na.omit(dat)
    rows.omit <- unclass(na.action(dat))

    nrow = dim(dat)[1]
    ncol = dim(dat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    odat <- dat
    dat = toMatrix(dat)

    duno = dat[, bws$iuno, drop = FALSE]
    dcon = dat[, bws$icon, drop = FALSE]
    dord = dat[, bws$iord, drop = FALSE]

    tbw <- bws

    mysd <- EssDee(dcon)
    nconfac <- nrow^(-1.0/(1.0+bws$ckerorder))
    ncatfac <- nrow^(-2.0/(1.0+bws$ckerorder))

    ## these are the points where we evaluate the CDF
    ## for now we default to the training points (hence cdf_on_train = TRUE below)
    if(!is.null(gdat)){
      gdat <- toFrame(gdat)
      if(any(is.na(gdat)))
        stop("na's not allowed to be present in cdf gdata")
      gdat <- toMatrix(gdat)

      guno = gdat[, bws$iuno, drop = FALSE]
      gord = gdat[, bws$iord, drop = FALSE]
      gcon = gdat[, bws$icon, drop = FALSE]
      cdf_on_train = FALSE
      nog = nrow(gdat)

    } else {
      if(do.full.integral) {
        cdf_on_train = TRUE
        nog = 0
        guno = data.frame()
        gord = data.frame()
        gcon = data.frame()
      } else {
        cdf_on_train = FALSE
        nog = ngrid
        probs <- seq(0,1,length.out = nog)
        ev <- odat[seq_len(nog),,drop = FALSE]
        for (i in seq_len(ncol(ev))) {
          ev[,i] <- uocquantile(odat[,i], probs)
        }

        ev <- toMatrix(ev)

        guno = ev[, bws$iuno, drop = FALSE]
        gord = ev[, bws$iord, drop = FALSE]
        gcon = ev[, bws$icon, drop = FALSE]
      }
    }

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

    if (bandwidth.compute){
      cont.start <- npContinuousSearchStartControls(
        scale.factor.init.lower,
        scale.factor.init.upper,
        scale.factor.init,
        scale.factor.search.lower,
        where = "npudistbw"
      )
      myopti = list(num_obs_train = dim(dat)[1],
        num_obs_eval = nog,
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
          cv.cdf = DBWM_CVLS),
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
        cdf_on_train = cdf_on_train,
        nuno = dim(duno)[2],
        nord = dim(dord)[2],
        ncon = dim(dcon)[2],
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
        nconfac = nconfac, ncatfac = ncatfac, memfac = memfac,
        scale.factor.lower.bound = scale.factor.search.lower)
      cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

      if (bws$method != "normal-reference"){
        myout <-
          .Call("C_np_distribution_bw",
                as.double(duno), as.double(dord), as.double(dcon),
                as.double(guno), as.double(gord), as.double(gcon), as.double(mysd),
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
        gbw = bws$ncon
        if (gbw > 0){
          gbw_idx <- seq_len(gbw)
          nbw[gbw_idx] = 1.587
          if(!bws$scaling)
            nbw[gbw_idx]=nbw[gbw_idx]*mysd*nconfac
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

    tbw <- dbandwidth(bw = tbw$bw,
                      bwmethod = tbw$method,
                      bwscaling = tbw$scaling,
                      bwtype = tbw$type,
                      ckertype = tbw$ckertype,
                      ckerorder = tbw$ckerorder,
                      ckerbound = tbw$ckerbound,
                      ckerlb = tbw$ckerlb,
                      ckerub = tbw$ckerub,
                      ukertype = c("aitchisonaitken"),
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

npudistbw.default <-
  function(dat = stop("invoked without input data 'dat'"),
           bws,
           gdat,
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
           do.full.integral,
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
           memfac,
           ngrid,
           nmulti,
           okertype,
           penalty.multiplier,
           remin,
           scale.init.categorical.sample,
           scale.factor.search.lower = NULL,
           small,
           tol,
           transform.bounds,
           ...){
    .npRmpi_require_active_slave_pool(where = "npudistbw()")
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
    if (!missing(okertype)) bw.args$okertype <- okertype
    tbw <- do.call(dbandwidth, bw.args)


    ## next grab dummies for actual bandwidth selection and perform call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("gdat","bandwidth.compute", "nmulti", "remin", "itmax",
               "do.full.integral", "ngrid", "ftol", "tol",
               "small", "lbc.dir", "dfc.dir", "cfac.dir", "initc.dir",
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.init.categorical.sample", "scale.factor.search.lower", "memfac",
               "transform.bounds",
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
      "Selecting distribution bandwidth",
      do.call(npudistbw.dbandwidth, bwsel.args)
    )

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }
