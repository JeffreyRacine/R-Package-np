npcdistbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npcdistbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat", "gydat"),
                                     eval_env = parent.frame())
    UseMethod("npcdistbw", target)
  }

npcdistbw.formula <-
  function(formula, data, subset, na.action, call, gdata = NULL, ...){
    orig.ts <- tryCatch({
        if (missing(data))
            .np_terms_ts_mask(terms_obj = terms(formula),
                              data = environment(formula),
                              eval_env = environment(formula))
        else .np_terms_ts_mask(terms_obj = terms(formula, data = data),
                               data = data,
                               eval_env = environment(formula))
    }, error = function(e) FALSE)

    has.gval <- !is.null(gdata)

    gmf <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    gm <- match(c("formula", "gdata"),
               names(gmf), nomatch = 0)

    mf <- mf[c(1,m)]
    gmf <- gmf[c(1,gm)]

    formula.call <- .np_bw_formula_from_call(call_obj = call, eval_env = parent.frame())
    if (!is.null(formula.call)) {
      mf[[2]] <- formula.call
      gmf[[2]] <- formula.call
    }


    mf[[1]] <- as.name("model.frame")
    gmf[[1]] <- as.name("model.frame")
    formula.obj <- .np_bw_resolve_formula(formula_obj = formula,
                                        formula_call = formula.call,
                                        eval_env = parent.frame())

    variableNames <- if(m[2] > 0) explodeFormula(formula.obj, data = data) else explodeFormula(formula.obj)

    ## make formula evaluable, then eval
    varsPlus <- lapply(variableNames, paste, collapse=" + ")
    mf[["formula"]] <- as.formula(paste(" ~ ", varsPlus[[1]]," + ",
                                        varsPlus[[2]]),
                                  env = environment(formula))
    gmf[["formula"]] <- mf[["formula"]]

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

    ydat <- mf[, variableNames[[1]], drop = FALSE]
    xdat <- mf[, variableNames[[2]], drop = FALSE]

    if (has.gval) {
      gmf.args <- as.list(gmf[-1L])
      names(gmf.args)[names(gmf.args) == "gdata"] <- "data"
      gmf <- do.call(stats::model.frame, gmf.args, envir = parent.frame())
      gydat <- gmf[, variableNames[[1]], drop = FALSE]
    }

    dots <- list(...)
    seed.args <- c(list(xdat = xdat, ydat = ydat),
                   if (has.gval) list(gydat = gydat) else list(),
                   dots)
    tbw <- do.call(npcdistbw, seed.args)

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$variableNames <- variableNames

    tbw
  }

.npRmpi_with_local_cdist_eval <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  old.local <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.autodispatch.context = TRUE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old.local), add = TRUE)
  old.mode <- .Call("C_np_set_local_regression_mode", TRUE, PACKAGE = "npRmpi")
  on.exit(.Call("C_np_set_local_regression_mode", old.mode, PACKAGE = "npRmpi"), add = TRUE)
  force(expr)
}

.npRmpi_npcdistbw_bounded_adaptive_requested <- function(bwtype = NULL,
                                                         cxkerbound = "none",
                                                         cykerbound = "none",
                                                         bandwidth.compute = TRUE) {
  bwtype.ch <- tolower(as.character(bwtype)[1L])
  cxkerbound.ch <- tolower(as.character(cxkerbound)[1L])
  cykerbound.ch <- tolower(as.character(cykerbound)[1L])
  isTRUE(bandwidth.compute) &&
    identical(bwtype.ch, "adaptive_nn") &&
    (!identical(cxkerbound.ch, "none") || !identical(cykerbound.ch, "none"))
}

npcdistbw.condbandwidth <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           gydat = NULL,
           bws,
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
    ydat = toFrame(ydat)
    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(dim(ydat)[2]+dim(xdat)[2])
    }
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
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)
    if (length(bws$ybw) != dim(ydat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of", "'ydat'"))

    if (length(bws$xbw) != dim(xdat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of", "'xdat'"))

    if (dim(ydat)[1] != dim(xdat)[1])
      stop(paste("number of rows of", "'ydat'", "does not match", "'xdat'"))

    if ((any(bws$iycon) &&
         !all(vapply(as.data.frame(ydat[, bws$iycon]), inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iyord) &&
         !all(vapply(as.data.frame(ydat[, bws$iyord]), inherits, logical(1), "ordered"))) ||
        (any(bws$iyuno) &&
         !all(vapply(as.data.frame(ydat[, bws$iyuno]), inherits, logical(1), "factor"))))
      stop(paste("supplied bandwidths do not match", "'ydat'", "in type"))

    if ((any(bws$ixcon) &&
         !all(vapply(as.data.frame(xdat[, bws$ixcon]), inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$ixord) &&
         !all(vapply(as.data.frame(xdat[, bws$ixord]), inherits, logical(1), "ordered"))) ||
        (any(bws$ixuno) &&
         !all(vapply(as.data.frame(xdat[, bws$ixuno]), inherits, logical(1), "factor"))))
      stop(paste("supplied bandwidths do not match", "'xdat'", "in type"))

    ##if (bws$type != 'fixed')
    ##stop("only fixed bandwidths currently supported with ccdf bandwidth selection")

    ## catch and destroy NA's
    goodrows <- seq_len(nrow(xdat))
    rows.omit <- unclass(na.action(na.omit(data.frame(xdat,ydat))))
    goodrows[rows.omit] <- 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    spec <- npCanonicalConditionalRegSpec(
      regtype = if (is.null(bws$regtype)) "lc" else bws$regtype,
      basis = if (is.null(bws$basis)) "glp" else bws$basis,
      degree = if (is.null(bws$degree)) NULL else bws$degree,
      bernstein.basis = isTRUE(bws$bernstein.basis),
      ncon = bws$xncon,
      where = "npcdistbw"
    )
    .npRmpi_require_active_slave_pool(where = "npcdistbw()")
    if (.npRmpi_npcdistbw_bounded_adaptive_requested(
      bwtype = bws$type,
      cxkerbound = bws$cxkerbound,
      cykerbound = bws$cykerbound,
      bandwidth.compute = bandwidth.compute
    ))
      stop("bounded adaptive_nn remains unsupported for npcdistbw() in npRmpi",
           call. = FALSE)
    use.local.compiled.adaptive.cvls <- bandwidth.compute &&
      identical(bws$method, "cv.ls") &&
      identical(bws$type, "adaptive_nn")
    keep_local_cvls_nn <- bandwidth.compute &&
      identical(bws$method, "cv.ls") &&
      (use.local.compiled.adaptive.cvls ||
       (identical(spec$regtype.engine, "lp") &&
        identical(bws$type, "generalized_nn")))
    keep_local_raw_degree1_cvls <- bandwidth.compute &&
      identical(bws$method, "cv.ls") &&
      identical(bws$type, "fixed") &&
      npIsRawDegreeOneConditionalSpec(spec, bws$xncon)
    if (.npRmpi_autodispatch_active() &&
        !keep_local_cvls_nn &&
        !keep_local_raw_degree1_cvls)
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows,,drop = FALSE]


    nrow = nrow(ydat)
    yncol = ncol(ydat)
    xncol = ncol(xdat)

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    oydat <- ydat

    ydat = toMatrix(ydat)

    yuno = ydat[, bws$iyuno, drop = FALSE]
    ycon = ydat[, bws$iycon, drop = FALSE]
    yord = ydat[, bws$iyord, drop = FALSE]

    xdat = toMatrix(xdat)

    xuno = xdat[, bws$ixuno, drop = FALSE]
    xcon = xdat[, bws$ixcon, drop = FALSE]
    xord = xdat[, bws$ixord, drop = FALSE]

    tbw <- bws
    spec <- npCanonicalConditionalRegSpec(
      regtype = if (is.null(tbw$regtype)) "lc" else tbw$regtype,
      basis = if (is.null(tbw$basis)) "glp" else tbw$basis,
      degree = if (is.null(tbw$degree)) NULL else tbw$degree,
      bernstein.basis = isTRUE(tbw$bernstein.basis),
      ncon = tbw$xncon,
      where = "npcdistbw"
    )
    tbw$regtype <- spec$regtype
    tbw$pregtype <- switch(spec$regtype,
                           lc = "Local-Constant",
                           ll = "Local-Linear",
                           lp = "Local-Polynomial")
    tbw$basis <- spec$basis
    tbw$degree <- spec$degree
    tbw$bernstein.basis <- spec$bernstein.basis
    tbw$regtype.engine <- spec$regtype.engine
    tbw$basis.engine <- spec$basis.engine
    tbw$degree.engine <- spec$degree.engine
    tbw$bernstein.basis.engine <- spec$bernstein.basis.engine
    reg.code <- if (identical(spec$regtype.engine, "lp")) REGTYPE_LP else REGTYPE_LC
    degree.code <- if (tbw$xncon > 0L) as.integer(spec$degree.engine) else integer(0)
    basis.code <- as.integer(npLpBasisCode(spec$basis.engine))
    bernstein.engine <- isTRUE(spec$bernstein.basis.engine)

    if(!is.null(gydat)){
      gydat <- toFrame(gydat)
      if(any(is.na(gydat)))
        stop("na's not allowed to be present in cdf gdata")

      gydat <- toMatrix(gydat)

      gyuno = gydat[, bws$iyuno, drop = FALSE]
      gyord = gydat[, bws$iyord, drop = FALSE]
      gycon = gydat[, bws$iycon, drop = FALSE]
      cdf_on_train = FALSE
      nog = nrow(gydat)

    } else {
      if(do.full.integral) {
        cdf_on_train = TRUE
        nog = 0

        gyuno = data.frame()
        gyord = data.frame()
        gycon = data.frame()

      } else {
        cdf_on_train = FALSE
        nog = ngrid
        probs <- seq(0,1,length.out = nog)

        evy <- oydat[seq_len(nog),,drop = FALSE]
        for (i in seq_len(ncol(evy))) {
          evy[,i] <- uocquantile(oydat[,i], probs)
        }

        evy <- toMatrix(evy)

        gyuno = evy[, bws$iyuno, drop = FALSE]
        gyord = evy[, bws$iyord, drop = FALSE]
        gycon = evy[, bws$iycon, drop = FALSE]

      }

    }

    mysd <- EssDee(data.frame(xcon,ycon))
    nconfac <- nrow^(-1.0/(2.0*bws$cxkerorder+bws$ncon))
    ncatfac <- nrow^(-2.0/(2.0*bws$cxkerorder+bws$ncon))

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

    if (bandwidth.compute){
      cont.start <- npContinuousSearchStartControls(
        scale.factor.init.lower,
        scale.factor.init.upper,
        scale.factor.init,
        scale.factor.search.lower,
        where = "npcdistbw"
      )
      myopti = list(num_obs_train = nrow,
        num_obs_grid = nog,
        iMultistart = IMULTI_TRUE,
        iNum_Multistart = nmulti,
        int_use_starting_values = (if (all(bws$ybw==0) && all(bws$xbw==0)) USE_START_NO else USE_START_YES),
        int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
        BANDWIDTH_den_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=(if (remin) RE_MIN_TRUE else RE_MIN_FALSE),
        int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
        bwmethod = switch(bws$method,
          cv.ls = CDBWM_CVLS),
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
        ynuno = dim(yuno)[2],
        ynord = dim(yord)[2],
        yncon = dim(ycon)[2],
        xnuno = dim(xuno)[2],
        xnord = dim(xord)[2],
        xncon = dim(xcon)[2],
        cdf_on_train = cdf_on_train,
        int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
        scale.init.categorical.sample=scale.init.categorical.sample,
        dfc.dir = dfc.dir,
        transform.bounds = transform.bounds)

      myoptd = list(ftol=ftol, tol=tol, small=small, memfac = memfac,
        lbc.dir = lbc.dir, cfac.dir = cfac.dir, initc.dir = initc.dir,
        lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir,
        lbc.init = cont.start$scale.factor.init.lower,
        hbc.init = cont.start$scale.factor.init.upper,
        cfac.init = cont.start$scale.factor.init,
        lbd.init = lbd.init, hbd.init = hbd.init, dfac.init = dfac.init,
        nconfac = nconfac, ncatfac = ncatfac,
        scale.factor.lower.bound = scale.factor.search.lower)

      cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
      cyker.bounds.c <- npKernelBoundsMarshal(bws$cykerlb[bws$iycon], bws$cykerub[bws$iycon])

      if (bws$method != "normal-reference"){
        myout <- npWithLocalLinearRawBasisSearchError(
          if (keep_local_cvls_nn || keep_local_raw_degree1_cvls) {
            .npRmpi_with_local_cdist_eval(
              .Call("C_np_distribution_conditional_bw",
                    as.double(yuno), as.double(yord), as.double(ycon),
                    as.double(xuno), as.double(xord), as.double(xcon),
                    as.double(gyuno), as.double(gyord), as.double(gycon),
                    as.double(mysd),
                    as.integer(myopti), as.double(myoptd),
                    as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                                bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                                bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
                    as.integer(nmulti),
                    as.integer(penalty_mode),
                    as.double(penalty.multiplier),
                    as.integer(degree.code),
                    as.integer(bernstein.engine),
                    as.integer(basis.code),
                    as.integer(reg.code),
                    as.double(cxker.bounds.c$lb),
                    as.double(cxker.bounds.c$ub),
                    as.double(cyker.bounds.c$lb),
                    as.double(cyker.bounds.c$ub),
                    PACKAGE="npRmpi")
            )
          } else {
            .Call("C_np_distribution_conditional_bw",
                  as.double(yuno), as.double(yord), as.double(ycon),
                  as.double(xuno), as.double(xord), as.double(xcon),
                  as.double(gyuno), as.double(gyord), as.double(gycon),
                  as.double(mysd),
                  as.integer(myopti), as.double(myoptd),
                  as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                              bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                              bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
                  as.integer(nmulti),
                  as.integer(penalty_mode),
                  as.double(penalty.multiplier),
                  as.integer(degree.code),
                  as.integer(bernstein.engine),
                  as.integer(basis.code),
                  as.integer(reg.code),
                  as.double(cxker.bounds.c$lb),
                  as.double(cxker.bounds.c$ub),
                  as.double(cyker.bounds.c$lb),
                  as.double(cyker.bounds.c$ub),
                  PACKAGE="npRmpi")
          },
          where = "npcdistbw",
          spec = spec,
          bwmethod = bws$method,
          ncon = tbw$xncon
        )
        total.time <- proc.time()[3] - elapsed.start
      } else {
        nbw = double(yncol+xncol)
        gbw = bws$yncon+bws$xncon
        if (gbw > 0){
          xcon_idx <- seq_len(bws$xncon)
          ycon_idx <- seq.int(from = bws$xncon + 1L, length.out = bws$yncon)
          if (length(xcon_idx) > 0L)
            nbw[xcon_idx] <- 1.06
          if (length(ycon_idx) > 0L)
            nbw[ycon_idx] <- 1.587
          if(!bws$scaling){
            gbw_idx <- seq_len(gbw)
            nbw[gbw_idx]=nbw[gbw_idx]*mysd*nconfac
          }
        }
        myout= list( bw = nbw, fval = c(NA,NA,NA) )
        total.time <- NA
      }

      yr = seq_len(yncol)
      xr = seq_len(xncol)
      rorder = numeric(yncol + xncol)

      ## bandwidths are passed back from the C routine in an unusual order
      ## xc, y[cuo], x[uo]

      rxcon = xr[bws$ixcon]
      rxuno = xr[bws$ixuno]
      rxord = xr[bws$ixord]

      rycon = yr[bws$iycon]
      ryuno = yr[bws$iyuno]
      ryord = yr[bws$iyord]


      ## rorder[c(rxcon,rycon,ryuno,ryord,rxuno,rxord)]=1:(yncol+xncol)

      tbw <- bws
      tbw$ybw[c(rycon,ryuno,ryord)] <- myout$bw[yr+bws$xncon]
      tbw$xbw[c(rxcon,rxuno,rxord)] <- myout$bw[setdiff(seq_len(yncol + xncol), yr + bws$xncon)]

      tbw$fval = myout$fval[1]
      tbw$ifval = myout$fval[2]
      tbw$initial.fval <- if (length(myout$fval) >= 3L) myout$fval[3L] else NA_real_
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
      tbw$fval.history <- myout$fval.history
      tbw$eval.history <- myout$eval.history
      tbw$invalid.history <- myout$invalid.history
      tbw$timing <- myout$timing
      tbw$total.time <- total.time
    }

    ## bandwidth metadata
    tbw$sfactor <- tbw$bandwidth <- list(x = tbw$xbw, y = tbw$ybw)

    apply_bw_meta <- function(tl, dfactor){
      for (nm in names(tl)) {
        idx <- tl[[nm]]
        if (length(idx) == 0L)
          next
        if (tbw$scaling) {
          tbw$bandwidth[[nm]][idx] <- tbw$bandwidth[[nm]][idx] * dfactor[[nm]]
        } else {
          tbw$sfactor[[nm]][idx] <- tbw$sfactor[[nm]][idx] / dfactor[[nm]]
        }
      }
    }

    if ((tbw$xnuno+tbw$ynuno) > 0){
      dfactor <- ncatfac
      dfactor <- list(x = dfactor, y = dfactor)

      tl <- list(x = tbw$xdati$iuno, y = tbw$ydati$iuno)

      apply_bw_meta(tl = tl, dfactor = dfactor)
    }

    if ((tbw$xnord+tbw$ynord) > 0){
      dfactor <- ncatfac
      dfactor <- list(x = dfactor, y = dfactor)

      tl <- list(x = tbw$xdati$iord, y = tbw$ydati$iord)

      apply_bw_meta(tl = tl, dfactor = dfactor)
    }


    if (tbw$ncon > 0){
      dfactor <- nconfac
      dfactor <- list(x = EssDee(xcon)*dfactor, y = EssDee(ycon)*dfactor)

      tl <- list(x = tbw$xdati$icon, y = tbw$ydati$icon)

      apply_bw_meta(tl = tl, dfactor = dfactor)
    }

    initial.fval <- tbw$initial.fval
    tbw <- condbandwidth(xbw = tbw$xbw,
                         ybw = tbw$ybw,
                         bwmethod = tbw$method,
                         bwscaling = tbw$scaling,
                         bwtype = tbw$type,
                         cxkertype = tbw$cxkertype,
                         cxkerorder = tbw$cxkerorder,
                         cxkerbound = tbw$cxkerbound,
                         cxkerlb = tbw$cxkerlb,
                         cxkerub = tbw$cxkerub,
                         uxkertype = tbw$uxkertype,
                         oxkertype = tbw$oxkertype,
                         cykertype = tbw$cykertype,
                         cykerorder = tbw$cykerorder,
                         cykerbound = tbw$cykerbound,
                         cykerlb = tbw$cykerlb,
                         cykerub = tbw$cykerub,
                         uykertype = tbw$uykertype,
                         oykertype = tbw$oykertype,
                        fval = tbw$fval,
                        ifval = tbw$ifval,
                        num.feval = tbw$num.feval,
                        num.feval.fast = tbw$num.feval.fast,
                        fval.history = tbw$fval.history,
                         eval.history = tbw$eval.history,
                         invalid.history = tbw$invalid.history,
                         nobs = tbw$nobs,
                         xdati = tbw$xdati,
                         ydati = tbw$ydati,
                         xnames = tbw$xnames,
                         ynames = tbw$ynames,
                         sfactor = tbw$sfactor,
                         bandwidth = tbw$bandwidth,
                         rows.omit = rows.omit,
                         nconfac = nconfac,
                         ncatfac = ncatfac,
                         sdev = mysd,
                         bandwidth.compute = bandwidth.compute,
                         timing = tbw$timing,
                         total.time = tbw$total.time,
                         regtype = tbw$regtype,
                         pregtype = tbw$pregtype,
                         basis = tbw$basis,
                         degree = tbw$degree,
                         bernstein.basis = tbw$bernstein.basis,
                         regtype.engine = tbw$regtype.engine,
                         basis.engine = tbw$basis.engine,
                         degree.engine = tbw$degree.engine,
                         bernstein.basis.engine = tbw$bernstein.basis.engine)
    tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)

    tbw <- .np_refresh_xy_bandwidth_metadata(tbw)
    tbw$initial.fval <- if (!is.null(initial.fval)) initial.fval else NA_real_

    tbw
  }

.npcdistbw_build_condbandwidth <- function(xdat,
                                           ydat,
                                           bws,
                                           bandwidth.compute,
                                           reg.args) {
  x.info <- untangle(xdat)
  y.info <- untangle(ydat)
  y.idx <- seq_len(ncol(ydat))
  x.idx <- seq_len(ncol(xdat))

  bw.args <- c(
    list(
      xbw = bws[length(y.idx) + x.idx],
      ybw = bws[y.idx],
      nobs = nrow(xdat),
      xdati = x.info,
      ydati = y.info,
      xnames = names(xdat),
      ynames = names(ydat),
      bandwidth.compute = bandwidth.compute
    ),
    reg.args
  )

  out <- do.call(condbandwidth, bw.args)
  if (!is.null(reg.args$scale.factor.search.lower))
    out$scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      reg.args$scale.factor.search.lower
    )
  out
}

.npcdistbw_run_fixed_degree <- function(xdat, ydat, bws, reg.args, opt.args) {
  tbw <- .npcdistbw_build_condbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args
  )

  do.call(npcdistbw.condbandwidth, c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args))
}

.npcdistbw_run_fixed_degree_bcast_payload <- function(xdat, ydat, bws, reg.args, opt.args) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.messages <- getOption("np.messages")
  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)

  options(npRmpi.autodispatch.disable = TRUE)
  if (!isTRUE(rank == 0L))
    options(np.messages = FALSE)

  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(np.messages = old.messages), add = TRUE)

  .npcdistbw_run_fixed_degree(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    reg.args = reg.args,
    opt.args = opt.args
  )
}

.npcdistbw_run_fixed_degree_collective <- function(xdat,
                                                   ydat,
                                                   bws,
                                                   reg.args,
                                                   opt.args,
                                                   comm = 1L) {
  if (.npRmpi_has_active_slave_pool(comm = comm) &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    mc <- substitute(
      get(".npcdistbw_run_fixed_degree_bcast_payload", envir = asNamespace("npRmpi"), inherits = FALSE)(
        XDAT,
        YDAT,
        BWS,
        REGARGS,
        OPTARGS
      ),
      list(
        XDAT = xdat,
        YDAT = ydat,
        BWS = bws,
        REGARGS = reg.args,
        OPTARGS = opt.args
      )
    )
    return(.npRmpi_bcast_cmd_expr(mc, comm = comm, caller.execute = TRUE))
  }

  .npcdistbw_run_fixed_degree(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    reg.args = reg.args,
    opt.args = opt.args
  )
}

.npcdistbw_eval_only <- function(xdat,
                                 ydat,
                                 gydat = NULL,
                                 bws,
                                 do.full.integral = FALSE,
                                 ngrid = 100L,
                                 invalid.penalty = c("baseline", "dbmax"),
                                 penalty.multiplier = 10,
                                 force.local = TRUE) {
  invalid.penalty <- match.arg(invalid.penalty)

  ydat <- toFrame(ydat)
  xdat <- toFrame(xdat)

  if (length(bws$ybw) != dim(ydat)[2])
    stop("length of bandwidth vector does not match number of columns of 'ydat'")
  if (length(bws$xbw) != dim(xdat)[2])
    stop("length of bandwidth vector does not match number of columns of 'xdat'")
  if (dim(ydat)[1] != dim(xdat)[1])
    stop("number of rows of 'ydat' does not match 'xdat'")

  goodrows <- seq_len(nrow(xdat))
  rows.omit <- unclass(na.action(na.omit(data.frame(xdat, ydat))))
  goodrows[rows.omit] <- 0

  xdat <- xdat[goodrows,, drop = FALSE]
  ydat <- ydat[goodrows,, drop = FALSE]

  oydat <- ydat
  ymat <- toMatrix(ydat)
  xmat <- toMatrix(xdat)

  yuno <- ymat[, bws$iyuno, drop = FALSE]
  ycon <- ymat[, bws$iycon, drop = FALSE]
  yord <- ymat[, bws$iyord, drop = FALSE]
  xuno <- xmat[, bws$ixuno, drop = FALSE]
  xcon <- xmat[, bws$ixcon, drop = FALSE]
  xord <- xmat[, bws$ixord, drop = FALSE]

  if (!is.null(gydat)) {
    gydat <- toFrame(gydat)
    if (any(is.na(gydat)))
      stop("na's not allowed to be present in cdf gdata")
    gmat <- toMatrix(gydat)
    gyuno <- gmat[, bws$iyuno, drop = FALSE]
    gyord <- gmat[, bws$iyord, drop = FALSE]
    gycon <- gmat[, bws$iycon, drop = FALSE]
    cdf_on_train <- FALSE
    nog <- nrow(gmat)
  } else if (isTRUE(do.full.integral)) {
    cdf_on_train <- TRUE
    nog <- 0L
    gyuno <- data.frame()
    gyord <- data.frame()
    gycon <- data.frame()
  } else {
    cdf_on_train <- FALSE
    nog <- npValidatePositiveInteger(ngrid, "ngrid")
    probs <- seq(0, 1, length.out = nog)
    evy <- oydat[seq_len(nog),, drop = FALSE]
    for (i in seq_len(ncol(evy)))
      evy[, i] <- uocquantile(oydat[, i], probs)
    evy <- toMatrix(evy)
    gyuno <- evy[, bws$iyuno, drop = FALSE]
    gyord <- evy[, bws$iyord, drop = FALSE]
    gycon <- evy[, bws$iycon, drop = FALSE]
  }

  mysd <- EssDee(data.frame(xcon, ycon))
  nrow <- nrow(ymat)
  nconfac <- nrow^(-1.0 / (2.0 * bws$cxkerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$cxkerorder + bws$ncon))

  penalty_mode <- if (invalid.penalty == "baseline") 1L else 0L
  reg.code <- if (identical(bws$regtype.engine, "lp")) REGTYPE_LP else REGTYPE_LC
  degree.code <- if (bws$xncon > 0L) as.integer(bws$degree.engine) else integer(0L)
  basis.code <- as.integer(npLpBasisCode(bws$basis.engine))
  bernstein.engine <- isTRUE(bws$bernstein.basis.engine)

  myopti <- list(
    num_obs_train = nrow,
    num_obs_grid = nog,
    iMultistart = IMULTI_FALSE,
    iNum_Multistart = 0L,
    int_use_starting_values = USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = 0L,
    int_RESTART_FROM_MIN = RE_MIN_FALSE,
    int_MINIMIZE_IO = IO_MIN_TRUE,
    bwmethod = CDBWM_CVLS,
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
      racineliyan = OKER_RLY),
    oykerneval = switch(bws$oykertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY),
    ynuno = dim(yuno)[2],
    ynord = dim(yord)[2],
    yncon = dim(ycon)[2],
    xnuno = dim(xuno)[2],
    xnord = dim(xord)[2],
    xncon = dim(xcon)[2],
    cdf_on_train = cdf_on_train,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    scale.init.categorical.sample = FALSE,
    dfc.dir = 0L,
    transform.bounds = FALSE
  )

  myoptd <- list(
    ftol = 0,
    tol = 0,
    small = 0,
    memfac = 0,
    lbc.dir = 0,
    cfac.dir = 0,
    initc.dir = 0,
    lbd.dir = 0,
    hbd.dir = 0,
    dfac.dir = 0,
    initd.dir = 0,
    lbc.init = 0,
    hbc.init = 0,
    cfac.init = 0,
    lbd.init = 0,
    hbd.init = 0,
    dfac.init = 0,
    nconfac = nconfac,
    ncatfac = ncatfac,
    scale.factor.lower.bound = npResolveScaleFactorLowerBound(
      npGetScaleFactorSearchLower(bws),
      argname = "bws$scale.factor.search.lower"
    )
  )

  cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
  cyker.bounds.c <- npKernelBoundsMarshal(bws$cykerlb[bws$iycon], bws$cykerub[bws$iycon])

  eval_call <- function() {
    .Call(
      "C_np_distribution_conditional_bw_eval",
      as.double(yuno),
      as.double(yord),
      as.double(ycon),
      as.double(xuno),
      as.double(xord),
      as.double(xcon),
      as.double(gyuno),
      as.double(gyord),
      as.double(gycon),
      as.double(mysd),
      as.integer(myopti),
      as.double(myoptd),
      as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                  bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                  bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
      as.integer(1L),
      as.integer(penalty_mode),
      as.double(penalty.multiplier),
      as.integer(degree.code),
      as.integer(bernstein.engine),
      as.integer(basis.code),
      as.integer(reg.code),
      as.double(cxker.bounds.c$lb),
      as.double(cxker.bounds.c$ub),
      as.double(cyker.bounds.c$lb),
      as.double(cyker.bounds.c$ub),
      PACKAGE = "npRmpi"
    )
  }
  out <- if (isTRUE(force.local) &&
             .npRmpi_has_active_slave_pool(comm = 1L) &&
             !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    .npRmpi_with_local_cdist_eval(eval_call())
  } else {
    eval_call()
  }

  list(
    objective = as.numeric(out$fval[1L]),
    num.feval = 1L,
    num.feval.fast = as.numeric(as.numeric(out$fast.history[1L]) > 0)
  )
}

.npcdistbw_run_fixed_degree <- function(xdat, ydat, bws, reg.args, opt.args) {
  tbw <- .npcdistbw_build_condbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args
  )

  do.call(npcdistbw.condbandwidth, c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args))
}

.npcdistbw_nomad_bw_setup <- function(xdat,
                                      ydat,
                                      template,
                                      bandwidth.scale.categorical = 1e4) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  xmat <- toMatrix(xdat)
  ymat <- toMatrix(ydat)
  xcon <- xmat[, template$ixcon, drop = FALSE]
  ycon <- ymat[, template$iycon, drop = FALSE]
  nrow <- nrow(xmat)
  nconfac <- nrow^(-1.0 / (2.0 * template$cxkerorder + template$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * template$cxkerorder + template$ncon))

  x_offset <- length(template$ybw)
  y_cont_flat <- which(template$iycon)
  x_cont_flat <- x_offset + which(template$ixcon)
  y_uno_flat <- which(template$iyuno)
  y_ord_flat <- which(template$iyord)
  x_uno_flat <- x_offset + which(template$ixuno)
  x_ord_flat <- x_offset + which(template$ixord)

  cat_upper_one <- function(values, kernel) {
    if (identical(kernel, "aitchisonaitken")) {
      nlev <- length(unique(values))
      return((nlev - 1) / nlev)
    }
    1
  }

  cat_upper <- c(
    if (length(y_uno_flat)) vapply(which(template$iyuno), function(i) cat_upper_one(ydat[[i]], template$uykertype), numeric(1L)) else numeric(0L),
    rep.int(1, length(y_ord_flat)),
    if (length(x_uno_flat)) vapply(which(template$ixuno), function(i) cat_upper_one(xdat[[i]], template$uxkertype), numeric(1L)) else numeric(0L),
    rep.int(1, length(x_ord_flat))
  )

  list(
    cont_flat = c(y_cont_flat, x_cont_flat),
    cont_scale = c(EssDee(ycon), EssDee(xcon)) * nconfac,
    cat_flat = c(y_uno_flat, y_ord_flat, x_uno_flat, x_ord_flat),
    ncatfac = ncatfac,
    bandwidth.scale.categorical = bandwidth.scale.categorical,
    cat_upper = cat_upper
  )
}

.npcdistbw_nomad_point_to_bw <- function(point, template, setup) {
  point <- as.numeric(point)
  ncont <- length(setup$cont_flat)
  ncat <- length(setup$cat_flat)
  bws <- numeric(length(template$ybw) + length(template$xbw))

  if (ncont > 0L) {
    gamma <- point[seq_len(ncont)]
    ext_bw <- gamma * setup$cont_scale
    bws[setup$cont_flat] <- if (isTRUE(template$scaling)) gamma else ext_bw
  }

  if (ncat > 0L) {
    lambda_scaled <- point[ncont + seq_len(ncat)]
    ext_bw <- lambda_scaled / setup$bandwidth.scale.categorical
    bws[setup$cat_flat] <- if (isTRUE(template$scaling)) ext_bw / setup$ncatfac else ext_bw
  }

  bws
}

.npcdistbw_nomad_bw_to_point <- function(bws, template, setup) {
  point <- numeric(length(setup$cont_flat) + length(setup$cat_flat))

  if (length(setup$cont_flat) > 0L) {
    raw <- bws[setup$cont_flat]
    point[seq_along(setup$cont_flat)] <- if (isTRUE(template$scaling)) {
      raw
    } else {
      raw / setup$cont_scale
    }
  }

  if (length(setup$cat_flat) > 0L) {
    raw <- bws[setup$cat_flat]
    ext_bw <- if (isTRUE(template$scaling)) raw * setup$ncatfac else raw
    point[length(setup$cont_flat) + seq_along(setup$cat_flat)] <- ext_bw * setup$bandwidth.scale.categorical
  }

  point
}

npRmpiNomadShadowSearchConditionalDistribution <- function(xdat,
                                                           ydat,
                                                           template,
                                                           setup,
                                                           reg.args,
                                                           opt.args,
                                                           degree.search,
                                                           x0,
                                                           bbin,
                                                           lb,
                                                           ub,
                                                           nomad.nmulti = 1L,
                                                           nomad.inner.nmulti = 0L,
                                                           random.seed = 42L,
                                                           use.runtime.bandwidth.progress = FALSE) {
  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
  old.messages <- getOption("np.messages")
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)

  if (!isTRUE(rank == 0L))
    options(np.messages = FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  on.exit(options(np.messages = old.messages), add = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)

  set.seed(as.integer(random.seed))
  if (isTRUE(rank == 0L))
    .np_nomad_baseline_note(degree.search$start.degree)

  bwdim <- length(setup$cont_flat) + length(setup$cat_flat)
  ndeg <- length(degree.search$start.degree)
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[bwdim + seq_len(ndeg)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npcdistbw_nomad_point_to_bw(point[seq_len(bwdim)], template = template, setup = setup)

    eval.reg.args <- reg.args
    eval.reg.args$regtype <- "lp"
    eval.reg.args$pregtype <- "Local-Polynomial"
    eval.reg.args$degree <- degree
    eval.reg.args$bernstein.basis <- degree.search$bernstein.basis
    eval.reg.args$regtype.engine <- "lp"
    eval.reg.args$degree.engine <- degree
    eval.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

    tbw <- .npcdistbw_build_condbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = eval.reg.args
    )

    out <- .npcdistbw_eval_only(
      xdat = xdat,
      ydat = ydat,
      gydat = opt.args$gydat,
      bws = tbw,
      do.full.integral = if (is.null(opt.args$do.full.integral)) FALSE else opt.args$do.full.integral,
      ngrid = if (is.null(opt.args$ngrid)) 100L else opt.args$ngrid,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier,
      force.local = FALSE
    )
    nomad.num.feval.total <<- nomad.num.feval.total + as.numeric(out$num.feval[1L])
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])

    list(
      objective = out$objective,
      degree = degree,
      num.feval = out$num.feval,
      num.feval.fast = out$num.feval.fast
    )
  }

  external.progress <- if (isTRUE(use.runtime.bandwidth.progress) &&
                           isTRUE(rank == 0L) &&
                           !is.null(.np_progress_runtime$bandwidth_state)) {
    .np_progress_runtime$bandwidth_state
  } else {
    NULL
  }
  search.engine.used <- if (identical(degree.search$engine, "nomad+powell")) {
    "nomad"
  } else {
    degree.search$engine
  }

  search.result <- .np_nomad_search(
    engine = search.engine.used,
    baseline_record = NULL,
    start_degree = degree.search$start.degree,
    x0 = x0,
    bbin = bbin,
    lb = lb,
    ub = ub,
    eval_fun = eval_fun,
    build_payload = function(point, best_record, solution, interrupted) {
      list(
        payload = NULL,
        objective = as.numeric(best_record$objective),
        powell.time = NA_real_
      )
    },
    direction = "min",
    objective_name = "fval",
    nmulti = nomad.nmulti,
    nomad.inner.nmulti = nomad.inner.nmulti,
    random.seed = random.seed,
    progress_state = external.progress,
    manage_progress_lifecycle = is.null(external.progress),
    bind_bandwidth_runtime = !is.null(external.progress),
    handoff_before_build = identical(degree.search$engine, "nomad+powell"),
    nomad.opts = list(
      DIRECTION_TYPE = "ORTHO 2N",
      QUAD_MODEL_SEARCH = "no",
      NM_SEARCH = "no",
      SPECULATIVE_SEARCH = "no",
      EVAL_OPPORTUNISTIC = "no"
    ),
    degree_spec = list(
      initial = degree.search$start.degree,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      user_supplied = degree.search$start.user
    )
  )
  if (!identical(search.engine.used, degree.search$engine))
    search.result$method <- degree.search$engine

  search.result$best_payload <- NULL
  search.result$powell.time <- NA_real_
  search.result$num.feval.total <- as.numeric(nomad.num.feval.total)
  search.result$num.feval.fast.total <- as.numeric(nomad.num.feval.fast.total)
  search.result
}

.npcdistbw_nomad_search <- function(xdat,
                                    ydat,
                                    bws,
                                    reg.args,
                                    opt.args,
                                    degree.search,
                                    nomad.inner.nmulti = 0L,
                                    random.seed = 42L) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$pregtype <- "Local-Polynomial"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template.reg.args$regtype.engine <- "lp"
  template.reg.args$degree.engine <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

  template <- .npcdistbw_build_condbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args
  )

  if (!identical(template$type, "fixed"))
    stop("automatic degree search with search.engine='nomad' currently requires bwtype='fixed'")
  setup <- .npcdistbw_nomad_bw_setup(xdat = xdat, ydat = ydat, template = template)
  bwdim <- length(setup$cont_flat) + length(setup$cat_flat)
  ndeg <- length(degree.search$start.degree)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(dim(ydat)[2]+dim(xdat)[2]) else npValidateNmulti(opt.args$nmulti[1L])
  cont_lower <- npGetScaleFactorSearchLower(
    template,
    argname = "template$scale.factor.search.lower"
  )
  bw_lower <- c(rep.int(cont_lower, length(setup$cont_flat)), rep.int(0, length(setup$cat_flat)))
  bw_upper <- c(rep.int(1e6, length(setup$cont_flat)), setup$cat_upper * setup$bandwidth.scale.categorical)

  x0 <- c(
    .np_nomad_complete_start_point(
      point = {
        raw <- c(template$ybw, template$xbw)
        if (all(raw == 0)) NULL else .npcdistbw_nomad_bw_to_point(raw, template = template, setup = setup)
      },
      lower = bw_lower,
      upper = bw_upper,
      ncont = length(setup$cont_flat)
    ),
    as.integer(degree.search$start.degree)
  )
  lb <- c(bw_lower, degree.search$lower)
  ub <- c(bw_upper, degree.search$upper)
  bbin <- c(rep.int(0L, bwdim), rep.int(1L, ndeg))
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[bwdim + seq_len(ndeg)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npcdistbw_nomad_point_to_bw(point[seq_len(bwdim)], template = template, setup = setup)

    eval.reg.args <- reg.args
    eval.reg.args$regtype <- "lp"
    eval.reg.args$pregtype <- "Local-Polynomial"
    eval.reg.args$degree <- degree
    eval.reg.args$bernstein.basis <- degree.search$bernstein.basis
    eval.reg.args$regtype.engine <- "lp"
    eval.reg.args$degree.engine <- degree
    eval.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

    tbw <- .npcdistbw_build_condbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = eval.reg.args
    )

    out <- .npcdistbw_eval_only(
      xdat = xdat,
      ydat = ydat,
      gydat = opt.args$gydat,
      bws = tbw,
      do.full.integral = if (is.null(opt.args$do.full.integral)) FALSE else opt.args$do.full.integral,
      ngrid = if (is.null(opt.args$ngrid)) 100L else opt.args$ngrid,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier
    )
    nomad.num.feval.total <<- nomad.num.feval.total + as.numeric(out$num.feval[1L])
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])

    list(
      objective = out$objective,
      degree = degree,
      num.feval = out$num.feval,
      num.feval.fast = out$num.feval.fast
    )
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    degree <- as.integer(best_record$degree)
    bw_vec <- .npcdistbw_nomad_point_to_bw(point[seq_len(bwdim)], template = template, setup = setup)
    powell.elapsed <- NA_real_

    build_direct_payload <- function() {
      final.reg.args <- reg.args
      final.reg.args$regtype <- "lp"
      final.reg.args$pregtype <- "Local-Polynomial"
      final.reg.args$degree <- degree
      final.reg.args$bernstein.basis <- degree.search$bernstein.basis
      final.reg.args$regtype.engine <- "lp"
      final.reg.args$degree.engine <- degree
      final.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

      tbw <- .npcdistbw_build_condbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bw_vec,
        bandwidth.compute = FALSE,
        reg.args = final.reg.args
      )
      tbw$fval <- as.numeric(best_record$objective)
      tbw$ifval <- as.numeric(best_record$objective)
      tbw$num.feval <- as.numeric(nomad.num.feval.total)
      tbw$num.feval.fast <- as.numeric(nomad.num.feval.fast.total)
      tbw$fval.history <- as.numeric(best_record$objective)
      tbw$eval.history <- if (isTRUE(nomad.num.feval.total > 0)) rep(1, max(1L, as.integer(nomad.num.feval.total))) else 1
      tbw$invalid.history <- 0
      tbw$timing <- NA_real_
      tbw$total.time <- NA_real_

      payload <- npcdistbw.condbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = tbw,
        bandwidth.compute = FALSE
      )
      if (!is.null(payload$method) && length(payload$method))
        payload$pmethod <- bwmToPrint(as.character(payload$method[1L]))
      payload
    }

    direct.payload <- build_direct_payload()
    if (is.null(direct.payload$timing.profile) && is.list(best_record$timing.profile))
      direct.payload$timing.profile <- best_record$timing.profile
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$pregtype <- "Local-Polynomial"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.reg.args$regtype.engine <- "lp"
      hot.reg.args$degree.engine <- degree
      hot.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
      hot.opt.args <- opt.args
      hot.opt.args$nmulti <- .np_nomad_powell_hotstart_nmulti("disable_multistart")
      powell.start <- proc.time()[3L]
      hot.payload <- .np_nomad_with_powell_progress(
        degree,
        .npcdistbw_run_fixed_degree_collective(
          xdat = xdat,
          ydat = ydat,
          bws = bw_vec,
          reg.args = hot.reg.args,
          opt.args = hot.opt.args
        )
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
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  if (.npRmpi_has_active_slave_pool(comm = 1L) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    search.template <- list(
      scaling = isTRUE(template$scaling),
      ybw = template$ybw,
      xbw = template$xbw
    )
    search.setup <- list(
      cont_flat = setup$cont_flat,
      cont_scale = setup$cont_scale,
      cat_flat = setup$cat_flat,
      ncatfac = setup$ncatfac,
      bandwidth.scale.categorical = setup$bandwidth.scale.categorical
    )
    search.degree <- list(
      engine = degree.search$engine,
      start.degree = degree.search$start.degree,
      candidates = degree.search$candidates,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      start.user = degree.search$start.user,
      bernstein.basis = degree.search$bernstein.basis
    )
    search.opt.args <- list(
      gydat = opt.args$gydat,
      do.full.integral = opt.args$do.full.integral,
      ngrid = opt.args$ngrid,
      penalty.multiplier = opt.args$penalty.multiplier
    )

    mc <- substitute(
      get("npRmpiNomadShadowSearchConditionalDistribution", envir = asNamespace("npRmpi"), inherits = FALSE)(
        XDAT,
        YDAT,
        TEMPLATE,
        SETUP,
        REGARGS,
        OPTARGS,
        DEGREESEARCH,
        X0,
        BBIN,
        LB,
        UB,
        NOMADNMULTI,
        INNERNMULTI,
        RSEED,
        RPROGRESS
      ),
      list(
        XDAT = xdat,
        YDAT = ydat,
        TEMPLATE = search.template,
        SETUP = search.setup,
        REGARGS = reg.args,
        OPTARGS = search.opt.args,
        DEGREESEARCH = search.degree,
        X0 = x0,
        BBIN = bbin,
        LB = lb,
        UB = ub,
        NOMADNMULTI = nomad.nmulti,
        INNERNMULTI = nomad.inner.nmulti,
        RSEED = random.seed,
        RPROGRESS = TRUE
      )
    )

    if (isTRUE(.npRmpi_autodispatch_called_from_bcast())) {
      search.result <- eval(mc, envir = environment())
    } else {
      search.result <- .npRmpi_bcast_cmd_expr(mc, comm = 1L, caller.execute = TRUE)
    }
    if (!is.null(search.result$num.feval.total))
      nomad.num.feval.total <- as.numeric(search.result$num.feval.total[1L])
    if (!is.null(search.result$num.feval.fast.total))
      nomad.num.feval.fast.total <- as.numeric(search.result$num.feval.fast.total[1L])
    best.solution <- NULL
    if (!is.null(search.result$best.restart) &&
        is.finite(search.result$best.restart) &&
        length(search.result$restart.results) >= as.integer(search.result$best.restart)) {
      best.solution <- search.result$restart.results[[as.integer(search.result$best.restart)]]
    }

    payload_result <- build_payload(
      point = search.result$best_point,
      best_record = search.result$best,
      solution = best.solution,
      interrupted = !isTRUE(search.result$completed)
    )

    if (is.list(payload_result) && !is.null(payload_result$payload)) {
      search.result$best_payload <- payload_result$payload
      if (!is.null(payload_result$powell.time))
        search.result$powell.time <- as.numeric(payload_result$powell.time[1L])
      if (!is.null(payload_result$objective) &&
          .np_degree_better(payload_result$objective, search.result$best$objective, direction = "min")) {
        search.result$best$objective <- as.numeric(payload_result$objective[1L])
      }
    } else {
      search.result$best_payload <- payload_result
    }

    return(.npRmpi_reconcile_nomad_search_timing(search.result))
  }

  .np_nomad_baseline_note(degree.search$start.degree)

  .np_nomad_search(
    engine = degree.search$engine,
    baseline_record = baseline.record,
    start_degree = degree.search$start.degree,
    x0 = x0,
    bbin = bbin,
    lb = lb,
    ub = ub,
    eval_fun = eval_fun,
    build_payload = build_payload,
    direction = "min",
    objective_name = "fval",
    nmulti = nomad.nmulti,
    nomad.inner.nmulti = nomad.inner.nmulti,
    random.seed = random.seed,
    nomad.opts = list(
      DIRECTION_TYPE = "ORTHO 2N",
      QUAD_MODEL_SEARCH = "no",
      NM_SEARCH = "no",
      SPECULATIVE_SEARCH = "no",
      EVAL_OPPORTUNISTIC = "no"
    ),
    degree_spec = list(
      initial = degree.search$start.degree,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      user_supplied = degree.search$start.user
    )
  )
}

.npcdistbw_degree_search_controls <- function(regtype,
                                              regtype.named,
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
                                              bernstein.named) {
  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  if (identical(degree.select, "manual"))
    return(NULL)
  search.engine <- .np_degree_search_engine_controls(search.engine)

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (ncon < 1L)
    stop("automatic degree search requires at least one continuous conditioning predictor")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = ncon,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  baseline.degree <- rep.int(0L, ncon)
  # Density/distribution degree search should anchor on the local-constant
  # baseline unless the user explicitly requests another start.
  default.start.degree <- baseline.degree
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, default.start.degree))
  } else {
    start.raw <- npValidateGlpDegree(regtype = "lp", degree = degree.start, ncon = ncon, argname = "degree.start")
    out.of.range <- vapply(seq_len(ncon), function(j) !(start.raw[j] %in% bounds$candidates[[j]]), logical(1))
    if (any(out.of.range))
      stop("degree.start must lie within the searched degree candidates for every continuous conditioning predictor")
    start.raw
  }

  list(
    method = if (identical(search.engine, "cell")) degree.select else search.engine,
    engine = search.engine,
    candidates = bounds$candidates,
    lower = bounds$lower,
    upper = bounds$upper,
    baseline.degree = baseline.degree,
    start.degree = start.degree,
    start.user = !is.null(degree.start),
    basis = if (missing(basis) || is.null(basis)) "glp" else as.character(basis[1L]),
    nobs = as.integer(nobs[1L]),
    restarts = npValidateNonNegativeInteger(degree.restarts, "degree.restarts"),
    max.cycles = npValidatePositiveInteger(degree.max.cycles, "degree.max.cycles"),
    verify = npValidateScalarLogical(degree.verify, "degree.verify"),
    bernstein.basis = bern.auto
  )
}

.npcdistbw_attach_degree_search <- function(bws, search_result) {
  metadata <- list(
    mode = search_result$method,
    verify = isTRUE(search_result$verify),
    completed = isTRUE(search_result$completed),
    certified = isTRUE(search_result$certified),
    interrupted = isTRUE(search_result$interrupted),
    baseline.degree = search_result$baseline$degree,
    baseline.fval = search_result$baseline$objective,
    best.degree = search_result$best$degree,
    best.fval = search_result$best$objective,
    nomad.time = search_result$nomad.time,
    powell.time = search_result$powell.time,
    optim.time = search_result$optim.time,
    n.unique = search_result$n.unique,
    n.visits = search_result$n.visits,
    n.cached = search_result$n.cached,
    grid.size = search_result$grid.size,
    best.restart = search_result$best.restart,
    restart.starts = search_result$restart.starts,
    restart.degree.starts = search_result$restart.degree.starts,
    restart.bandwidth.starts = search_result$restart.bandwidth.starts,
    restart.start.info = search_result$restart.start.info,
    restart.results = search_result$restart.results,
    trace = search_result$trace
  )

  if (!is.null(search_result$nomad.time))
    bws$nomad.time <- as.numeric(search_result$nomad.time[1L])
  if (!is.null(search_result$powell.time))
    bws$powell.time <- as.numeric(search_result$powell.time[1L])
  if (!is.null(search_result$optim.time) && is.finite(search_result$optim.time))
    bws$total.time <- as.numeric(search_result$optim.time[1L])
  if (is.null(bws$timing.profile) && is.list(search_result$best$timing.profile))
    bws$timing.profile <- search_result$best$timing.profile
  bws <- .np_attach_nomad_restart_summary(bws, search_result)
  bws$degree.search <- metadata
  bws
}

npcdistbw.NULL <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           bws, ...){
    dots <- list(...)
    if (.npRmpi_npcdistbw_bounded_adaptive_requested(
      bwtype = dots$bwtype,
      cxkerbound = if (is.null(dots$cxkerbound)) "none" else dots$cxkerbound,
      cykerbound = if (is.null(dots$cykerbound)) "none" else dots$cykerbound,
      bandwidth.compute = if (is.null(dots$bandwidth.compute)) TRUE else dots$bandwidth.compute
    ))
      stop("bounded adaptive_nn remains unsupported for npcdistbw() in npRmpi",
           call. = FALSE)
    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain y names and 'toFrame'
    ydat <- toFrame(ydat)

    ## do bandwidths

    bws = double(ncol(ydat)+ncol(xdat))

    tbw <- npcdistbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw
  }

npcdistbw.default <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           gydat,
           bws,
           bandwidth.compute = TRUE,
           bwmethod,
           bwscaling,
           bwtype,
           cfac.dir,
           scale.factor.init,
           cxkerbound,
           cxkerlb,
           cxkerorder,
           cxkertype,
           cxkerub,
           cykerbound,
           cykerlb,
           cykerorder,
           cykertype,
           cykerub,
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
           oxkertype,
           oykertype,
           penalty.multiplier,
           remin,
           scale.init.categorical.sample,
           scale.factor.search.lower = NULL,
           small,
           tol,
           transform.bounds,
           uxkertype,
           regtype = c("lc", "ll", "lp"),
           basis = c("glp", "additive", "tensor"),
           degree = NULL,
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
           bernstein.basis = FALSE,
           ## dummy arguments for condbandwidth() function call
           ...){
    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain y names and 'toFrame'
    ydat <- toFrame(ydat)

    x.info <- untangle(xdat)
    y.info <- untangle(ydat)

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
      where = "npcdistbw"
    )

    if (isTRUE(nomad.shortcut$enabled)) {
      if ("degree" %in% mc.names)
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$bwtype, c("fixed", "generalized_nn", "adaptive_nn")))[1L], "fixed"))
        stop("nomad=TRUE currently requires bwtype='fixed'")
      if ("degree.select" %in% mc.names &&
          identical(as.character(match.arg(nomad.shortcut$values$degree.select, c("manual", "coordinate", "exhaustive")))[1L], "manual"))
        stop("nomad=TRUE requires automatic degree search; use degree.select='coordinate' or 'exhaustive'")
      if ("search.engine" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$search.engine, c("nomad+powell", "cell", "nomad")))[1L] %in%
              c("nomad", "nomad+powell")))
        stop("nomad=TRUE requires search.engine='nomad' or 'nomad+powell'")
      if ("bernstein.basis" %in% mc.names &&
          !isTRUE(npValidateGlpBernstein(regtype = "lp",
                                        bernstein.basis = nomad.shortcut$values$bernstein.basis)))
        stop("nomad=TRUE currently requires bernstein.basis=TRUE")
      if ("degree.verify" %in% mc.names &&
          isTRUE(npValidateScalarLogical(nomad.shortcut$values$degree.verify, "degree.verify")))
        stop("nomad=TRUE currently requires degree.verify=FALSE")
    }

    regtype.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "regtype")
    basis.named <- any(mc.names == "basis")
    degree.named <- any(mc.names == "degree")
    bernstein.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "bernstein.basis")

    regtype <- if (!is.null(nomad.shortcut$values$regtype)) {
      match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp"))
    } else {
      "lc"
    }
    if (identical(regtype, "lc") && (basis.named || degree.named || bernstein.named))
      stop("regtype='lc' does not accept basis/degree/bernstein.basis; use regtype='lp' for local-polynomial controls")
    if (identical(regtype, "ll")) {
      if (degree.named)
        stop("regtype='ll' uses canonical LP(degree=1, basis='glp'); remove 'degree' or use regtype='lp'")
      if (basis.named && !identical(match.arg(basis), "glp"))
        stop("regtype='ll' uses canonical basis='glp'; use regtype='lp' for alternate LP bases")
      if (bernstein.named && isTRUE(bernstein.basis))
        stop("regtype='ll' uses canonical bernstein.basis=FALSE; use regtype='lp' for Bernstein LP")
    }

    bernstein.value <- if (!is.null(nomad.shortcut$values$bernstein.basis)) {
      nomad.shortcut$values$bernstein.basis
    } else {
      bernstein.basis
    }
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    degree.setup <- npSetupGlpDegree(
      regtype = regtype,
      degree = degree,
      ncon = sum(x.info$icon),
      degree.select = degree.select.value
    )
    spec <- npCanonicalConditionalRegSpec(
      regtype = regtype,
      basis = basis,
      degree = degree.setup,
      bernstein.basis = bernstein.value,
      ncon = sum(x.info$icon),
      where = "npcdistbw"
    )
    pregtype <- switch(spec$regtype,
                       lc = "Local-Constant",
                       ll = "Local-Linear",
                       lp = "Local-Polynomial")

    search.mc.names <- names(mc)
    lp.dot.args <- list(...)
    random.seed.value <- .np_degree_extract_random_seed(lp.dot.args)
    search.engine.value <- if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell"
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)
    degree.min.value <- nomad.shortcut$values$degree.min
    degree.max.value <- nomad.shortcut$values$degree.max
    degree.start.value <- if ("degree.start" %in% search.mc.names) degree.start else NULL
    degree.restarts.value <- if ("degree.restarts" %in% search.mc.names) degree.restarts else 0L
    degree.max.cycles.value <- if ("degree.max.cycles" %in% search.mc.names) degree.max.cycles else 20L
    degree.verify.value <- if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE
    degree.search <- .npcdistbw_degree_search_controls(
      regtype = regtype,
      regtype.named = regtype.named,
      ncon = sum(x.info$icon),
      nobs = NROW(xdat),
      basis = if (basis.named) basis else "glp",
      degree.select = degree.select.value,
      search.engine = search.engine.value,
      degree.min = degree.min.value,
      degree.max = degree.max.value,
      degree.start = degree.start.value,
      degree.restarts = degree.restarts.value,
      degree.max.cycles = degree.max.cycles.value,
      degree.verify = degree.verify.value,
      bernstein.basis = bernstein.value,
      bernstein.named = bernstein.named
    )
    nomad.inner <- .np_nomad_validate_inner_multistart(
      call_names = search.mc.names,
      dot.args = lp.dot.args,
      nomad.nmulti = nomad.nmulti,
      regtype = regtype,
      automatic.degree.search = !is.null(degree.search),
      search.engine = if (is.null(degree.search)) "" else degree.search$engine
    )
    nomad.inner.named <- nomad.inner$named
    nomad.inner.nmulti <- nomad.inner$nmulti

    if (!is.null(degree.search)) {
      spec$bernstein.basis <- degree.search$bernstein.basis
      spec$bernstein.basis.engine <- degree.search$bernstein.basis
    }

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    mc.names <- names(mc)
    margs <- c("bwmethod", "bwscaling", "bwtype", "cxkertype", "cxkerorder",
               "cxkerbound", "cxkerlb", "cxkerub",
               "cykertype", "cykerorder", "cykerbound", "cykerlb", "cykerub",
               "uxkertype", "oxkertype", "oykertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    y.idx <- seq_len(length(ydat))
    x.idx <- seq_len(length(xdat))
    bw.args <- list(
      xbw = bws[length(ydat) + x.idx],
      ybw = bws[y.idx],
      uykertype = "aitchisonaitken",
      nobs = nrow(xdat),
      xdati = x.info,
      ydati = y.info,
      xnames = names(xdat),
      ynames = names(ydat),
      bandwidth.compute = bandwidth.compute,
      regtype = spec$regtype,
      pregtype = pregtype,
      basis = spec$basis,
      degree = spec$degree,
      bernstein.basis = spec$bernstein.basis,
      regtype.engine = spec$regtype.engine,
      basis.engine = spec$basis.engine,
      degree.engine = spec$degree.engine,
      bernstein.basis.engine = spec$bernstein.basis.engine
    )
    if (any.m) {
      nms <- mc.names[m]
      bw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    reg.args <- bw.args[setdiff(names(bw.args), c("xbw", "ybw", "nobs", "xdati", "ydati", "xnames", "ynames", "bandwidth.compute"))]
    tbw <- do.call(condbandwidth, bw.args)
    .npRmpi_require_active_slave_pool(where = "npcdistbw()")
    use.local.compiled.adaptive.cvls <- bandwidth.compute &&
      identical(tbw$method, "cv.ls") &&
      identical(tbw$type, "adaptive_nn")
    keep_local_cvls_nn <- bandwidth.compute &&
      identical(tbw$method, "cv.ls") &&
      (use.local.compiled.adaptive.cvls ||
       (identical(tbw$regtype.engine, "lp") &&
        identical(tbw$type, "generalized_nn")))
    keep_local_raw_degree1_cvls <- bandwidth.compute &&
      identical(tbw$method, "cv.ls") &&
      identical(tbw$type, "fixed") &&
      npIsRawDegreeOneConditionalSpec(spec, tbw$xncon)
    if (.npRmpi_autodispatch_active() &&
        !keep_local_cvls_nn &&
        !keep_local_raw_degree1_cvls &&
        is.null(degree.search))
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    ## next grab dummies for actual bandwidth selection and perform call

    mc.names <- names(mc)
    margs <- c("gydat", "nmulti", "remin", "itmax", "do.full.integral", "ngrid", "ftol",
               "tol", "small", "memfac",
               "lbc.dir", "dfc.dir", "cfac.dir","initc.dir",
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.factor.search.lower",
               "scale.init.categorical.sample",
               "transform.bounds",
               "invalid.penalty",
               "penalty.multiplier")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (any.m) {
      nms <- mc.names[m]
      opt.args <- mget(nms, envir = environment(), inherits = FALSE)
    } else {
      opt.args <- list()
    }
    opt.args <- c(list(bandwidth.compute = bandwidth.compute), opt.args)
    reg.args$scale.factor.search.lower <- scale.factor.search.lower
    opt.args$scale.factor.search.lower <- scale.factor.search.lower

    if (!is.null(degree.search)) {
      if (identical(degree.search$engine, "cell")) {
        eval_fun <- function(degree.vec) {
          cell.reg.args <- reg.args
          cell.reg.args$regtype <- "lp"
          cell.reg.args$pregtype <- "Local-Polynomial"
          cell.reg.args$degree <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
          cell.reg.args$regtype.engine <- "lp"
          cell.reg.args$degree.engine <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
          cell.bws <- .npcdistbw_run_fixed_degree(
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            reg.args = cell.reg.args,
            opt.args = opt.args
          )
          list(
            objective = as.numeric(cell.bws$fval[1L]),
            payload = cell.bws,
            num.feval = if (!is.null(cell.bws$num.feval)) as.numeric(cell.bws$num.feval[1L]) else NA_real_
          )
        }

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
          objective_name = "fval"
        )
      } else {
        search.result <- .npcdistbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          reg.args = reg.args,
          opt.args = opt.args,
          degree.search = degree.search,
          nomad.inner.nmulti = nomad.inner.nmulti,
          random.seed = random.seed.value
        )
      }
      tbw <- .npcdistbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
    } else {
      tbw <- .npcdistbw_build_condbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bws,
        bandwidth.compute = bandwidth.compute,
        reg.args = reg.args
      )
      bwsel.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting conditional distribution bandwidth",
        do.call(npcdistbw.condbandwidth, bwsel.args)
      )
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)
  }
