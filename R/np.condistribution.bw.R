npcdistbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
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

npcdistbw.condbandwidth <- 
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           gydat = NULL,
           bws,
           bandwidth.compute = TRUE,
           cfac.dir = 2.5*(3.0-sqrt(5)),
           cfac.init = 0.5,
           dfac.dir = 0.25*(3.0-sqrt(5)),
           dfac.init = 0.375,
           dfc.dir = 3,
           do.full.integral = FALSE,
           ftol = 1.490116e-07,
           hbc.init = 2.0,
           hbd.dir = 1,
           hbd.init = 0.9,
           initc.dir = 1.0,
           initd.dir = 1.0,
           invalid.penalty = c("baseline","dbmax"),
           itmax = 10000,
           lbc.dir = 0.5,
           lbc.init = 0.1,
           lbd.dir = 0.1,
           lbd.init = 0.1,
           memfac = 500.0,
           ngrid = 100,
           nmulti,
           penalty.multiplier = 10,
           remin = TRUE,
           scale.init.categorical.sample = FALSE,
           small = 1.490116e-05,
           tol = 1.490116e-04,
           transform.bounds = FALSE,
           ...){
    ydat = toFrame(ydat)
    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- min(5,(dim(ydat)[2]+dim(xdat)[2]))
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
    nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
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
    if (bandwidth.compute &&
        identical(spec$regtype.engine, "lp") &&
        identical(bws$method, "cv.ls") &&
        !identical(bws$type, "fixed"))
      stop(
        "public npcdistbw() LP/LL cv.ls route is temporarily disabled pending low-memory shadow CV remediation"
      )

    .npRmpi_require_active_slave_pool(where = "npcdistbw()")
    if (.npRmpi_autodispatch_active())
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
      myopti = list(num_obs_train = nrow,
        num_obs_grid = nog,
        iMultistart = (if (nmulti==0) IMULTI_FALSE else IMULTI_TRUE),
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
        lbc.init = lbc.init, hbc.init = hbc.init, cfac.init = cfac.init, 
        lbd.init = lbd.init, hbd.init = hbd.init, dfac.init = dfac.init, 
        nconfac = nconfac, ncatfac = ncatfac)

      cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
      cyker.bounds.c <- npKernelBoundsMarshal(bws$cykerlb[bws$iycon], bws$cykerub[bws$iycon])

      if (bws$method != "normal-reference"){
        total.time <-
          system.time(myout <-
                        .Call("C_np_distribution_conditional_bw",
                              as.double(yuno), as.double(yord), as.double(ycon),
                              as.double(xuno), as.double(xord), as.double(xcon),
                              as.double(gyuno), as.double(gyord), as.double(gycon),
                              as.double(mysd),
                              as.integer(myopti), as.double(myoptd),
                              as.double(c(bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
                                          bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
                                          bws$xbw[bws$ixuno], bws$xbw[bws$ixord])),
                              as.integer(max(1, nmulti)),
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
                              PACKAGE="npRmpi"))[1]
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
        myout= list( bw = nbw, fval = c(NA,NA) )
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
           
    tbw
  }

npcdistbw.NULL <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           bws, ...){
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
           cfac.init,
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
           hbc.init,
           hbd.dir,
           hbd.init,
           initc.dir,
           initd.dir,
           invalid.penalty,
           itmax,
           lbc.dir,
           lbc.init,
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
           small,
           tol,
           transform.bounds,
           uxkertype,
           regtype = c("lc", "ll", "lp"),
           basis = c("glp", "additive", "tensor"),
           degree = NULL,
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
    regtype.named <- any(mc.names == "regtype")
    basis.named <- any(mc.names == "basis")
    degree.named <- any(mc.names == "degree")
    bernstein.named <- any(mc.names == "bernstein.basis")

    regtype <- if (regtype.named) match.arg(regtype) else "lc"
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

    spec <- npCanonicalConditionalRegSpec(
      regtype = regtype,
      basis = basis,
      degree = degree,
      bernstein.basis = bernstein.basis,
      ncon = sum(x.info$icon),
      where = "npcdistbw"
    )
    pregtype <- switch(spec$regtype,
                       lc = "Local-Constant",
                       ll = "Local-Linear",
                       lp = "Local-Polynomial")

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
    tbw <- do.call(condbandwidth, bw.args)
    if (bandwidth.compute &&
        identical(tbw$regtype.engine, "lp") &&
        identical(tbw$method, "cv.ls") &&
        !identical(tbw$type, "fixed"))
      stop(
        "public npcdistbw() LP/LL cv.ls route is temporarily disabled pending low-memory shadow CV remediation"
      )
    .npRmpi_require_active_slave_pool(where = "npcdistbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))
                        
    ## next grab dummies for actual bandwidth selection and perform call

    mc.names <- names(mc)
    margs <- c("gydat", "bandwidth.compute", "nmulti", "remin", "itmax", "do.full.integral", "ngrid", "ftol",
               "tol", "small", "memfac",
               "lbc.dir", "dfc.dir", "cfac.dir","initc.dir", 
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir", 
               "lbc.init", "hbc.init", "cfac.init", 
               "lbd.init", "hbd.init", "dfac.init", 
               "scale.init.categorical.sample",
               "transform.bounds",
               "invalid.penalty",
               "penalty.multiplier")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    bwsel.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
    if (any.m) {
      nms <- mc.names[m]
      bwsel.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    tbw <- do.call(npcdistbw.condbandwidth, bwsel.args)

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }
