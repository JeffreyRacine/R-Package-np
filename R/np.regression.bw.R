npregbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    dots <- mc$...

    if (length(dots) == 0L)
      stop("invoked without arguments")

    dot.names <- names(dots)

    first.val <- eval(dots[[1L]], envir = parent.frame())
    if (inherits(first.val, "formula"))
      return(UseMethod("npregbw", first.val))

    if (!is.null(dot.names) && any(dot.names == "formula")) {
      formula.val <- eval(dots[[which(dot.names == "formula")[1L]]], envir = parent.frame())
      return(UseMethod("npregbw", formula.val))
    }

    if (!is.null(dot.names) && any(dot.names == "bws")) {
      bws.val <- eval(dots[[which(dot.names == "bws")[1L]]], envir = parent.frame())
      return(UseMethod("npregbw", bws.val))
    }

    UseMethod("npregbw", first.val)
  }

npregbw.formula <-
  function(formula, data, subset, na.action, call, ...){

    orig.ts <- if (missing(data))
      sapply(eval(attr(terms(formula), "variables"), environment(formula)), inherits, "ts")
    else sapply(eval(attr(terms(formula), "variables"), data, environment(formula)), inherits, "ts")

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    if(all(orig.ts)){
      args <- (as.list(attr(terms(formula), "variables"))[-1])
      formula <- terms(formula)
      attr(formula, "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
      mf[["formula"]] <- formula
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(terms(formula), "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      formula <- terms(formula)
      attr(formula, "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
      mf[["formula"]] <- formula
    }
      
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]
    
    dots <- list(...)
    tbw <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat), dots))

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)
    
    tbw
  }

npregbw.NULL <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npregbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_as_generic_call("npregbw", match.call()),
        parent.frame()))

    xdat <- toFrame(xdat)

    bws = double(dim(xdat)[2])
    
    tbw <- npregbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)
    
    tbw
  }

npregbw.rbandwidth <- 
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws, bandwidth.compute = TRUE,
           nmulti, remin = TRUE, itmax = 10000,
           ftol = 1.490116e-07, tol = 1.490116e-04, small = 1.490116e-05,
           lbc.dir = 0.5, dfc.dir = 3, cfac.dir = 2.5*(3.0-sqrt(5)),initc.dir = 1.0, 
           lbd.dir = 0.1, hbd.dir = 1, dfac.dir = 0.25*(3.0-sqrt(5)), initd.dir = 1.0, 
           lbc.init = 0.1, hbc.init = 2.0, cfac.init = 0.5, 
           lbd.init = 0.1, hbd.init = 0.9, dfac.init = 0.375, 
          scale.init.categorical.sample = FALSE,
          transform.bounds = FALSE,
          invalid.penalty = c("baseline","dbmax"),
          penalty.multiplier = 10,
          ...){
    elapsed.start <- proc.time()[3]
    .npRmpi_require_active_slave_pool(where = "npregbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_as_generic_call("npregbw", match.call()),
        parent.frame()))

    xdat <- toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- min(5,dim(xdat)[2])
    }

    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (length(bws$bw) != dim(xdat)[2])
      stop("length of bandwidth vector does not match number of columns of 'xdat'")

    ccon = unlist(lapply(xdat[,bws$icon, drop = FALSE],class))
    if ((any(bws$icon) && !all((ccon == "integer") | (ccon == "numeric"))) ||
        (any(bws$iord) && !all(sapply(xdat[,bws$iord, drop = FALSE],inherits, "ordered"))) ||
        (any(bws$iuno) && !all(sapply(xdat[,bws$iuno, drop = FALSE],inherits, "factor"))))
      stop("supplied bandwidths do not match 'xdat' in type")

    if (dim(xdat)[1] != length(ydat))
      stop("number of regression data and response data do not match")

    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]
    
    nrow = dim(xdat)[1]
    ncol = dim(xdat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat = toMatrix(xdat)

    runo = xdat[, bws$iuno, drop = FALSE]
    rcon = xdat[, bws$icon, drop = FALSE]
    rord = xdat[, bws$iord, drop = FALSE]

    tbw <- bws
    tbw$basis <- npValidateLpBasis(regtype = tbw$regtype,
                                   basis = tbw$basis)
    tbw$degree <- npValidateGlpDegree(regtype = tbw$regtype,
                                          degree = tbw$degree,
                                          ncon = tbw$ncon)
    tbw$bernstein.basis <- npValidateGlpBernstein(regtype = tbw$regtype,
                                                bernstein.basis = tbw$bernstein.basis)

    mysd <- EssDee(rcon)
    nconfac <- nrow^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- nrow^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    invalid.penalty <- match.arg(invalid.penalty)
    penalty_mode <- ifelse(invalid.penalty == "baseline", 1L, 0L)

    reg.c <- npRegtypeToC(regtype = tbw$regtype,
                          degree = tbw$degree,
                          ncon = tbw$ncon,
                          context = "npregbw")
    npCheckRegressionDesignCondition(reg.code = reg.c$code,
                                     xcon = rcon,
                                     basis = tbw$basis,
                                     degree = tbw$degree,
                                     bernstein.basis = tbw$bernstein.basis,
                                     where = "npregbw")
    if (identical(tbw$regtype, "lp") &&
        identical(reg.c$code, REGTYPE_GLP) &&
        !isTRUE(transform.bounds))
      transform.bounds <- TRUE
    degree.c <- if (tbw$ncon > 0) {
      as.integer(if (is.null(reg.c$degree)) rep.int(0L, tbw$ncon) else reg.c$degree)
    } else {
      integer(1)
    }

    if (bandwidth.compute){
      myopti = list(num_obs_train = dim(xdat)[1], 
        iMultistart = ifelse(nmulti==0,IMULTI_FALSE,IMULTI_TRUE),
        iNum_Multistart = nmulti,
        int_use_starting_values = ifelse(all(bws$bw==0),USE_START_NO, USE_START_YES),
        int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
        BANDWIDTH_reg_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=ifelse(remin,RE_MIN_TRUE,RE_MIN_FALSE), 
        int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
        bwmethod = switch(bws$method,
          cv.aic = BWM_CVAIC,
          cv.ls = BWM_CVLS),
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
        nuno = bws$nuno,
        nord = bws$nord,
        ncon = bws$ncon,
        regtype = reg.c$code,
        int_do_tree = ifelse(options('np.tree'), DO_TREE_YES, DO_TREE_NO),
        scale.init.categorical.sample = scale.init.categorical.sample,
        dfc.dir = dfc.dir,
        transform.bounds = transform.bounds)
      
      myoptd = list(ftol=ftol, tol=tol, small=small,
        lbc.dir = lbc.dir, cfac.dir = cfac.dir, initc.dir = initc.dir, 
        lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir, 
        lbc.init = lbc.init, hbc.init = hbc.init, cfac.init = cfac.init, 
        lbd.init = lbd.init, hbd.init = hbd.init, dfac.init = dfac.init, 
        nconfac = nconfac, ncatfac = ncatfac)

        cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

        system.time(myout <-
          .Call("C_np_regression_bw",
                as.double(runo), as.double(rord), as.double(rcon), as.double(ydat),
                as.double(mysd),
                as.integer(myopti), as.double(myoptd),
                as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
                as.integer(max(1, nmulti)),
                as.integer(penalty_mode),
                as.double(penalty.multiplier),
                as.integer(degree.c),
                as.integer(isTRUE(tbw$bernstein.basis)),
                as.integer(npLpBasisCode(tbw$basis)),
                as.double(cker.bounds.c$lb),
                as.double(cker.bounds.c$ub),
                PACKAGE = "npRmpi"))[1]
      

      rorder = numeric(ncol)
      rorder[c((1:ncol)[bws$icon], (1:ncol)[bws$iuno], (1:ncol)[bws$iord])]=1:ncol

      tbw$bw <- myout$bw[rorder]
      tbw$fval <- myout$fval[1]
      tbw$ifval <- myout$fval[2]
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
      tbw$num.feval.fallback <- myout$fallback.history[1]
      tbw$fval.history <- myout$fval.history
      tbw$eval.history <- myout$eval.history
      tbw$invalid.history <- myout$invalid.history
      tbw$timing <- myout$timing
    }

    tbw$total.time <- proc.time()[3] - elapsed.start

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
    


    tbw <- rbandwidth(bw = tbw$bw,
                      regtype = tbw$regtype,
                      basis = tbw$basis,
                      degree = tbw$degree,
                      bernstein.basis = tbw$bernstein.basis,
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
                      num.feval.fallback = tbw$num.feval.fallback,
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
                      total.time = tbw$total.time)
    tbw
  }

npregbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws,
           bandwidth.compute = TRUE, nmulti,
           remin, itmax, ftol, tol, small,
           lbc.dir, dfc.dir, cfac.dir, initc.dir, 
           lbd.dir, hbd.dir, dfac.dir, initd.dir, 
           lbc.init, hbc.init, cfac.init, 
           lbd.init, hbd.init, dfac.init,
           scale.init.categorical.sample,
           transform.bounds = FALSE,
           invalid.penalty = c("baseline","dbmax"),
           penalty.multiplier = 10,
           ## dummy arguments for later passing into rbandwidth()
           regtype, basis, degree, bernstein.basis, bwmethod, bwscaling, bwtype,
           ckertype, ckerorder, ckerbound, ckerlb, ckerub, ukertype, okertype,
           ...){
    .npRmpi_require_active_slave_pool(where = "npregbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(
        .npRmpi_autodispatch_as_generic_call("npregbw", match.call()),
        parent.frame()))
    npRejectLegacyLpArgs(names(list(...)), where = "npregbw")

    xdat <- toFrame(xdat)

    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    rb.args <- list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      xnames = names(xdat),
      ynames = deparse(substitute(ydat)),
      bandwidth.compute = bandwidth.compute
    )

    if (!missing(regtype)) rb.args$regtype <- regtype
    if (!missing(basis)) rb.args$basis <- basis
    if (!missing(degree)) rb.args$degree <- degree
    if (!missing(bernstein.basis)) rb.args$bernstein.basis <- bernstein.basis
    if (!missing(bwmethod)) rb.args$bwmethod <- bwmethod
    if (!missing(bwscaling)) rb.args$bwscaling <- bwscaling
    if (!missing(bwtype)) rb.args$bwtype <- bwtype
    if (!missing(ckertype)) rb.args$ckertype <- ckertype
    if (!missing(ckerorder)) rb.args$ckerorder <- ckerorder
    if (!missing(ckerbound)) rb.args$ckerbound <- ckerbound
    if (!missing(ckerlb)) rb.args$ckerlb <- ckerlb
    if (!missing(ckerub)) rb.args$ckerub <- ckerub
    if (!missing(ukertype)) rb.args$ukertype <- ukertype
    if (!missing(okertype)) rb.args$okertype <- okertype

    tbw <- do.call(rbandwidth, rb.args)

    opt.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
    if (!missing(bandwidth.compute)) opt.args$bandwidth.compute <- bandwidth.compute
    if (!missing(nmulti)) opt.args$nmulti <- nmulti
    if (!missing(remin)) opt.args$remin <- remin
    if (!missing(itmax)) opt.args$itmax <- itmax
    if (!missing(ftol)) opt.args$ftol <- ftol
    if (!missing(tol)) opt.args$tol <- tol
    if (!missing(small)) opt.args$small <- small
    if (!missing(lbc.dir)) opt.args$lbc.dir <- lbc.dir
    if (!missing(dfc.dir)) opt.args$dfc.dir <- dfc.dir
    if (!missing(cfac.dir)) opt.args$cfac.dir <- cfac.dir
    if (!missing(initc.dir)) opt.args$initc.dir <- initc.dir
    if (!missing(lbd.dir)) opt.args$lbd.dir <- lbd.dir
    if (!missing(hbd.dir)) opt.args$hbd.dir <- hbd.dir
    if (!missing(dfac.dir)) opt.args$dfac.dir <- dfac.dir
    if (!missing(initd.dir)) opt.args$initd.dir <- initd.dir
    if (!missing(lbc.init)) opt.args$lbc.init <- lbc.init
    if (!missing(hbc.init)) opt.args$hbc.init <- hbc.init
    if (!missing(cfac.init)) opt.args$cfac.init <- cfac.init
    if (!missing(lbd.init)) opt.args$lbd.init <- lbd.init
    if (!missing(hbd.init)) opt.args$hbd.init <- hbd.init
    if (!missing(dfac.init)) opt.args$dfac.init <- dfac.init
    if (!missing(scale.init.categorical.sample))
      opt.args$scale.init.categorical.sample <- scale.init.categorical.sample
    if (!missing(transform.bounds)) opt.args$transform.bounds <- transform.bounds
    if (!missing(invalid.penalty)) opt.args$invalid.penalty <- invalid.penalty
    if (!missing(penalty.multiplier)) opt.args$penalty.multiplier <- penalty.multiplier

    tbw <- do.call(npregbw.rbandwidth, opt.args)

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
    
  }
