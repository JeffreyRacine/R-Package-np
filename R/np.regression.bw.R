npregbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat"),
                                     eval_env = parent.frame())
    UseMethod("npregbw", target)
  }

npregbw.formula <-
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

    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]
    
    tbw <- npregbw(xdat = xdat, ydat = ydat, ...)

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
           bws,
           bandwidth.compute = TRUE,
           cfac.dir = 2.5*(3.0-sqrt(5)),
           cfac.init = 0.5,
           dfac.dir = 0.25*(3.0-sqrt(5)),
           dfac.init = 0.375,
           dfc.dir = 3,
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
           nmulti,
           penalty.multiplier = 10,
           remin = TRUE,
           scale.init.categorical.sample = FALSE,
           small = 1.490116e-05,
           tol = 1.490116e-04,
           transform.bounds = FALSE,
           ...){
    elapsed.start <- proc.time()[3]

    xdat <- toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- min(5,dim(xdat)[2])
    }
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
    nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
    .np_progress_bandwidth_set_total(nmulti)

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (length(bws$bw) != dim(xdat)[2])
      stop("length of bandwidth vector does not match number of columns of 'xdat'")

    npValidateRegressionNnLowerBound(
      bws,
      where = "npregbw",
      allow.zero.placeholder = TRUE
    )

    if ((any(bws$icon) &&
         !all(vapply(xdat[, bws$icon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iord) &&
         !all(vapply(xdat[, bws$iord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$iuno) &&
         !all(vapply(xdat[, bws$iuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'xdat' in type")

    if (dim(xdat)[1] != length(ydat))
      stop("number of regression data and response data do not match")

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- attr(na.omit(data.frame(xdat,ydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]
    
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
    penalty_mode <- (if (invalid.penalty == "baseline") 1L else 0L)

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
    degree.c <- if (tbw$ncon > 0) {
      as.integer(if (is.null(reg.c$degree)) rep.int(0L, tbw$ncon) else reg.c$degree)
    } else {
      integer(1)
    }

    if (bandwidth.compute){
      myopti = list(num_obs_train = dim(xdat)[1], 
        iMultistart = (if (nmulti==0) IMULTI_FALSE else IMULTI_TRUE),
        iNum_Multistart = nmulti,
        int_use_starting_values = (if (all(bws$bw==0)) USE_START_NO else USE_START_YES),
        int_LARGE_SF = (if (bws$scaling) SF_NORMAL else SF_ARB),
        BANDWIDTH_reg_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=(if (remin) RE_MIN_TRUE else RE_MIN_FALSE), 
        int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE, 
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
        int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
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
                PACKAGE = "np"))[1]
      

      rorder = numeric(ncol)
      ord_idx <- seq_len(ncol)
      rorder[c(ord_idx[bws$icon], ord_idx[bws$iuno], ord_idx[bws$iord])] <- ord_idx

      tbw$bw <- myout$bw[rorder]
      tbw$fval <- myout$fval[1]
      tbw$ifval <- myout$fval[2]
      tbw$num.feval <- sum(myout$eval.history[is.finite(myout$eval.history)])
      tbw$num.feval.fast <- myout$fast.history[1]
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

.npregbw_build_rbandwidth <- function(xdat,
                                      ydat,
                                      bws,
                                      bandwidth.compute,
                                      reg.args,
                                      yname) {
  bw.args <- c(
    list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      xnames = names(xdat),
      ynames = yname,
      bandwidth.compute = bandwidth.compute
    ),
    reg.args
  )
  do.call(rbandwidth, bw.args)
}

.npregbw_eval_only <- function(xdat,
                               ydat,
                               bws,
                               invalid.penalty = c("baseline", "dbmax"),
                               penalty.multiplier = 10) {
  invalid.penalty <- match.arg(invalid.penalty)

  xdat <- toFrame(xdat)
  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector")

  if (length(bws$bw) != dim(xdat)[2])
    stop("length of bandwidth vector does not match number of columns of 'xdat'")

  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE

  xdat <- xdat[keep.rows,, drop = FALSE]
  ydat <- ydat[keep.rows]

  if (is.factor(ydat))
    ydat <- dlev(ydat)[as.integer(ydat)]
  else
    ydat <- as.double(ydat)

  xmat <- toMatrix(xdat)
  runo <- xmat[, bws$iuno, drop = FALSE]
  rcon <- xmat[, bws$icon, drop = FALSE]
  rord <- xmat[, bws$iord, drop = FALSE]

  mysd <- EssDee(rcon)
  nrow <- dim(xmat)[1L]
  nconfac <- nrow^(-1.0 / (2.0 * bws$ckerorder + bws$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * bws$ckerorder + bws$ncon))

  penalty_mode <- if (match.arg(invalid.penalty) == "baseline") 1L else 0L
  reg.c <- npRegtypeToC(regtype = bws$regtype,
                        degree = bws$degree,
                        ncon = bws$ncon,
                        context = "npregbw")
  degree.c <- if (bws$ncon > 0) as.integer(bws$degree) else integer(1L)

  myopti <- list(
    num_obs_train = nrow,
    iMultistart = IMULTI_FALSE,
    iNum_Multistart = 0L,
    int_use_starting_values = USE_START_YES,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_reg_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    itmax = 0L,
    int_RESTART_FROM_MIN = RE_MIN_FALSE,
    int_MINIMIZE_IO = IO_MIN_TRUE,
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
      racineliyan = OKER_RLY),
    nuno = bws$nuno,
    nord = bws$nord,
    ncon = bws$ncon,
    regtype = reg.c$code,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    scale.init.categorical.sample = FALSE,
    dfc.dir = 0,
    transform.bounds = FALSE
  )

  myoptd <- list(
    ftol = 0,
    tol = 0,
    small = 0,
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
    ncatfac = ncatfac
  )

  cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

  out <- .Call(
    "C_np_regression_bw_eval",
    as.double(runo),
    as.double(rord),
    as.double(rcon),
    as.double(ydat),
    as.double(mysd),
    as.integer(myopti),
    as.double(myoptd),
    as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
    as.integer(1L),
    as.integer(penalty_mode),
    as.double(penalty.multiplier),
    as.integer(degree.c),
    as.integer(isTRUE(bws$bernstein.basis)),
    as.integer(npLpBasisCode(bws$basis)),
    as.double(cker.bounds.c$lb),
    as.double(cker.bounds.c$ub),
    PACKAGE = "np"
  )

  list(
    objective = as.numeric(out$fval[1L]),
    num.feval = 1L
  )
}

.npregbw_nomad_controls <- function(search.engine) {
  .np_degree_search_engine_controls(search.engine)
}

.npregbw_nomad_bw_setup <- function(xdat, template, bandwidth.scale.categorical = 1e4) {
  xdat <- toFrame(xdat)
  xmat <- toMatrix(xdat)
  rcon <- xmat[, template$icon, drop = FALSE]
  mysd <- EssDee(rcon)
  nrow <- dim(xmat)[1L]
  nconfac <- nrow^(-1.0 / (2.0 * template$ckerorder + template$ncon))
  ncatfac <- nrow^(-2.0 / (2.0 * template$ckerorder + template$ncon))

  cont_scale <- mysd * nconfac
  cont_idx <- which(template$icon)
  cat_idx <- c(which(template$iuno), which(template$iord))

  cat_upper <- numeric(length(cat_idx))
  for (k in seq_along(cat_idx)) {
    i <- cat_idx[k]
    if (template$iuno[i] && identical(template$ukertype, "aitchisonaitken")) {
      nlev <- length(unique(xdat[[i]]))
      cat_upper[k] <- (nlev - 1) / nlev
    } else {
      cat_upper[k] <- 1
    }
  }

  list(
    cont_scale = cont_scale,
    cont_idx = cont_idx,
    cat_idx = cat_idx,
    ncatfac = ncatfac,
    bandwidth.scale.categorical = bandwidth.scale.categorical,
    cat_upper = cat_upper
  )
}

.npregbw_nomad_point_to_bw <- function(point, template, setup) {
  point <- as.numeric(point)
  bws <- numeric(length(template$bw))
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)

  if (ncon > 0L) {
    gamma <- point[seq_len(ncon)]
    ext_bw <- gamma * setup$cont_scale
    bws[setup$cont_idx] <- if (isTRUE(template$scaling)) gamma else ext_bw
  }

  if (ncat > 0L) {
    lambda_scaled <- point[ncon + seq_len(ncat)]
    ext_bw <- lambda_scaled / setup$bandwidth.scale.categorical
    bws[setup$cat_idx] <- if (isTRUE(template$scaling)) ext_bw / setup$ncatfac else ext_bw
  }

  bws
}

.npregbw_nomad_bw_to_point <- function(bws, template, setup) {
  point <- numeric(length(setup$cont_idx) + length(setup$cat_idx))

  if (length(setup$cont_idx) > 0L) {
    raw <- bws[setup$cont_idx]
    point[seq_along(setup$cont_idx)] <- if (isTRUE(template$scaling)) {
      raw
    } else {
      raw / setup$cont_scale
    }
  }

  if (length(setup$cat_idx) > 0L) {
    raw <- bws[setup$cat_idx]
    ext_bw <- if (isTRUE(template$scaling)) raw * setup$ncatfac else raw
    point[length(setup$cont_idx) + seq_along(setup$cat_idx)] <- ext_bw * setup$bandwidth.scale.categorical
  }

  point
}

.npregbw_nomad_search <- function(xdat,
                                  ydat,
                                  bws,
                                  reg.args,
                                  opt.args,
                                  yname,
                                  degree.search,
                                  random.seed = 42L) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template <- .npregbw_build_rbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args,
    yname = yname
  )
  if (!identical(template$type, "fixed"))
    stop("automatic degree search with search.engine='nomad' currently requires bwtype='fixed'")

  setup <- .npregbw_nomad_bw_setup(xdat = xdat, template = template)
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) 1L else max(1L, as.integer(opt.args$nmulti[1L]))

  bw_lower <- c(rep.int(1e-2, ncon), rep.int(0, ncat))
  bw_upper <- c(rep.int(1e6, ncon), setup$cat_upper * setup$bandwidth.scale.categorical)

  point.start <- if (all(template$bw == 0)) {
    NULL
  } else {
    .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup)
  }
  x0 <- c(
    .np_nomad_complete_start_point(
      point = point.start,
      lower = bw_lower,
      upper = bw_upper,
      ncont = ncon
    ),
    as.integer(degree.search$start.degree)
  )
  lb <- c(bw_lower, degree.search$lower)
  ub <- c(bw_upper, degree.search$upper)
  bbin <- c(rep.int(0L, ncon + ncat), rep.int(1L, ncon))
  baseline.record <- NULL

  .np_nomad_baseline_note(degree.search$start.degree)

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[ncon + ncat + seq_len(ncon)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)

    eval.reg.args <- reg.args
    eval.reg.args$regtype <- "lp"
    eval.reg.args$degree <- degree
    eval.reg.args$bernstein.basis <- degree.search$bernstein.basis

    tbw <- .npregbw_build_rbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = eval.reg.args,
      yname = yname
    )

    out <- .npregbw_eval_only(
      xdat = xdat,
      ydat = ydat,
      bws = tbw,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier
    )

    list(
      objective = out$objective,
      degree = degree,
      num.feval = out$num.feval
    )
  }

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

      tbw <- .npregbw_build_rbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bw_vec,
        bandwidth.compute = FALSE,
        reg.args = final.reg.args,
        yname = yname
      )
      tbw$fval <- as.numeric(best_record$objective)
      tbw$ifval <- as.numeric(best_record$objective)
      tbw$num.feval <- if (!is.null(solution$bbe)) as.numeric(solution$bbe) else as.numeric(best_record$num.feval)
      tbw$num.feval.fast <- 0
      tbw$fval.history <- as.numeric(best_record$objective)
      tbw$eval.history <- if (!is.null(solution$bbe)) rep(1, max(1L, as.integer(solution$bbe))) else 1
      tbw$invalid.history <- 0
      tbw$timing <- NA_real_
      tbw$total.time <- NA_real_

      npregbw.rbandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = tbw,
        bandwidth.compute = FALSE
      )
    }

    direct.payload <- build_direct_payload()
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      .np_nomad_powell_note(degree)
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.opt.args <- opt.args
      hot.opt.args$nmulti <- 0L
      powell.start <- proc.time()[3L]
      hot.payload <- local({
        .np_progress_bandwidth_set_context(
          sprintf("Powell hot start deg %s", .np_degree_format_degree(degree))
        )
        on.exit(.np_progress_bandwidth_set_context(NULL), add = TRUE)
        .npregbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          bws = bw_vec,
          reg.args = hot.reg.args,
          opt.args = hot.opt.args,
          yname = yname
        )
      })
      powell.elapsed <- proc.time()[3L] - powell.start
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min")) {
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

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
    random.seed = random.seed,
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

.npregbw_run_fixed_degree <- function(xdat, ydat, bws, reg.args, opt.args, yname) {
  tbw <- .npregbw_build_rbandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args,
    yname = yname
  )

  do.call(npregbw.rbandwidth, c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args))
}

.npregbw_degree_search_controls <- function(regtype,
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
  search.engine <- .npregbw_nomad_controls(search.engine)

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (ncon < 1L)
    stop("automatic degree search requires at least one continuous regressor")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = ncon,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  if (!isTRUE(bern.auto) && any(bounds$upper > 3L))
    stop("automatic degree search with bernstein.basis=FALSE currently requires degree.max <= 3")

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
      stop("degree.start must lie within the searched degree candidates for every continuous predictor")
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

.npregbw_attach_degree_search <- function(bws, search_result) {
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
  bws <- .np_attach_nomad_restart_summary(bws, search_result)
  bws$degree.search <- metadata
  bws
}

npregbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws,
           bandwidth.compute = TRUE,
           basis,
           bernstein.basis,
           bwmethod,
           bwscaling,
           bwtype,
           cfac.dir,
           cfac.init,
           ckerbound,
           ckerlb,
           ckerorder,
           ckertype,
           ckerub,
           degree,
           degree.select = c("manual", "coordinate", "exhaustive"),
           search.engine = c("nomad+powell", "cell", "nomad"),
           degree.min = NULL,
           degree.max = NULL,
           degree.start = NULL,
           degree.restarts = 0L,
           degree.max.cycles = 20L,
           degree.verify = FALSE,
           dfac.dir,
           dfac.init,
           dfc.dir,
           ftol,
           hbc.init,
           hbd.dir,
           hbd.init,
           initc.dir,
           initd.dir,
           invalid.penalty = c("baseline","dbmax"),
           itmax,
           lbc.dir,
           lbc.init,
           lbd.dir,
           lbd.init,
           nmulti,
           okertype,
           penalty.multiplier = 10,
           regtype,
           remin,
           scale.init.categorical.sample,
           small,
           tol,
           transform.bounds = FALSE,
           ukertype,
           ...){

    xdat <- toFrame(xdat)
    yname <- deparse(substitute(ydat))
    lp.dot.args <- list(...)
    npRejectLegacyLpArgs(names(lp.dot.args), where = "npregbw")
    random.seed.value <- .np_degree_extract_random_seed(lp.dot.args)

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("regtype", "basis", "degree", "bernstein.basis", "bwmethod", "bwscaling", "bwtype",
               "ckertype", "ckerorder", "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    bw.args <- list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      xnames = names(xdat),
      ynames = yname,
      bandwidth.compute = bandwidth.compute
    )
    if (any.m) {
      nms <- mc.names[m]
      bw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("nmulti", "remin", "itmax", "ftol", "tol",
               "small",
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

    if (any.m) {
      nms <- mc.names[m]
      opt.args <- mget(nms, envir = environment(), inherits = FALSE)
    } else {
      opt.args <- list()
    }

    reg.args <- bw.args[setdiff(names(bw.args), c("bw", "nobs", "xdati", "ydati", "xnames", "ynames", "bandwidth.compute"))]
    opt.args <- c(list(bandwidth.compute = bandwidth.compute), opt.args)

    ncon <- sum(untangle(xdat)$icon)
    search.mc.names <- names(match.call(expand.dots = FALSE))
    regtype.value <- if ("regtype" %in% search.mc.names) regtype else "lc"
    bernstein.value <- if ("bernstein.basis" %in% search.mc.names) bernstein.basis else TRUE
    degree.select.value <- if ("degree.select" %in% search.mc.names) degree.select else "manual"
    search.engine.value <- if ("search.engine" %in% search.mc.names) search.engine else "nomad+powell"
    degree.min.value <- if ("degree.min" %in% search.mc.names) degree.min else NULL
    degree.max.value <- if ("degree.max" %in% search.mc.names) degree.max else NULL
    degree.start.value <- if ("degree.start" %in% search.mc.names) degree.start else NULL
    degree.restarts.value <- if ("degree.restarts" %in% search.mc.names) degree.restarts else 0L
    degree.max.cycles.value <- if ("degree.max.cycles" %in% search.mc.names) degree.max.cycles else 20L
    degree.verify.value <- if ("degree.verify" %in% search.mc.names) degree.verify else FALSE
    degree.search <- .npregbw_degree_search_controls(
      regtype = regtype.value,
      regtype.named = "regtype" %in% search.mc.names,
      ncon = ncon,
      nobs = NROW(xdat),
      basis = if ("basis" %in% search.mc.names) basis else "glp",
      degree.select = degree.select.value,
      search.engine = search.engine.value,
      degree.min = degree.min.value,
      degree.max = degree.max.value,
      degree.start = degree.start.value,
      degree.restarts = degree.restarts.value,
      degree.max.cycles = degree.max.cycles.value,
      degree.verify = degree.verify.value,
      bernstein.basis = bernstein.value,
      bernstein.named = "bernstein.basis" %in% search.mc.names
    )

    if (!is.null(degree.search)) {
      if (identical(degree.search$engine, "cell")) {
        eval_fun <- function(degree.vec) {
          cell.reg.args <- reg.args
          cell.reg.args$regtype <- "lp"
          cell.reg.args$degree <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
          cell.bws <- .npregbw_run_fixed_degree(
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            reg.args = cell.reg.args,
            opt.args = opt.args,
            yname = yname
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
        search.result <- .npregbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          reg.args = reg.args,
          opt.args = opt.args,
          yname = yname,
          degree.search = degree.search,
          random.seed = random.seed.value
        )
      }
      tbw <- .npregbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
      mc <- match.call(expand.dots = FALSE)
      environment(mc) <- parent.frame()
      tbw$call <- mc
      return(tbw)
    }

    tbw <- .npregbw_build_rbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = bws,
      bandwidth.compute = bandwidth.compute,
      reg.args = reg.args,
      yname = yname
    )
    bwsel.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
    tbw <- .np_progress_select_bandwidth_enhanced(
      "Selecting regression bandwidth",
      do.call(npregbw.rbandwidth, bwsel.args)
    )

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
    
  }
