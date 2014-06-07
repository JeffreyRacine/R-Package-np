npqreg <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npqreg",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npqreg",bws$call)
        else if (!is.call(bws))
          UseMethod("npqreg",bws)
        else
          UseMethod("npqreg",NULL)
      } else {
        UseMethod("npqreg", NULL)
      }
    } else {
      UseMethod("npqreg", NULL)
    }
  }

npqreg.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))

    tydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
    txdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]

    if ((has.eval <- !is.null(newdata))) {
      tt <- drop.terms(tt, match(bws$variableNames$response, attr(tt, 'term.labels')))
      umf <- emf <- model.frame(tt, data = newdata)
      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }

    tbw <-
    eval(parse(text=paste("npqreg(txdat = txdat, tydat = tydat,",
                 ifelse(has.eval,"exdat = exdat,",""), "bws = bws, ...)")))

    tbw$omit <- attr(umf,"na.action")
    tbw$rows.omit <- as.vector(tbw$omit)
    tbw$nobs.omit <- length(tbw$rows.omit)

    tbw$quantile <- napredict(tbw$omit, tbw$quantile)
    tbw$quanterr <- napredict(tbw$omit, tbw$quanterr)

    if(tbw$gradients){
        tbw$quantgrad <- napredict(tbw$omit, tbw$quantgrad)
    }

    return(tbw)
  }

npqreg.call <-
  function(bws, ...) {
    npqreg(txdat = eval(bws$call[["xdat"]], environment(bws$call)),
           tydat = eval(bws$call[["ydat"]], environment(bws$call)),
           bws = bws, ...)
  }

npqreg.conbandwidth <-
  function(bws, ...){
    stop("incorrect bandwidth type: expected conditional distribution bandwidths instead of conditional density bandwidths")
  }

npqreg.condbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           tau = 0.5,
           gradients = FALSE,
           ftol = 1.490116e-07, tol = 1.490116e-04,
           small = 1.490116e-05, itmax = 10000,
           lbc.dir = 0.5, dfc.dir = 3, cfac.dir = 2.5*(3.0-sqrt(5)),initc.dir = 1.0, 
           lbd.dir = 0.1, hbd.dir = 1, dfac.dir = 0.25*(3.0-sqrt(5)), initd.dir = 1.0, 
           ...){

    no.ex = missing(exdat)

    txdat = toFrame(txdat)
    tydat = toFrame(tydat)

    if (dim(tydat)[2] != 1)
      stop("'tydat' has more than one column")

    if (!no.ex){
      exdat = toFrame(exdat)
      
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
    }

    if(gradients)
      stop("gradients not currently supported for this object")
    
    if (length(bws$xbw) != length(txdat))
      stop("length of bandwidth vector does not match number of columns of 'txdat'")

    if (length(bws$ybw) != 1)
      stop("length of bandwidth vector does not match number of columns of 'tydat'")

    if (any(bws$iyord) | any(bws$iyuno) | coarseclass(tydat[,1]) != "numeric")
      stop("'tydat' is not continuous")
    
    xccon = unlist(lapply(txdat[,bws$ixcon, drop = FALSE],class))
    if ((any(bws$ixcon) && !all((xccon == class(integer(0))) | (xccon == class(numeric(0))))) ||
        (any(bws$ixord) && !all(unlist(lapply(txdat[,bws$ixord, drop = FALSE],class)) ==
                                class(ordered(0)))) ||
        (any(bws$ixuno) && !all(unlist(lapply(txdat[,bws$ixuno, drop = FALSE],class)) ==
                                class(factor(0)))))
      stop("supplied bandwidths do not match 'txdat' in type")

    ## catch and destroy NA's
    goodrows = 1:dim(txdat)[1]
    rows.omit = attr(na.omit(data.frame(txdat,tydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Training data has no rows without NAs")

    txdat = txdat[goodrows,,drop = FALSE]
    tydat = tydat[goodrows,,drop = FALSE]

    if (!no.ex)
      exdat = na.omit(exdat)
    
    tnrow = dim(txdat)[1]
    enrow = ifelse(no.ex,tnrow,dim(exdat)[1])

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati)
    tydat <- adjustLevels(tydat, bws$ydati)
    
    if (!no.ex){
      exdat <- adjustLevels(exdat, bws$xdati)
    }

    ## grab the evaluation data before it is converted to numeric
    if(no.ex){
      txeval <- txdat
    } else {
      txeval <- exdat
    }

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.
    
    tydat = toMatrix(tydat)

    txdat = toMatrix(txdat)

    txuno = txdat[, bws$ixuno, drop = FALSE]
    txcon = txdat[, bws$ixcon, drop = FALSE]
    txord = txdat[, bws$ixord, drop = FALSE]

    if (!no.ex){
      exdat = toMatrix(exdat)

      exuno = exdat[, bws$ixuno, drop = FALSE]
      excon = exdat[, bws$ixcon, drop = FALSE]
      exord = exdat[, bws$ixord, drop = FALSE]
    } else {
      exuno = data.frame()
      excon = data.frame()
      exord = data.frame()
    }

    myopti = list(
      num_obs_train = tnrow,
      num_obs_eval = enrow,
      int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
      BANDWIDTH_den_extern = switch(bws$type,
        fixed = BW_FIXED,
        generalized_nn = BW_GEN_NN,
        adaptive_nn = BW_ADAP_NN),
      int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE),
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
        liracine = OKER_LR),
      oykerneval = switch(bws$oykertype,
        wangvanryzin = OKER_WANG,
        liracine = OKER_NLR),
      num_yuno = bws$ynuno,
      num_yord = bws$ynord,
      num_ycon = bws$yncon,
      num_xuno = bws$xnuno,
      num_xord = bws$xnord,
      num_xcon = bws$xncon,
      no.ex = no.ex,
      gradients = gradients,
      itmax = itmax,
      xmcv.numRow = attr(bws$xmcv, "num.row"),
      nmulti = itmax,
      dfc.dir = dfc.dir)

    myoptd = list(
      ftol = ftol,
      tol = tol,
      small = small,
      lbc.dir = lbc.dir, cfac.dir = cfac.dir,initc.dir = initc.dir, 
      lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir)
    
    myout=
      .C("np_quantile_conditional",
         as.double(tydat),
         as.double(txuno), as.double(txord), as.double(txcon),
         as.double(exuno), as.double(exord), as.double(excon),
         as.double(tau),
         as.double(c(bws$xbw[bws$ixcon],bws$ybw[bws$iycon],
                     bws$ybw[bws$iyuno],bws$ybw[bws$iyord],
                     bws$xbw[bws$ixuno],bws$xbw[bws$ixord])),
         as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
         as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
         as.integer(myopti),
         as.double(myoptd),
         yq = double(enrow),
         yqerr = double(enrow),
         yqgrad = double(enrow*bws$xndim*gradients),
         PACKAGE="np" )[c("yq","yqerr", "yqgrad")]

    ##need to untangle yqgrad

    if(gradients){
      myout$yqgrad = matrix(data=myout$yqgrad, nrow = enrow, ncol = bws$xndim, byrow = FALSE) 
      rorder = numeric(bws$xndim)
      rorder[c((1:bws$xndim)[bws$ixcon], (1:bws$xndim)[bws$ixuno], (1:bws$xndim)[bws$ixord])]=1:bws$xndim
      myout$yqgrad = myout$yqgrad[, rorder, drop = FALSE]

    } else {
      myout$yqgrad = NA
    }


    qregression(bws = bws,
                xeval = txeval,
                tau = tau,
                quantile = myout$yq,
                quanterr = myout$yqerr,
                quantgrad = myout$yqgrad,
                ntrain = tnrow,
                trainiseval = no.ex,
                gradients = gradients)
  }


npqreg.default <- function(bws, txdat, tydat, ...){
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

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npcdistbw)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat')
  nstxy <- c('xdat','ydat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }
    
  tbw <- eval.parent(sc.bw)

  ## convention: drop 'bws' and up to two unnamed arguments (including bws)
  if(no.bws){
    tx.str <- ",txdat = txdat"
    ty.str <- ",tydat = tydat"
  } else {
    tx.str <- ifelse(txdat.named, ",txdat = txdat","")
    ty.str <- ifelse(tydat.named, ",tydat = tydat","")    
    if((!bws.named) && (!txdat.named)){
      ty.str <- ifelse(tydat.named, ",tydat = tydat",
                       ifelse(no.tydat,"",",tydat"))
    }
  }
  
  eval(parse(text=paste("npqreg(bws = tbw", tx.str, ty.str, ",...)")))
}

