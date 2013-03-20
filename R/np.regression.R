npreg <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
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
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))

    tydat <- model.response(tmf)
    txdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]

    if ((has.eval <- !is.null(newdata))) {
      if (!(has.ey <- succeedWithResponse(tt, newdata)))
        tt <- delete.response(tt)
      
      umf <- emf <- model.frame(tt, data = newdata)

      if (has.ey)
        eydat <- model.response(emf)
      
      exdat <- emf[, attr(attr(emf, "terms"),"term.labels"), drop = FALSE]
    }

    ev <- eval(parse(text=paste("npreg(txdat = txdat, tydat = tydat,",
                       ifelse(has.eval,paste("exdat = exdat,",ifelse(has.ey,"eydat = eydat,","")),""),
                       "bws = bws, ...)")))
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    return(ev)
  }

npreg.call <-
  function(bws, ...) {
    ev <- npreg(txdat = eval(bws$call[["xdat"]], environment(bws$call)),
                tydat = eval(bws$call[["ydat"]], environment(bws$call)),
                bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npreg.rbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat, eydat, gradients = FALSE, residuals = FALSE,
           ...){

    no.ex = missing(exdat)
    no.ey = missing(eydat)

    txdat = toFrame(txdat)

    if (!(is.vector(tydat) | is.factor(tydat)))
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
        if (!(is.vector(eydat) | is.factor(eydat)))
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

    ccon = unlist(lapply(txdat[,bws$icon, drop = FALSE],class))
    if ((any(bws$icon) && !all((ccon == class(integer(0))) | (ccon == class(numeric(0))))) ||
        (any(bws$iord) && !all(unlist(lapply(txdat[,bws$iord, drop = FALSE],class)) ==
                               class(ordered(0)))) ||
        (any(bws$iuno) && !all(unlist(lapply(txdat[,bws$iuno, drop = FALSE],class)) ==
                               class(factor(0)))))
      stop("supplied bandwidths do not match 'txdat' in type")

    if (dim(txdat)[1] != length(tydat))
      stop("number of explanatory data 'txdat' and dependent data 'tydat' do not match")

    ## catch and destroy NA's
    goodrows = 1:dim(txdat)[1]
    rows.omit = attr(na.omit(data.frame(txdat,tydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Training data has no rows without NAs")

    txdat = txdat[goodrows,,drop = FALSE]
    tydat = tydat[goodrows]

    ## no.ex = missing(exdat)
    ## no.ey = missing(eydat)

    if (!no.ex){
      goodrows = 1:dim(exdat)[1]
      rows.omit = eval(parse(text=paste('attr(na.omit(data.frame(exdat',
                               ifelse(no.ey,"",",eydat"),')), "na.action")')))

      goodrows[rows.omit] = 0

      exdat = exdat[goodrows,,drop = FALSE]
      if (!no.ey)
        eydat = eydat[goodrows]

      if (all(goodrows==0))
        stop("Evaluation data has no rows without NAs")
    }

    ## evaluate residuals before data conversion ...

    if (residuals){
      resid <- tydat - npreg(txdat = txdat, tydat = tydat, bws = bws)$mean
    }


    tnrow = dim(txdat)[1]
    enrow = ifelse(no.ex,tnrow,dim(exdat)[1])
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
      int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
      BANDWIDTH_reg_extern = switch(bws$type,
        fixed = BW_FIXED,
        generalized_nn = BW_GEN_NN,
        adaptive_nn = BW_ADAP_NN),
      int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
      kerneval = switch(bws$ckertype,
        gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
        epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
        uniform = CKER_UNI),
      ukerneval = switch(bws$ukertype,
        aitchisonaitken = UKER_AIT,
        liracine = UKER_LR),
      okerneval = switch(bws$okertype,
        wangvanryzin = OKER_WANG,
        liracine = OKER_LR),
      ey_is_ty = no.ey,
      do_grad = gradients,
      regtype = switch(bws$regtype,
        lc = REGTYPE_LC,
        ll = REGTYPE_LL),
      no.ex = no.ex,
      mcv.numRow = attr(bws$xmcv, "num.row"))
    

    myout=
      .C("np_regression",
         as.double(tuno), as.double(tord), as.double(tcon), as.double(tydat),
         as.double(euno),  as.double(eord),  as.double(econ), as.double(eydat),
         as.double(c(bws$bw[bws$icon],bws$bw[bws$iuno],bws$bw[bws$iord])),
         as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
         as.integer(myopti),
         mean = double(enrow),
         merr = double(enrow),
         g = double(ifelse(gradients,enrow*ncol,0)),
         gerr = double(ifelse(gradients,enrow*ncol,0)),
         xtra = double(6),
         PACKAGE="npRmpi" )[c("mean","merr", "g", "gerr", "xtra")]

    if (gradients){
      myout$g = matrix(data=myout$g, nrow = enrow, ncol = ncol, byrow = FALSE) 
      rorder = numeric(ncol)
      rorder[c((1:ncol)[bws$icon], (1:ncol)[bws$iuno], (1:ncol)[bws$iord])]=1:ncol
      myout$g = as.matrix(myout$g[,rorder])

      myout$gerr = matrix(data=myout$gerr, nrow = enrow, ncol = ncol, byrow = FALSE) 
      myout$gerr = as.matrix(myout$gerr[,rorder])
    }


    ev <- eval(parse(text = paste("npregression(bws = bws,",
                       "eval = teval,",
                       "mean = myout$mean, merr = myout$merr,",
                       ifelse(gradients,
                              "grad = myout$g, gerr = myout$gerr,",""),
                       ifelse(residuals, "resid = resid,", ""),
                       "ntrain = tnrow,",
                       "trainiseval = no.ex,",
                       "gradients = gradients,",
                       "residuals = residuals,",
                       "xtra = myout$xtra, rows.omit = rows.omit)")))

    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npreg.default <- function(bws, txdat, tydat, ...){
  sc.names <- names(sys.call())

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

  mc <- match.call()

  tx.str <- ifelse(txdat.named, "xdat = txdat,",
                   ifelse(no.txdat, "", "txdat,"))
  ty.str <- ifelse(tydat.named, "ydat = tydat,",
                   ifelse(no.tydat, "", "tydat,"))
  
  tbw <- eval(parse(text = paste("npregbw(",
                      ifelse(bws.named,                             
                             paste(tx.str, ty.str,
                                   "bws = bws, bandwidth.compute = FALSE,"),
                             paste(ifelse(no.bws, "", "bws,"), tx.str, ty.str)),
                      "call = mc, ...",")",sep="")))

  ## xnames = names(txdat)
  ##tbw <- updateBwNameMetadata(nameList =
  ##                            list(ynames = deparse(substitute(tydat))),
  ##                            bws = tbw)

  repair.args <- c("data", "subset", "na.action")
  
  m.par <- match(repair.args, names(mc), nomatch = 0)
  m.child <- match(repair.args, names(tbw$call), nomatch = 0)

  if(any(m.child > 0)) {
    tbw$call[m.child] <- mc[m.par]
  }

  ## next we repair arguments portion of the call
  m.bws.par <- match(c("bws","txdat","tydat"), names(mc), nomatch = 0)
  m.bws.child <- match(c("bws","txdat","tydat"), as.character(tbw$call), nomatch = 0)
  m.bws.union <- (m.bws.par > 0) & (m.bws.child > 0)
  
  tbw$call[m.bws.child[m.bws.union]] <- mc[m.bws.par[m.bws.union]]

  environment(tbw$call) <- parent.frame()

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
  
  ev <- eval(parse(text=paste("npreg(bws = tbw", tx.str, ty.str, ",...)")))

  ev$call <- match.call(expand.dots = FALSE)
  environment(ev$call) <- parent.frame()
  return(ev)
}

