npudist <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$tdat))
          UseMethod("npudist",bws$formula)
        else if (!is.null(bws$call) && is.null(args$tdat))
          UseMethod("npudist",bws$call)
        else if (!is.call(bws))
          UseMethod("npudist",bws)
        else
          UseMethod("npudist",NULL)
      } else {
        UseMethod("npudist", NULL)
      }
    } else {
      UseMethod("npudist", NULL)
    }
  }


npudist.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))

    tdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]

    if ((has.eval <- !is.null(newdata))) {
      umf <- emf <- model.frame(tt, data = newdata)

      edat <- emf[, attr(attr(emf, "terms"),"term.labels"), drop = FALSE]
    }
    
    ev <- eval(parse(text=paste("npudist(tdat = tdat,",
                 ifelse(has.eval,"edat = edat,",""), "bws = bws, ...)")))

    ev$omit <- attr(umf,"na.action")
    ev$rows.omit <- as.vector(ev$omit)
    ev$nobs.omit <- length(ev$rows.omit)

    ev$dist <- napredict(ev$omit, ev$dist)
    ev$derr <- napredict(ev$omit, ev$derr)

    return(ev)
  }

npudist.call <-
  function(bws, ...) {
    npudist(tdat = eval(bws$call[["dat"]], environment(bws$call)),
          bws = bws, ...)
  }

npudist.dbandwidth <-
  function(bws,
           tdat = stop("invoked without training data 'tdat'"),
           edat, ...){

    no.e = missing(edat)

    tdat = toFrame(tdat)

    if (!no.e)
      edat = toFrame(edat)

    if (!(no.e || tdat %~% edat ))
      stop("tdat and edat are not similar data frames!")

    if (length(bws$bw) != length(tdat))
      stop("length of bandwidth vector does not match number of columns of 'tdat'")

    ccon = unlist(lapply(as.data.frame(tdat[,bws$icon]),class))
    if ((any(bws$icon) && !all((ccon == class(integer(0))) | (ccon == class(numeric(0))))) ||
        (any(bws$iord) && !all(unlist(lapply(as.data.frame(tdat[,bws$iord]),class)) ==
                               class(ordered(0)))) ||
        (any(bws$iuno) && !all(unlist(lapply(as.data.frame(tdat[,bws$iuno]),class)) ==
                               class(factor(0)))))
      stop("supplied bandwidths do not match 'tdat' in type")

    tdat = na.omit(tdat)
    rows.omit <- unclass(na.action(tdat))

    if (!no.e){
      edat = na.omit(edat)
      rows.omit <- unclass(na.action(edat))
    }

    tnrow = nrow(tdat)
    enrow = ifelse(no.e,tnrow,nrow(edat))

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    tdat <- adjustLevels(tdat, bws$xdati)
    
    if (!no.e)
      edat <- adjustLevels(edat, bws$xdati, allowNewCells = TRUE)

    ## grab the evaluation data before it is converted to numeric
    if(no.e)
      teval <- tdat
    else
      teval <- edat

    ## put the unordered, ordered, and continuous data in their own objects
    ## data that is not a factor is continuous.
    
    tdat = toMatrix(tdat)

    tuno = tdat[, bws$iuno, drop = FALSE]
    tcon = tdat[, bws$icon, drop = FALSE]
    tord = tdat[, bws$iord, drop = FALSE]

    if (!no.e){
      edat = toMatrix(edat)

      euno = edat[, bws$iuno, drop = FALSE]
      econ = edat[, bws$icon, drop = FALSE]
      eord = edat[, bws$iord, drop = FALSE]

    } else {
      euno = data.frame()
      eord = data.frame()
      econ = data.frame()
    }

    myopti = list(
      num_obs_train = tnrow,
      num_obs_eval = enrow,
      num_uno = bws$nuno,
      num_ord = bws$nord,
      num_con = bws$ncon,
      int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
      BANDWIDTH_den_extern = switch(bws$type,
        fixed = BW_FIXED,
        generalized_nn = BW_GEN_NN,
        adaptive_nn = BW_ADAP_NN),
      int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE),
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
        liracine = OKER_NLR),
      no.e = no.e,
      mcv.numRow = attr(bws$xmcv, "num.row"),
      densOrDist = NP_DO_DIST,
      old.dist = FALSE,
      int_do_tree = ifelse(options('np.tree'), DO_TREE_YES, DO_TREE_NO))

    
    myout=
      .C("np_density", as.double(tuno), as.double(tord), as.double(tcon),
         as.double(euno),  as.double(eord),  as.double(econ), 
         as.double(c(bws$bw[bws$icon],bws$bw[bws$iuno],bws$bw[bws$iord])),
         as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
         as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
         as.integer(myopti),
         dist = double(enrow),
         derr = double(enrow),
         log_likelihood = double(1),
         PACKAGE="np" )[c("dist","derr", "log_likelihood")]

    ev <- npdistribution(bws=bws, eval=teval, dist = myout$dist,
                         derr = myout$derr, ntrain = tnrow, trainiseval = no.e,
                         rows.omit = rows.omit)
    return(ev)
  }

npudist.default <- function(bws, tdat, ...){
  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  tdat.named <- any(sc.names == "tdat")

  no.bws <- missing(bws)
  no.tdat <- missing(tdat)

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(tdat.named)
    tdat <- toFrame(tdat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npudistbw)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('tdat')
  nstxy <- c('dat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }
    
  tbw <- eval.parent(sc.bw)

  ## convention: first argument is always dropped, second, if present, propagated
  eval(parse(text=paste("npudist(bws = tbw",
               ifelse(no.tdat, "",
                      ifelse(tdat.named, ",tdat = tdat",",tdat")),
               ",...)")))
}
