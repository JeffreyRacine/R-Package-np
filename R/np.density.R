npudens <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$tdat))
          UseMethod("npudens",bws$formula)
        else if (!is.null(bws$call) && is.null(args$tdat))
          UseMethod("npudens",bws$call)
        else if (!is.call(bws))
          UseMethod("npudens",bws)
        else
          UseMethod("npudens",NULL)
      } else {
        UseMethod("npudens", NULL)
      }
    } else {
      UseMethod("npudens", NULL)
    }
  }

npudens.formula <-
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

    ev <- 
      eval(parse(text=paste("npudens(tdat = tdat,",
                   ifelse(has.eval,"edat = edat,",""), "bws = bws, ...)")))
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    return(ev)
  }

npudens.call <-
  function(bws, ...) {
    npudens(bws, tdat = eval(bws$call[["dat"]], environment(bws$call)),
            ...)
  }

npudens.bandwidth <-
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

  tdat <- na.omit(tdat)
  rows.omit <- unclass(na.action(tdat))

  if (!no.e){
    edat <- na.omit(edat)
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
    kerneval = switch(bws$ckertype,
      gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
      uniform = CKER_UNI),
    no.e = no.e,
    mcv.numRow = attr(bws$xmcv, "num.row"),
    densOrDist = NP_DO_DENS)
  
  myout=
    .C("np_density", as.double(tuno), as.double(tord), as.double(tcon),
       as.double(euno),  as.double(eord),  as.double(econ), 
       as.double(c(bws$bw[bws$icon],bws$bw[bws$iuno],bws$bw[bws$iord])),
       as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
       as.integer(myopti),
       dens = double(enrow),
       derr = double(enrow),
       log_likelihood = double(1),
       PACKAGE="npRmpi" )[c("dens","derr", "log_likelihood")]

  ev <- npdensity(bws=bws, eval=teval, dens = myout$dens,
                  derr = myout$derr, ll = myout$log_likelihood,
                  ntrain = tnrow, trainiseval = no.e,
                  rows.omit = rows.omit)
  return(ev)
}

npudens.default <- function(bws, tdat, ...){
  sc.names <- names(sys.call())

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

  mc <- match.call()
  
  tbw <- eval(parse(text = paste("npudensbw(",
                      ifelse(tdat.named, "dat = tdat",
                             ifelse(no.tdat,"","tdat")),
                      ifelse(no.tdat,"",","),
                      ifelse(bws.named,"bws = bws, bandwidth.compute = FALSE",
                             ifelse(no.bws,"","bws")),
                      ifelse(no.bws,"",","),                      
                      "call = mc, ...",")",sep="")))

  ## need to do some surgery on the call to
  ## allow it to work with the formula interface

  repair.args <- c("data", "subset", "na.action")
  
  m.par <- match(repair.args, names(mc), nomatch = 0)
  m.child <- match(repair.args, names(tbw$call), nomatch = 0)

  if(any(m.child > 0)) {
    tbw$call[m.child] <- mc[m.par]
  }

  ## next we repair 'bws' portion of the call
  m.bws.par <- match(c("bws","tdat"), names(mc), nomatch = 0)
  m.bws.child <- match(c("bws","tdat"), as.character(tbw$call), nomatch = 0)
  m.bws.union <- (m.bws.par > 0) & (m.bws.child > 0)
  
  tbw$call[m.bws.child[m.bws.union]] <- mc[m.bws.par[m.bws.union]]
  
  environment(tbw$call) <- parent.frame()

  ## convention: first argument is always dropped, second, if present, propagated
  eval(parse(text=paste("npudens(bws = tbw",
               ifelse(no.tdat, "",
                      ifelse(tdat.named, ",tdat = tdat",",tdat")),
               ",...)")))
}

