npplreg <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npplreg",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npplreg",bws$call)
        else if (!is.call(bws))
          UseMethod("npplreg",bws)
        else
          UseMethod("npplreg",NULL)
      } else {
        UseMethod("npplreg", NULL)
      }
    } else {
      UseMethod("npplreg", NULL)
    }
  }

npplreg.formula <-
  function(bws, data = NULL, newdata = NULL, ...){
    
    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))
    
    tydat <- model.response(tmf)
    txdat <- tmf[, bws$chromoly[[2]], drop = FALSE]
    tzdat <- tmf[, bws$chromoly[[3]], drop = FALSE]

    if ((has.eval <- !is.null(newdata))) {
      if (!(has.ey <- succeedWithResponse(tt, newdata)))
        tt <- delete.response(tt)
      
      umf <- emf <- model.frame(tt, data = newdata)

      if (has.ey)
        eydat <- model.response(emf)

      exdat <- emf[, bws$chromoly[[2]], drop = FALSE]
      ezdat <- emf[, bws$chromoly[[3]], drop = FALSE]
    }

    ev <-
      eval(parse(text=paste("npplreg(txdat = txdat, tydat = tydat, tzdat = tzdat,",
                   ifelse(has.eval,paste("exdat = exdat, ezdat = ezdat,",
                                         ifelse(has.ey,"eydat = eydat,","")),""),
                   "bws = bws, ...)")))
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    ev
  }

npplreg.call <-
  function(bws, ...) {
    npplreg(txdat = eval(bws$call[["xdat"]], environment(bws$call)),
            tydat = eval(bws$call[["ydat"]], environment(bws$call)),
            tzdat = eval(bws$call[["zdat"]], environment(bws$call)),
            bws = bws, ...)
  }


npplreg.plbandwidth <- 
  function(bws,
           txdat = stop("training data txdat missing"),
           tydat = stop("training data tydat missing"),
           tzdat = stop("training data tzdat missing"),
           exdat, eydat, ezdat, residuals = FALSE, ...){

    txdat = toFrame(txdat)
    tzdat = toFrame(tzdat)
    
    ## catch and destroy NA's, part 1
    goodrows = 1:dim(txdat)[1]
    rows.omit = attr(na.omit(data.frame(txdat,tydat,tzdat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Training data has no rows without NAs")

    txdat = txdat[goodrows,,drop = FALSE]
    tydat = tydat[goodrows]
    tzdat = tzdat[goodrows,,drop = FALSE]

    no.exz = missing(exdat)
    no.ey = missing(eydat)

    if (!no.exz){
      exdat = toFrame(exdat)
      ezdat = toFrame(ezdat)
      if (!no.ey)
        eydat = as.double(eydat)

      ## c& d NA's, part 2

      goodrows = 1:dim(exdat)[1]
      rows.omit = eval(parse(text=paste('attr(na.omit(data.frame(exdat,',
                               ifelse(no.ey,"","eydat,"),'ezdat)), "na.action")')))
      goodrows[rows.omit] = 0

      if (all(goodrows==0))
        stop("Evaluation data has no rows without NAs")

      exdat = exdat[goodrows,,drop = FALSE]
      if (!no.ey)
        eydat = eydat[goodrows]
      ezdat = ezdat[goodrows,,drop = FALSE]
    }

    ## tmp.ty and tmp.ey are the numeric representations of tydat and eydat
    if (is.factor(tydat)){
      tmp.ty <- adjustLevels(data.frame(tydat), bws$bw$yzbw$ydati)
      tmp.ty <- (bws$bw$yzbw$ydati$all.dlev[[1]])[as.integer(tmp.ty)]
    } else {
      tmp.ty <- as.double(tydat)
    }

    if (!no.ey){
      if (is.factor(tydat)){
        tmp.ey <- adjustLevels(data.frame(eydat), bws$bw$yzbw$ydati)
        tmp.ey <- (bws$bw$yzbw$ydati$all.dlev[[1]])[as.integer(tmp.ey)]
      } else {
        tmp.ey <- as.double(eydat)
      }
    }
    
    ## y on z
    mmy = npreg(txdat = tzdat, tydat = tydat, bws = bws$bw$yzbw)

    resy <- tmp.ty - mmy$mean

    if (!no.exz)
      mmy.eval = npreg(txdat = tzdat, tydat = tydat, exdat = ezdat, bws = bws$bw$yzbw)

    
    ## x on z
    nrow = nrow(txdat)
    nrow.eval = ifelse(no.exz,0,nrow(exdat))
    ncol = ncol(txdat)
    B = double(ncol)
    resx = matrix(data = 0, nrow = nrow, ncol = ncol)
    resx.eval = matrix(data = 0, nrow = nrow.eval, ncol = ncol)

    for (i in 1:ncol){
      mm = npreg(txdat=tzdat, tydat=txdat[,i], bws = bws$bw[[i+1]])

      if (is.factor(txdat[1,i])){
        tmp.dat <- adjustLevels(txdat[,i, drop=FALSE], bws$bw[[i+1]]$ydati)
        resx[,i] <- (bws$bw[[i+1]]$ydati$all.dlev[[1]])[as.integer(tmp.dat[,1])] - mm$mean
      } else {
        resx[,i] <- txdat[,i] - mm$mean
      }

      if(!no.exz) {
        mm = npreg(txdat=tzdat, tydat=txdat[,i], exdat=ezdat, bws = bws$bw[[i+1]])

        if (is.factor(txdat[1,i])){
          tmp.dat <- adjustLevels(exdat[,i, drop=FALSE], bws$bw[[i+1]]$ydati)
          resx.eval[,i] <- (bws$bw[[i+1]]$ydati$all.dlev[[1]])[as.integer(tmp.dat[,1])] - mm$mean
        } else {
          resx.eval[,i] <- exdat[,i] - mm$mean
        }
      }
    }

    B = coef((model = lm(resy ~ resx - 1)))

    ## computes the standard errors of B using the model matrix
    ## and the MSE of the training data predictions

    Bvcov = sum((tmp.ty-(mmy$mean + resx  %*% B))^2)/
      (dim(txdat)[1]-dim(txdat)[2]-dim(tzdat)[2])*
      solve(t(model.matrix(model))%*%model.matrix(model))

    Berr = sqrt(diag(Bvcov))

    if (!no.ey) {
      ply = mmy.eval$mean + resx.eval %*% B
      RSQ = RSQfunc(tmp.ey,ply)
      MSE = MSEfunc(tmp.ey,ply)
      MAE = MAEfunc(tmp.ey,ply)
      MAPE = MAPEfunc(tmp.ey,ply)
      CORR = CORRfunc(tmp.ey,ply)
      SIGN = SIGNfunc(tmp.ey,ply)

    } else {
      ply =  mmy$mean + resx %*% B
      RSQ = RSQfunc(tmp.ty,ply)
      MSE = MSEfunc(tmp.ty,ply)
      MAE = MAEfunc(tmp.ty,ply)
      MAPE = MAPEfunc(tmp.ty,ply)
      CORR = CORRfunc(tmp.ty,ply)
      SIGN = SIGNfunc(tmp.ty,ply)

      if (!no.exz)
        ply = mmy.eval$mean + resx.eval %*% B
    }

    ev <- eval(parse(text = paste("plregression(bws = bws,",
                       "xcoef = B, xcoeferr = Berr, xcoefvcov = Bvcov,",
                       "evalx =  if (no.exz) txdat else exdat,",
                       "evalz =  if (no.exz) tzdat else ezdat,",
                       "mean = ply, ntrain = nrow,",
                       ifelse(residuals, "resid = tmp.ty - ply,", ""),
                       "trainiseval = no.exz,",
                       "residuals = residuals,",
                       "xtra=c(RSQ,MSE,MAE,MAPE,CORR,SIGN))")))

    
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }


npplreg.default <- function(bws, txdat, tydat, tzdat, ...) {
  sc.names <- names(sys.call())

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")
  tzdat.named <- any(sc.names == "tzdat")

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  no.tzdat <- missing(tzdat)

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  if(tydat.named)
    tzdat <- toFrame(tzdat)

  mc <- match.call()

  tx.str <- ifelse(txdat.named, "xdat = txdat,",
                   ifelse(no.txdat, "", "txdat,"))
  ty.str <- ifelse(tydat.named, "ydat = tydat,",
                   ifelse(no.tydat, "", "tydat,"))
  tz.str <- ifelse(tzdat.named, "zdat = tzdat,",
                   ifelse(no.tzdat, "", "tzdat,"))
  
  tbw <- eval(parse(text = paste("npplregbw(",
                      ifelse(bws.named,                             
                             paste(tx.str, ty.str, tz.str,
                                   "bws = bws, bandwidth.compute = FALSE,"),
                             paste(ifelse(no.bws, "", "bws,"),
                                   tx.str, ty.str, tz.str)),
                      "call = mc, ...",")",sep="")))

  ## need to do some surgery on the call to
  ## allow it to work with the formula interface

  repair.args <- c("data", "subset", "na.action")
  
  m.par <- match(repair.args, names(mc), nomatch = 0)
  m.child <- match(repair.args, names(tbw$call), nomatch = 0)

  if(any(m.child > 0)) {
    tbw$call[m.child] <- mc[m.par]
  }

  ## next we repair arguments portion of the call
  m.bws.par <- match(c("bws","txdat","tydat","tzdat"), names(mc), nomatch = 0)
  m.bws.child <- match(c("bws","txdat","tydat","tzdat"), as.character(tbw$call), nomatch = 0)
  m.bws.union <- (m.bws.par > 0) & (m.bws.child > 0)
  
  tbw$call[m.bws.child[m.bws.union]] <- mc[m.bws.par[m.bws.union]]

  environment(tbw$call) <- parent.frame()

  ## convention: drop 'bws' and up to three unnamed arguments (including bws)
  ## for simplicity, we don't allow for inconsistent
  ## mixes of named/unnamed arguments
  ## so bws is named or unnamed, and t[xyz]dat collectively either
  ## named or unnamed
  
  if(no.bws){
    tx.str <- ",txdat = txdat"
    ty.str <- ",tydat = tydat"
    tz.str <- ",tzdat = tzdat"
  } else {
    tx.str <- ifelse(txdat.named, ",txdat = txdat","")
    ty.str <- ifelse(tydat.named, ",tydat = tydat","")
    tz.str <- ifelse(tzdat.named, ",tzdat = tzdat","")
    if((!bws.named) && (!txdat.named)){
      tz.str <- ifelse(no.tzdat,"",",tzdat")
    }
  }
  
  eval(parse(text=paste("npplreg(bws = tbw", tx.str, ty.str, tz.str, ",...)")))
}

