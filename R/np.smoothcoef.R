npscoef <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npscoef",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npscoef",bws$call)
        else if (!is.call(bws))
          UseMethod("npscoef",bws)
        else
          UseMethod("npscoef",NULL)
      } else {
        UseMethod("npscoef", NULL)
      }
    } else {
      UseMethod("npscoef", NULL)
    }
  }

npscoef.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf)

    tydat <- model.response(tmf)
    txdat <- tmf[, bws$chromoly[[2]], drop = FALSE]
    if (!(miss.z <- !(length(bws$chromoly) == 3)))
      tzdat <- tmf[, bws$chromoly[[3]], drop = FALSE]
    

    if ((has.eval <- !is.null(newdata))) {
      if (!(has.ey <- succeedWithResponse(tt, newdata)))
        tt <- delete.response(tt)
      
      umf <- emf <- model.frame(tt, data = newdata)

      if (has.ey)
        eydat <- model.response(emf)
      
      exdat <- emf[, bws$chromoly[[2]], drop = FALSE]
      if (!miss.z)
        ezdat <- emf[, bws$chromoly[[3]], drop = FALSE]
    }

    ev <-
      eval(parse(text=paste("npscoef(txdat = txdat, tydat = tydat,",
                   ifelse(miss.z, '','tzdat = tzdat,'),
                   ifelse(has.eval,paste("exdat = exdat,", 
                                         ifelse(has.ey,"eydat = eydat,",""),
                                         ifelse(miss.z,'', 'ezdat = ezdat,')),""),
                   "bws = bws, ...)")))
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    ev
  }

npscoef.call <-
  function(bws, ...) {
    eval(parse(text = paste('npscoef(txdat = eval(bws$call[["xdat"]], environment(bws$call)),',
                 'tydat = eval(bws$call[["ydat"]], environment(bws$call)),',
                 ifelse(!is.null(bws$zdati), 'tzdat = eval(bws$call[["zdat"]], environment(bws$call)),',''),
                 'bws = bws, ...)')))
  }


npscoef.default <- function(bws, txdat, tydat, tzdat, ...) {
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
  
  tbw <- eval(parse(text = paste("npscoefbw(",
                      ifelse(bws.named,                             
                             paste(tx.str, ty.str, tz.str,
                                   "bws = bws, bandwidth.compute = FALSE,"),
                             paste(ifelse(no.bws, "", "bws,"),
                                   tx.str, ty.str, tz.str)),
                      "call = mc, ...",")",sep="")))

  ## tbw <-
  ##  updateBwNameMetadata(nameList = list(ynames = deparse(substitute(tydat))),
  ##                       bws = tbw)

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

  ## because of some ambiguities in how the function might be called
  ## we only drop up to two unnamed arguments, when sometimes dropping three
  ## would be appropriate. 
  ## also, for simplicity, we don't allow for inconsistent
  ## mixes of named/unnamed arguments
  ## so bws is named or unnamed, and t[xyz]dat collectively either
  ## named or unnamed

  tz.str <- ifelse(tzdat.named, ",tzdat = tzdat",
                   ifelse(no.tzdat, "", "tzdat"))

  if(no.bws){
    tx.str <- ",txdat = txdat"
    ty.str <- ",tydat = tydat"
  } else {
    tx.str <- ifelse(txdat.named, ",txdat = txdat","")
    ty.str <- ifelse(tydat.named, ",tydat = tydat","")
    if((!bws.named) && (!txdat.named)){
      ty.str <- ifelse(tydat.named, ",tydat = tydat",
                       ifelse(no.tydat, "", "tydat"))
    }
  }
  
  eval(parse(text=paste("npscoef(bws = tbw", tx.str, ty.str, tz.str, ",...)")))

}

npscoef.scbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           tzdat = NULL,
           exdat,
           eydat,
           ezdat,
           residuals = FALSE,
           errors = TRUE,
           iterate = TRUE,
           maxiter = 100,
           tol = .Machine$double.eps,
           leave.one.out = FALSE,
           betas = FALSE, ...){
    
    miss.z <- missing(tzdat)
    
    miss.ex = missing(exdat)
    miss.ey = missing(eydat)
    

    ## if miss.ex then if !miss.ey then ey and tx must match, to get
    ## oos errors alternatively if miss.ey you get is errors if
    ## !miss.ex then if !miss.ey then ey and ex must match, to get
    ## oos errors alternatively if miss.ey you get NO errors since we
    ## don't evaluate on the training data
    
    txdat <- toFrame(txdat)

    if (!(is.vector(tydat) | is.factor(tydat)))
      stop("'tydat' must be a vector or a factor")

    if (!miss.z)
      tzdat <- toFrame(tzdat)


    if (!miss.ex){
      exdat <- toFrame(exdat)

      if (!miss.z) 
        ezdat <- toFrame(ezdat)
      
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
      
      if (!miss.ey){
        if (dim(exdat)[1] != length(eydat))
          stop("number of evaluation data 'exdat' and dependent data 'eydat' do not match")
      }
      
    } else if(!miss.ey) {
      if (dim(txdat)[1] != length(eydat))
        stop("number of training data 'txdat' and dependent data 'eydat' do not match")
    }

    if(iterate && !is.null(bws$bw.fitted) && !miss.ex){
      warning("iteration is not supported for out of sample evaluations; using overall bandwidths")
      iterate = FALSE
    }

    ## catch and destroy NA's
    goodrows = 1:dim(txdat)[1]
    rows.omit =
      eval(parse(text = paste("attr(na.omit(data.frame(txdat, tydat",
                   ifelse(miss.z,'',',tzdat'),')), "na.action")')))
    
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Training data has no rows without NAs")

    txdat = txdat[goodrows,,drop = FALSE]
    tydat = tydat[goodrows]
    if (!miss.z)
      tzdat <- tzdat[goodrows,, drop = FALSE]


    if (!miss.ex){
      goodrows = 1:dim(exdat)[1]
      rows.omit = eval(parse(text=paste('attr(na.omit(data.frame(exdat',
                               ifelse(miss.ey,"",",eydat"),
                               ifelse(miss.z, "",",ezdat"),
                               ')), "na.action")')))

      goodrows[rows.omit] = 0

      exdat = exdat[goodrows,,drop = FALSE]
      if (!miss.ey)
        eydat = eydat[goodrows]
      if (!miss.z)
        ezdat <- ezdat[goodrows,, drop = FALSE]


      if (all(goodrows==0))
        stop("Evaluation data has no rows without NAs")
    }

    ## convert tydat, eydat to numeric, from a factor with levels from the y-data
    ## used during bandwidth selection.

    if (is.factor(tydat)){
      tydat <- adjustLevels(as.data.frame(tydat), bws$ydati)[,1]
      tydat <- (bws$ydati$all.dlev[[1]])[as.integer(tydat)]
    }
    else
      tydat <- as.double(tydat)


    if (miss.ey)
      eydat <- double()
    else {
      if (is.factor(eydat)){
        eydat <- adjustLevels(as.data.frame(eydat), bws$ydati)[,1]
        eydat <- (bws$ydati$all.dlev[[1]])[as.integer(eydat)]
      }
      else
        eydat <- as.double(eydat)
    }

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati)
    if (!miss.z)
      tzdat <- adjustLevels(tzdat, bws$zdati)
      
    if (!miss.ex){
      exdat <- adjustLevels(exdat, bws$xdati)
      if (!miss.z)
        ezdat <- adjustLevels(ezdat, bws$zdati)
    }

    ## grab the evaluation data before it is converted to numeric
    if(miss.ex){
      teval <- txdat
      if (!miss.z)
        teval <- list(exdat = txdat, ezdat = tzdat)
    } else {
      teval <- exdat
      if (!miss.z)
        teval <- list(exdat = exdat, ezdat = ezdat)
    }


    ## put the unordered, ordered, and continuous data in their own objects
    ## data that is not a factor is continuous.
    
    txdat <- toMatrix(txdat)

    if (!miss.ex){
      exdat <- toMatrix(exdat)
    }

    ## from this point on txdat and exdat have been recast as matrices

    ## construct 'W' matrix
    W.train <- W <- as.matrix(data.frame(1,txdat))
    yW <- as.matrix(data.frame(tydat,1,txdat))

    if (miss.z){
      tzdat <- txdat
      if (!miss.ex)
        ezdat <- exdat
    }
      
    ## need to conserve + propagate desired bandwidth properties

    tyw <- eval(parse(text=paste("npksum(txdat = tzdat, tydat = tydat, weights = W,",
                    ifelse(miss.ex, "", "exdat = ezdat,"),
                    "bws = bws)$ksum")))

    tww <- eval(parse(text=paste("npksum(txdat = tzdat, tydat = W, weights = W,",
                    ifelse(miss.ex, "", "exdat = ezdat,"),
                    "bws = bws)$ksum")))

# jracine - this ought to work but barfs on evaluation data. Reverting to old code.
#    
#    tww <- eval(parse(text=paste("npksum(txdat = tzdat, tydat = yW, weights = yW,",
#                    ifelse(miss.ex, "", "exdat = ezdat,"),
#                    "bws = bws)$ksum")))

#    tyw <- tww[-1,1,]
#    tww <- tww[-1,-1,]

#    stop(cat("dim tww = ", dim(tww),"dim tyw = ", dim(tyw)))

    tnrow <- nrow(txdat)
    enrow <- ifelse(miss.ex, nrow(txdat), nrow(exdat))

    if (!miss.ex)
      W <- as.matrix(data.frame(1,exdat))

## ridging jracine Jan 28 2009

    maxPenalty <- sqrt(.Machine$double.xmax)
    coef.mat <- matrix(maxPenalty,ncol(W),enrow)
    epsilon <- 1.0/enrow
    ridge <- double(enrow)
    doridge <- !logical(enrow)

    nc <- ncol(tww[,,1])

    
    ridger <- function(i) {
      doridge[i] <<- FALSE
      ridge.val <- ridge[i]*tyw[,i][1]/
        (ifelse(tww[,,i][1,1]>=0, 1, -1)*max(.Machine$double.eps,abs(tww[,,i][1,1])))
      tryCatch(solve(tww[,,i]+diag(rep(ridge[i],nc)),
                     tyw[,i]+c(ridge.val,rep(0,nc-1))),
               error = function(e){
                 ridge[i] <<- ridge[i]+epsilon
                 doridge[i] <<- TRUE
                 return(rep(maxPenalty,nc))
               })
    }

    while(any(doridge)){
      iloo <- (1:enrow)[doridge]
      coef.mat[,iloo] <- sapply(iloo, ridger)
    }

##    coef.mat <- sapply(1:enrow, function(i) { solve(tww[,,i], tyw[,i]) })

    if (do.iterate <- (iterate && !is.null(bws$bw.fitted) && miss.ex)){
      resid <- tydat - sapply(1:enrow, function(i) { W[i,, drop = FALSE] %*% coef.mat[,i] })
      
      i = 0
      max.err <- .Machine$double.xmax
      aydat <- abs(tydat) + .Machine$double.eps

      n.part <- (ncol(txdat)+1)
 
      while((max.err > tol) & ((i <- i + 1) <= maxiter)){
        resid.old <- resid
        for(j in 1:n.part){
          ## estimate partial residuals
          partial <- W[,j] * coef.mat[j,] + resid
          
          ## use to calculate new beta implicitly

          ksum1 <- npksum(txdat = tzdat,
                          tydat = partial * W[,j],
                          bws = bws, leave.one.out = leave.one.out)$ksum
          ksum2 <- npksum(txdat = tzdat,
                          tydat = W[,j]^2,
                          bws = bws, leave.one.out = leave.one.out)$ksum

          coef.mat[j,] <- ksum1/
            (ifelse(ksum2>=0, 1, -1)*max(.Machine$double.eps,abs(ksum2)))

          ## estimate new full residuals 
          resid <- partial - W[,j] * coef.mat[j,]
          ## repeat for consistency ?
        }
        max.err <- max(abs(resid.old - resid)/aydat)
      }
      if (max.err > tol)
        warning(paste("backfit iterations did not converge. max err= ", max.err,", tol= ", tol,", maxiter= ", maxiter, sep=''))
      mean <- tydat - resid
    } else {
      mean <- sapply(1:enrow, function(i) { W[i,, drop = FALSE] %*% coef.mat[,i] })
    }

    if (!miss.ey) {
      RSQ = RSQfunc(eydat, mean)
      MSE = MSEfunc(eydat, mean)
      MAE = MAEfunc(eydat, mean)
      MAPE = MAPEfunc(eydat, mean)
      CORR = CORRfunc(eydat, mean)
      SIGN = SIGNfunc(eydat, mean)
    } else if(miss.ex) {
      RSQ = RSQfunc(tydat, mean)
      MSE = MSEfunc(tydat, mean)
      MAE = MAEfunc(tydat, mean)
      MAPE = MAPEfunc(tydat, mean)
      CORR = CORRfunc(tydat, mean)
      SIGN = SIGNfunc(tydat, mean)
    }
      
    if(errors | (residuals & miss.ex)){
      tyw <- npksum(txdat = tzdat, tydat = tydat, weights = W.train, bws = bws)$ksum
      tm <- npksum(txdat = tzdat, tydat = W.train, weights = W.train, bws = bws)$ksum

      mean.fit <- rep(maxPenalty,nrow(txdat))
      epsilon <- 1.0/nrow(txdat)
      ridge.tm <- double(nrow(txdat))
      doridge <- !logical(nrow(txdat))

      nc <- ncol(tm[,,1])

      ridger <- function(i) {
        doridge[i] <<- FALSE
        ridge.val <- ridge.tm[i]*tyw[,i][1]/
          (ifelse(tm[,,i][1,1]>=0, 1, -1)*max(.Machine$double.eps,abs(tm[,,i][1,1])))
        W.train[i,, drop = FALSE] %*% tryCatch(solve(tm[,,i]+diag(rep(ridge.tm[i],nc)),
                      tyw[,i]+c(ridge.val,rep(0,nc-1))),
                      error = function(e){
                        ridge.tm[i] <<- ridge.tm[i]+epsilon
                        doridge[i] <<- TRUE
                        return(rep(maxPenalty,nc))
                      })
      }

      while(any(doridge)){
        ii <- (1:nrow(txdat))[doridge]
        mean.fit[ii] <- sapply(ii, ridger)
      }

      u2.W <- (resid <- tydat - mean.fit)^2
    }

    if(errors){
      ## kernel^2 integrals
      k <- (int.kernels[switch(bws$ckertype,
                              gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
                              epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
                              uniform = CKER_UNI)+1])^length(bws$bw)
      


      u2.W <- sapply(1:tnrow, function(i) { W.train[i,, drop=FALSE]*u2.W[i] })
      u2.W <- t(u2.W)

      V.hat <- eval(parse(text = paste("npksum(txdat = tzdat, tydat = W.train, weights = u2.W,",
                            ifelse(!miss.ex, "exdat = ezdat,", ""), "bws = bws)$ksum")))

      ## asymptotics rely on positive definite nature of tww (ie. M.eval) and V.hat
      ## so choleski decomposition is used to assure their veracity
      merr <- sqrt(sapply(1:enrow, function(i){
        cm <- chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc))))
        k*(W[i,,drop=FALSE] %*% cm %*% V.hat[,,i] %*% cm %*% t(W[i,,drop=FALSE]))
      }))

    }

    eval(parse(text=paste("smoothcoefficient(bws = bws, eval = teval",
                 ", mean = mean,",
                 ifelse(errors & !do.iterate,"merr = merr,",""),
                 ifelse(betas, "beta = t(coef.mat),",""),
                 ifelse(residuals, "resid = resid,",""),
                 "residuals = residuals, betas = betas,",
                 "ntrain = nrow(txdat), trainiseval = miss.ex,",
                 ifelse(miss.ey && !miss.ex, "",
                        "xtra=c(RSQ,MSE,MAE,MAPE,CORR,SIGN)"),")")))

  }

