# We implement Ichimura's single index model and Klein and Spady's
# single index model using npksum() and the nlm() minimization
# routine in R. These semiparametric models are used to reduce
# dimensionality to a one-dimensional nonparametric estimator, though
# at the potential cost of misspecification.

# Note - klein and spady's estimator requires a binary y with 0/1
# values only - TRISTEN XXX kindly check that if someone selects this
# approach their y is 0/1 and if not return an informative message.

# Note - X must have at least two columns (otherwise foolish to reduce
# dimensionality) XXX - TRISTEN... kindly check that X has at least
# two columns.

# Note also that identification requires that X contain at least one
# continuous variable. XXX - TRISTEN - kindly ensure that this check
# is made and an appropriate error message is passed back.

# Note also that we will use the so-called scale normalization, i.e.,
# that beta_1=1 (no need to estimate) which reduces search by 1
# parameter (this is obviously restricted search subject to beta_1=1).

# Define the index function model... it is a simple local constant
# estimator of y on a linear index X\beta where beta_1 is presumed to
# be 1 by restriction though, at this stage, the user may feed in any
# value they so desire.

# TRISTEN: I would like this function to accept arguments that
# npksum will automatically grab... also a pretty print function
# similar to regression would rock stating `mean in mean, beta in
# beta' goodness of fit and so forth...

npindex <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npindex",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npindex",bws$call)
        else if (!is.call(bws))
          UseMethod("npindex",bws)
        else
          UseMethod("npindex",NULL)
      } else {
        UseMethod("npindex", NULL)
      }
    } else {
      UseMethod("npindex", NULL)
    }
  }

npindex.formula <-
  function(bws,
           data = NULL,
           newdata = NULL,
           ckertype = c("gaussian", "epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),           
           ...){

    ckertype = match.arg(ckertype)

    if(missing(ckerorder))
      ckerorder = 2
    else if (ckertype == "uniform")
      warning("ignoring kernel order specified with uniform kernel type")
    else {
      kord = eval(formals()$ckerorder) 
      if (!any(kord == ckerorder))
        stop("ckerorder must be one of ", paste(kord,collapse=" "))
    }

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

    ev <-
    eval(parse(text=paste("npindex(txdat = txdat, tydat = tydat,",
                 ifelse(has.eval,paste("exdat = exdat,",ifelse(has.ey,"eydat = eydat,","")),""),
                 "bws = bws, ckertype = ckertype, ckerorder = ckerorder,...)")))
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    ev
  }

npindex.call <-
  function(bws, ...) {
    npindex(txdat = eval(bws$call[["xdat"]], environment(bws$call)),
            tydat = eval(bws$call[["ydat"]], environment(bws$call)),
            bws = bws,
            ckertype = bws$ckertype,
            ckerorder = bws$ckerorder,
            ...)
  }

npindex.default <- function(bws,
                            txdat,
                            tydat,
                            ckertype = c("gaussian", "epanechnikov","uniform"), 
                            ckerorder = c(2,4,6,8),           
                            ...){

  ckertype = match.arg(ckertype)

  if(missing(ckerorder))
    ckerorder = 2
  else if (ckertype == "uniform")
    warning("ignoring kernel order specified with uniform kernel type")
  else {
    kord = eval(formals()$ckerorder) 
    if (!any(kord == ckerorder))
      stop("ckerorder must be one of ", paste(kord,collapse=" "))
  }

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

  tbw <- eval(parse(text = paste("npindexbw(",
                      ifelse(bws.named,                             
                             paste(tx.str, ty.str,
                                   "bws = bws, bandwidth.compute = FALSE, ckertype = ckertype, ckerorder = ckerorder,"),
                             paste(ifelse(no.bws, "", "bws,"), tx.str, ty.str)),
                      "call = mc, ...",")",sep="")))

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

  eval(parse(text=paste("npindex(bws = tbw", tx.str, ty.str, ", ckertype = ckertype, ckerorder = ckerorder,...)")))
#  eval(parse(text=paste("npindex(bws = tbw", tx.str, ty.str, ",...)")))  

}

npindex.sibandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           eydat,
           gradients = FALSE,
           residuals = FALSE,
           errors = FALSE,
           boot.num = 399,
           ckertype = c("gaussian", "epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),           
           ...) {

    ckertype = match.arg(ckertype)

    if(missing(ckerorder))
      ckerorder = 2
    else if (ckertype == "uniform")
      warning("ignoring kernel order specified with uniform kernel type")
    else {
      kord = eval(formals()$ckerorder) 
      if (!any(kord == ckerorder))
        stop("ckerorder must be one of ", paste(kord,collapse=" "))
    }

    no.ex = missing(exdat)
    no.ey = missing(eydat)
    
    ## if no.ex then if !no.ey then ey and tx must match, to get
    ## oos errors alternatively if no.ey you get is errors if
    ## !no.ex then if !no.ey then ey and ex must match, to get
    ## oos errors alternatively if no.ey you get NO errors since we
    ## don't evaluate on the training data
    
    txdat = toFrame(txdat)

    if (!(is.vector(tydat) | is.factor(tydat)))
      stop("'tydat' must be a vector or a factor")

    tydat =
      if (is.factor(tydat))
        as.numeric(levels(tydat))[as.integer(tydat)]
      else
        as.double(tydat)

    if (!no.ex){
      exdat = toFrame(exdat)
      
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
      
      if (!no.ey){
        if (dim(exdat)[1] != length(eydat))
          stop("number of evaluation data 'exdat' and dependent data 'eydat' do not match")
      }
      
    } else if(!no.ey) {
      if (dim(txdat)[1] != length(eydat))
        stop("number of training data 'txdat' and dependent data 'eydat' do not match")
    }

    ## catch and destroy NA's
    goodrows = 1:dim(txdat)[1]
    rows.omit = attr(na.omit(data.frame(txdat,tydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Training data has no rows without NAs")

    txdat = txdat[goodrows,,drop = FALSE]
    tydat = tydat[goodrows]


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
        eydat <- adjustLevels(data.frame(eydat), bws$ydati)[,1]
        eydat <- (bws$ydati$all.dlev[[1]])[as.integer(eydat)]
      }
      else
        eydat <- as.double(eydat)
    }

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati)
      
    if (!no.ex)
      exdat <- adjustLevels(exdat, bws$xdati)

    ## grab the evaluation data before it is converted to numeric
    if(no.ex)
      teval <- txdat
    else
      teval <- exdat

    ## put the unordered, ordered, and continuous data in their own objects
    ## data that is not a factor is continuous.
    
    txdat = toMatrix(txdat)

    if (!no.ex){
      exdat = toMatrix(exdat)
    }


    ## from this point on txdat and exdat have been recast as matrices

    ## TRISTEN - how do we cast txdat here so that matrix multiplication
    ## is permitted? I currently use cbind() but data.frame barfs out...

    ## First, create the scalar index (n \times 1 vector)

    index <- txdat %*% bws$beta

    if(no.ex) {
      index.eval <- index 
      exdat <- txdat
      eydat <- tydat
    } else {
      index.eval <- exdat %*% bws$beta
    }

    ## Next, if no gradients are requested, use (faster) npksum

    if(gradients==FALSE) {

      index.mean <- npksum(txdat=index,
                           tydat=tydat,
                           exdat=index.eval,
                           bws=bws$bw,
                           ckertype=ckertype,
                           ckerorder=ckerorder)$ksum/
                             npksum(txdat=index,
                                    exdat=index.eval,
                                    bws=bws$bw,
                                    ckertype=ckertype,
                                    ckerorder=ckerorder)$ksum

      if(!no.ex & (no.ey | residuals)){
        ## want to evaluate on training data for in sample errors even
        ## if evaluation x's are different from training but no y's
        ## are specified
        
        index.tmean <- npksum(txdat=index,
                              tydat=tydat,
                              bws=bws$bw,
                              ckertype=ckertype,
                              ckerorder=ckerorder)$ksum/
                                npksum(txdat=index,
                                       bws=bws$bw,
                                       ckertype=ckertype,
                                       ckerorder=ckerorder)$ksum
      }

    } else if(gradients==TRUE) {

      model <- npreg(txdat=index,
                     tydat=tydat,
                     exdat=index.eval,
                     bws=bws$bw,
                     regtype="lc",
                     ckertype=ckertype,
                     ckerorder=ckerorder,
                     gradients=TRUE)

      index.mean <- model$mean
      
      ## index.grad is a matrix, one column for each variable, each
      ## equal to its coefficient beta_i times the first derivative of
      ## the local-constant model
      
      index.grad <- as.matrix(model$grad)%*%t(as.vector(bws$beta))

      if(!no.ex & (no.ey | residuals)){

        ## Want to evaluate on training data for in sample errors even
        ## if evaluation x's are different from training but no y's
        ## are specified. Also, needed for variance-covariance matrix
        ## (uses on ly the training data)

        model <- npreg(txdat=index,
                       tydat=tydat,
                       bws=bws$bw,
                       regtype="lc",
                       ckertype=ckertype,
                       ckerorder=ckerorder,
                       gradients=TRUE)

        index.tmean <- model$mean

        index.tgrad <- as.matrix(model$grad)%*%t(as.vector(bws$beta))        

      }


    }

    if (no.ex) {
      index.tmean <- index.mean
    }

    if (no.ex & gradients) {
      index.tgrad <- index.grad
    }

    ## 5/3/2010, jracine, added vcov methods... thanks to Juan Carlos
    ## Escanciano <jescanci@indiana.edu> for pushing me on this for
    ## the Klein and Spady estimator... use index.tmean, index.tgrad
    ## (training X) - need gradients == TRUE in order for this to
    ## work.

    if(bws$method == "ichimura" & gradients == TRUE) {

      ## First row & column of covariance matrix `Bvcov' are zero due
      ## to identification condition that beta_1=0. Note the n n^{-1}
      ## n in V^{-1}\Sigma V^{-1} and the \sqrt{n} in the
      ## normalization of \hat\beta will cancel.

      q <- ncol(txdat)
      Bvcov <- matrix(0,q,q)
      dimnames(Bvcov) <- list(bws$xnames,bws$xnames)

      ## Use the weight matrix so we can compute all expectations with
      ## only one call to npksum (the kernel arguments x\beta do not
      ## change, only the j for X_{ij} in E(X_{ij}|X_i'\beta)

      W <- txdat[,-1,drop=FALSE]

      tyindex <- npksum(txdat = index,
                        tydat = rep(1,length(tydat)),
                        weights = W,
                        bws = bws$bw,
                        ckertype=ckertype,
                        ckerorder=ckerorder)$ksum

      tindex <- npksum(txdat = index,
                       bws = bws$bw,
                       ckertype=ckertype,
                       ckerorder=ckerorder)$ksum

      ## Need to trap case where k-1=1... ksum will return a 1 D
      ## array, need a 1 x n matrix

      if(length(dim(tyindex))==1) tyindex <- matrix(tyindex,nrow=1,ncol=dim(tyindex))

      ## xmex = X_i-\hat E(X_i|X_i'\beta), dimension k\times n.

      xmex <- sapply(1:length(tydat),function(i){W[i,]-tyindex[,i]/tindex[i]})

      ## Need to trap case where k-1=1..., sapply will return a
      ## vector, need a 1 x n matrix

      if(is.vector(xmex)) xmex <- matrix(xmex,nrow=1,ncol=length(xmex))

      ## g^{(1)}=dg/d\beta, first beta normalized to one so this
      ## simplifies computation (beta's drop out)

      dg.db.sq <- (W*index.tgrad[,1])^2

      dg.db.sq.xmex <- sapply(1:length(tydat),function(i){dg.db.sq[i,]*xmex[,i]})      

      uhat <- tydat - index.tmean ## Training y and training mean

      Vinv <- solve(dg.db.sq.xmex%*%t(xmex))

      Sigma <- ((uhat^2)*dg.db.sq.xmex)%*%t(xmex)

      Bvcov[-1,-1] <- Vinv %*% Sigma %*% Vinv
    
      dimnames(Bvcov) <- list(bws$xnames,bws$xnames)      

      ## Now export this in an S3 method...

    } else if(bws$method == "kleinspady" & gradients == TRUE) {

      ## We divide by P(1-P) so test for P=0 or 1...
      
      keep <- which(index.tmean < 1 & index.tmean > 0)
      dg.db <- txdat[,-1,drop=FALSE]*index.tgrad[,1]

      ## First row & column of covariance matrix are zero due to
      ## identification condition that beta_1=0. Note the n^{-1} in
      ## the E and the \sqrt{n} in the normalization of \hat\beta will
      ## cancel.

      q <- ncol(txdat)
      Bvcov <- matrix(0,q,q)
      Bvcov[-1,-1] <- solve(t(dg.db[keep,])%*%(dg.db[keep,]/(index.tmean[keep]*
        (1-index.tmean[keep]))))

      dimnames(Bvcov) <- list(bws$xnames,bws$xnames)      

      ## Now export this in an S3 method...

    }

    ## TRISTEN XXX - for continuous y we want to return the fitted model
    ## along with the measures of goodness of fit RSQ, MSE, and other
    ## measures of goodness of fit. For discrete y (0/1), npconmode()
    ## type measures. But, I don't want to implement asymptotic standard
    ## errors however, these can be readily bootstrapped, so perhaps we
    ## can add a simple bootstrap section that bootstraps beta and the
    ## mean (and gradients and avgderiv)?

    if (gradients){
      boofun = function(data, indices){
        rindex = txdat[indices,] %*% bws$beta
        model = npreg(regtype = 'lc', gradients = TRUE,
          txdat = rindex,
          tydat = tydat[indices],
          exdat = index.eval,
          bws = bws$bw,
          ckertype=ckertype,
          ckerorder=ckerorder,
          )[c('mean','grad')]

        c(model$mean, model$grad, mean(model$grad))
      }

    } else {
      boofun = function(data, indices){
        rindex = txdat[indices,] %*% bws$beta
        npksum(txdat = rindex,
               tydat = tydat[indices],
               exdat = index.eval,
               bws = bws$bw,
               ckertype=ckertype,
               ckerorder=ckerorder)$ksum/
                 npksum(txdat = rindex,
                        exdat = index.eval,
                        bws=bws$bw,
                        ckertype=ckertype,
                        ckerorder=ckerorder)$ksum
      }
    }

    if (errors){

      boot.out = suppressWarnings(boot(data.frame(txdat,tydat), boofun, R = boot.num))

      index.merr = matrix(data = 0, ncol = 1, nrow = length(index.eval))
      index.merr[,] = sqrt(diag(cov(boot.out$t[,1:length(index.eval)])))

      if (gradients) {
        index.gerr = matrix(data = 0, ncol = ncol(txdat), nrow = length(index.eval))
        index.gerr[,] = sqrt(diag(cov(boot.out$t[,(length(index.eval)+1):(2*length(index.eval))])))

        for (i in ncol(txdat))
          index.gerr[,i] = abs(bws$beta[i])*index.gerr[,i]


        index.mgerr = sd(boot.out$t[,2*length(index.eval)+1])
        index.mgerr = abs(bws$beta)*index.mgerr
      }
    }
    ## goodness of fit

    if(bws$method == "ichimura") {
      if (!no.ey) {
        RSQ = RSQfunc(eydat,index.mean)
        MSE = MSEfunc(eydat,index.mean)
        MAE = MAEfunc(eydat,index.mean)
        MAPE = MAPEfunc(eydat,index.mean)
        CORR = CORRfunc(eydat,index.mean)
        SIGN = SIGNfunc(eydat,index.mean)
      } else {
        RSQ = RSQfunc(tydat,index.tmean)
        MSE = MSEfunc(tydat,index.tmean)
        MAE = MAEfunc(tydat,index.tmean)
        MAPE = MAPEfunc(tydat,index.tmean)
        CORR = CORRfunc(tydat,index.tmean)
        SIGN = SIGNfunc(tydat,index.tmean)
      }
      strgof = "xtra=c(RSQ,MSE,MAE,MAPE,CORR,SIGN),"
      strres = ifelse(residuals, "resid = tydat - index.tmean,","")
    } else if(bws$method == "kleinspady") {
      index.pred =
        if (!no.ey) round(index.mean)
        else round(index.tmean)
      
      confusion.matrix =
        table(if (!no.ey) eydat else tydat,
              index.pred, dnn=c("Actual", "Predicted"))
      
      CCR.overall <- sum(diag(confusion.matrix))/sum(confusion.matrix)
      CCR.byoutcome <- diag(confusion.matrix)/rowSums(confusion.matrix)


      fit.mcfadden <- confusion.matrix/sum(confusion.matrix)
      
      fit.mcfadden <- sum(diag(fit.mcfadden)) -
        (sum(fit.mcfadden^2)-sum(diag(fit.mcfadden)^2))


      strgof = "confusion.matrix = confusion.matrix, CCR.overall = CCR.overall,
           CCR.byoutcome =  CCR.byoutcome, fit.mcfadden = fit.mcfadden,"
      strres = ""
    }

    eval(parse(text=paste(
                 "singleindex(bws = bws, index = index.eval, mean = index.mean,",
                 ifelse(errors,"merr = index.merr,",""),
                 ifelse(gradients,"grad = index.grad, mean.grad = colMeans(index.grad), betavcov = Bvcov,",""),
                 ifelse(errors & gradients,"gerr = index.gerr, mean.gerr = index.mgerr,",""),
                 strres,
                 "ntrain = nrow(txdat),", strgof,
                 "trainiseval = no.ex, residuals = residuals, gradients = gradients)")))
  
  }
