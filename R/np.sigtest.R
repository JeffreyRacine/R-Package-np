# This function implements an individual test of significance for both
# discrete (Racine, hart, Li, 2006, ER) and continuous variables
# (Racine, 1997, JBES). It accepts a data frame for explanatory data
# (mixed datatypes allowed), a vector for y for a regression model, an
# npregbw object, and a set of indices for the columns of X for which
# the test is to be run (default = all).

npsigtest <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$xdat))
          UseMethod("npsigtest",bws$formula)
        else if (!is.null(bws$call) && is.null(args$xdat) && (class(bws) != "npregression"))
          UseMethod("npsigtest",bws$call)
        else if (!is.call(bws))
          UseMethod("npsigtest",bws)
        else
          UseMethod("npsigtest",NULL)
      } else {
        UseMethod("npsigtest", NULL)
      }
    } else {
      UseMethod("npsigtest", NULL)
    }
  }

npsigtest.formula <-
  function(bws, data = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]

    ev <- npsigtest(xdat = xdat, ydat = ydat, bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    return(ev)
  }

npsigtest.call <-
  function(bws, ...) {
    ev <- npsigtest(xdat = eval(bws$call[["xdat"]], environment(bws$call)),
                    ydat = eval(bws$call[["ydat"]], environment(bws$call)),
                    bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npsigtest.npregression <-
  function(bws, ...){
    ev <- npsigtest(bws$bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npsigtest.rbandwidth <- function(bws,
                                 xdat = stop("data xdat missing"),
                                 ydat = stop("data ydat missing"),
                                 boot.num = 399,
                                 boot.method = c("iid","wild","wild-rademacher","pairwise"),
                                 boot.type = c("I","II"),
                                 pivot = TRUE,
                                 joint = FALSE,
                                 index = seq(1,ncol(xdat)),
                                 random.seed = 42,
                                 ...) {

  xdat <- toFrame(xdat)

  if(boot.num < 9) stop("number of bootstrap replications must be >= 9")

  ## catch and destroy NA's
  goodrows <- 1:dim(xdat)[1]
  rows.omit <- attr(na.omit(data.frame(xdat,ydat)), "na.action")
  goodrows[rows.omit] <- 0

  if (all(goodrows==0))
    stop("Data has no rows without NAs")

  xdat <- xdat[goodrows,,drop = FALSE]
  ydat <- ydat[goodrows]

  if (is.factor(ydat))
    stop("dependent variable must be continuous.")

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  boot.type <- match.arg(boot.type)
  boot.method <- match.arg(boot.method)

  if(boot.type=="II") {
    ## Store a copy of the bandwidths passed in
    bws.original <- bws
  }

  num.obs <- nrow(xdat)

  ## Test for valid entries in index

  if(any(index < 1 || index>NCOL(xdat))) stop(paste("invalid index provided: index entries must lie between 1 and ",NCOL(xdat),sep=""))
  if(length(unique(index))<length(unique)) stop("index contains repeated values (must be unique)")

  if(!joint) {

    In <- numeric(length(index))
    P <- numeric(length(index))

  }

  ## Some constants for the wild bootstrap

  a <- -0.6180339887499  # (1-sqrt(5))/2
  b <- 1.6180339887499   # (1+sqrt(5))/2
  P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))

  ## A vector for storing the resampled statistics

  In.vec <- numeric(boot.num)

  console <- newLineConsole()

  if(joint==TRUE) {

    ## Joint test

    In.mat = matrix(data = 0, ncol = 1, nrow = boot.num)

    if(boot.type=="II") {

      ## Reset bw vals to original as the ith component of bws gets
      ## overwritten when index changes so needs to be set to its
      ## original value

      bws <- bws.original

    }

    ## Note - xdat must be a data frame

    ## Construct In, the average value of the squared derivatives of
    ## the jth element, discrete or continuous

    npreg.out <- npreg(txdat = xdat,
                       tydat = ydat,
                       bws = bws,
                       gradients = TRUE,
                       ...)

    In <- if(!pivot) {
      mean(npreg.out$grad[,index]^2)
    } else {
      ## Temporarily trap NaN XXX
      npreg.out$gerr[is.nan(npreg.out$gerr)] <- .Machine$double.xmax
      mean((npreg.out$grad[,index]/NZD(npreg.out$gerr[,index]))^2)
    }

    if(boot.method != "pairwise") {

      ## Compute scale and mean of unrestricted residuals

      ei.unres <- scale(residuals(npreg(bws=bws)))
      ei.unres.scale <- attr(ei.unres,"scaled:scale")
      ei.unres.center <- attr(ei.unres,"scaled:center")      

      ## We now construct mhat.xi holding constant the variable whose
      ## significance is being tested at its median. First, make a copy
      ## of the data frame xdat
      
      xdat.eval <- xdat
      
      ## Impose the null by evaluating the conditional mean holding
      ## xdat[,i] constant at its median (numeric) or mode
      ## (factor/ordered) using uocquantile()

      for(i in index) xdat.eval[,i] <- uocquantile(xdat[,i], 0.5)
      
      mhat.xi <-  npreg(txdat = xdat,
                        tydat = ydat,
                        exdat = xdat.eval,
                        bws = bws,
                        ...)$mean

      ## Rescale and recenter the residuals under the null to those
      ## under the alternative
      
      ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
      
      ## Recenter the residuals to have mean zero

      ei <- ei - mean(ei)
      
    }
    
    for(i.star in 1:boot.num) {
      
      if(boot.type=="I") {
        msg <- paste("Bootstrap replication ",
                     i.star,
                     "/",
                     boot.num,
                     sep="")
      } else {
        msg <- paste("Bootstrap rep. ",
                     i.star,
                     "/",
                     boot.num,
                     sep="")
      }

      console <- printPush(msg = msg, console)

      if(boot.method == "iid") {

        ydat.star <- mhat.xi + sample(ei, replace=TRUE)

      } else if(boot.method == "wild") {

        ## Conduct a wild bootstrap. We generate a sample for ydat
        ## (ydat.star) drawn from the conditional mean evaluated
        ## holding the variable tested at its median, and add to that
        ## a wild bootstrap draw from the original disturbance vector

        ydat.star <- mhat.xi + ei*ifelse(rbinom(num.obs, 1, P.a) == 1, a, b)

      } else if(boot.method == "wild-rademacher") {

        ## Conduct a wild bootstrap. We generate a sample for ydat
        ## (ydat.star) drawn from the conditional mean evaluated
        ## holding the variable tested at its median, and add to that
        ## a wild bootstrap draw from the original disturbance vector

        ydat.star <- mhat.xi + ei*ifelse(rbinom(num.obs, 1, P.a) == 1, -1, 1)

      } else if(boot.method =="pairwise") {

        ## Leave variable being tested untouched, resample remaining
        ## pairs of y,X thereby breaking any systematic relationship
        ## between variable being tested in y
        boot.index <- sample(1:num.obs,replace=TRUE)
        ydat.star <- ydat[boot.index]
        xdat.star <- xdat[boot.index,]
        for(i in index) xdat.star[,i] <- xdat[,i]

      }

      if(boot.type=="II") {

        ## For Bootstrap II method, starting values are taken from
        ## bandwidths passed in (bws.original). We then conduct
        ## cross-validation for the bootstrap sample and use only the
        ## new bw for variable i along with the original bandwidths
        ## for the remaining variables

        if(boot.method == "pairwise") {

          bws.boot <- npregbw(xdat = xdat.star,
                              ydat = ydat.star,
                              bws = bws.original,
                              ...)

        } else {

          bws.boot <- npregbw(xdat = xdat,
                              ydat = ydat.star,
                              bws = bws.original,
                              ...)

        }

        ## Copy the new cross-validated bandwidth for variable i into
        ## bw.original and use this below.

        bws <- bws.original

        bws$bw[index] <- bws.boot$bw[index]

      }

      if(boot.method == "pairwise") {

        npreg.boot <- npreg(txdat = xdat.star,
                            tydat = ydat.star,
                            bws = bws,
                            gradients = TRUE,
                            ...)

      } else {

        npreg.boot <- npreg(txdat = xdat,
                            tydat = ydat.star,
                            bws = bws,
                            gradients = TRUE,
                            ...)

      }

      In.vec[i.star] <- if(!pivot) {
        mean(npreg.boot$grad[,index]^2)
      } else {
        ## Temporarily trap NaN XXX
        npreg.boot$gerr[is.nan(npreg.boot$gerr)] <- .Machine$double.xmax
        mean((npreg.boot$grad[,index]/NZD(npreg.boot$gerr[,index]))^2)
      }

      console <- printPop(console)
    }

    ## Compute the P-value

    P <- mean(ifelse(In.vec>In,1,0))

    In.mat[,1] = In.vec

  } else {

    ## Individual test

    ## ii is the counter for successive elements of In and P...

    In.mat = matrix(data = 0, ncol = length(index), nrow = boot.num)

    ii <- 0

    for(i in index) {
      
      ## Increment counter...
      
      ii <- ii + 1
      
      if(boot.type=="II") {
        
        ## Reset bw vals to original as the ith component of bws gets
        ## overwritten when index changes so needs to be set to its
        ## original value
        
        bws <- bws.original
        
      }
      
      ## Note - xdat must be a data frame
      
      ## Construct In, the average value of the squared derivatives of
      ## the jth element, discrete or continuous
      
      npreg.out <- npreg(txdat = xdat,
                         tydat = ydat,
                         bws = bws,
                         gradients = TRUE,
                         ...)
      
      In[ii] <- if(!pivot) {
        mean(npreg.out$grad[,i]^2)
      } else {
        ## Temporarily trap NaN XXX
        npreg.out$gerr[is.nan(npreg.out$gerr)] <- .Machine$double.xmax
        mean((npreg.out$grad[,i]/NZD(npreg.out$gerr[,i]))^2)
      }
      
      if(boot.method != "pairwise") {

        ## Compute scale and mean of unrestricted residuals

        ei.unres <- scale(residuals(npreg(bws=bws)))
        ei.unres.scale <- attr(ei.unres,"scaled:scale")
        ei.unres.center <- attr(ei.unres,"scaled:center")      

        ## We now construct mhat.xi holding constant the variable whose
        ## significance is being tested at its median. First, make a copy
        ## of the data frame xdat
        
        xdat.eval <- xdat
        
        ## Impose the null by evaluating the conditional mean holding
        ## xdat[,i] constant at its median (numeric) or mode
        ## (factor/ordered) using uocquantile()
        
        xdat.eval[,i] <- uocquantile(xdat[,i], 0.5)
        
        mhat.xi <-  npreg(txdat = xdat,
                          tydat = ydat,
                          exdat = xdat.eval,
                          bws = bws,
                          ...)$mean
        
        ## Rescale and recenter the residuals under the null to those
        ## under the alternative
        
        ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
        
        ## Recenter the residuals to have mean zero
        
        ei <- ei - mean(ei)
        
      }
      
      for(i.star in 1:boot.num) {
        
        if(boot.type=="I") {
          msg <- paste("Bootstrap replication ",
                       i.star,
                       "/",
                       boot.num,
                       " for variable ",
                       i,
                       " of (",
                       paste(index,collapse=","),
                       ")... ",
                       sep="")
        } else {
          msg <- paste("Bootstrap rep. ",
                       i.star,
                       "/",
                       boot.num,
                       " for variable ",
                       i,
                       " of (",
                       paste(index,collapse=","),
                       ")... ",
                       sep="")
        }
        
        console <- printPush(msg = msg, console)
        
        if(boot.method == "iid") {
          
          ydat.star <- mhat.xi + sample(ei, replace=TRUE)
          
        } else if(boot.method == "wild") {
          
          ## Conduct a wild bootstrap. We generate a sample for ydat
          ## (ydat.star) drawn from the conditional mean evaluated
          ## holding the variable tested at its median, and add to that
          ## a wild bootstrap draw from the original disturbance vector
          
          ydat.star <- mhat.xi + ei*ifelse(rbinom(num.obs, 1, P.a) == 1, a, b)
          
        } else if(boot.method == "wild-rademacher") {
          
          ## Conduct a wild bootstrap. We generate a sample for ydat
          ## (ydat.star) drawn from the conditional mean evaluated
          ## holding the variable tested at its median, and add to that
          ## a wild bootstrap draw from the original disturbance vector
          
          ydat.star <- mhat.xi + ei*ifelse(rbinom(num.obs, 1, P.a) == 1, -1, 1)
          
        } else if(boot.method =="pairwise") {
          
          ## Leave variable being tested untouched, resample remaining
          ## pairs of y,X thereby breaking any systematic relationship
          ## between variable being tested in y
          boot.index <- sample(1:num.obs,replace=TRUE)
          ydat.star <- ydat[boot.index]
          xdat.star <- xdat
          xdat.star[,-i] <- xdat[boot.index,-i]
          
        }
        
        if(boot.type=="II") {
          
          ## For Bootstrap II method, starting values are taken from
          ## bandwidths passed in (bws.original). We then conduct
          ## cross-validation for the bootstrap sample and use only the
          ## new bw for variable i along with the original bandwidths
          ## for the remaining variables
          
          if(boot.method == "pairwise") {
            
            bws.boot <- npregbw(xdat = xdat.star,
                                ydat = ydat.star,
                                bws = bws.original,
                                ...)
            
          } else {
            
            bws.boot <- npregbw(xdat = xdat,
                                ydat = ydat.star,
                                bws = bws.original,
                                ...)
            
          }
          
          ## Copy the new cross-validated bandwidth for variable i into
          ## bw.original and use this below.
          
          bws <- bws.original
          
          bws$bw[i] <- bws.boot$bw[i]
          
        }
        
        if(boot.method == "pairwise") {
          
          npreg.boot <- npreg(txdat = xdat.star,
                              tydat = ydat.star,
                              bws = bws,
                              gradients = TRUE,
                              ...)
          
        } else {
          
          npreg.boot <- npreg(txdat = xdat,
                              tydat = ydat.star,
                              bws = bws,
                              gradients = TRUE,
                              ...)
          
        }
        
        In.vec[i.star] <- if(!pivot) {
          mean(npreg.boot$grad[,i]^2)
        } else {
          ## Temporarily trap NaN XXX
          npreg.boot$gerr[is.nan(npreg.boot$gerr)] <- .Machine$double.xmax
          mean((npreg.boot$grad[,i]/NZD(npreg.boot$gerr[,i]))^2)
        }
        
        console <- printPop(console)

      }
      
      ## Compute the P-value
      
      P[ii] <- mean(ifelse(In.vec>In[ii],1,0))
      
      In.mat[,ii] = In.vec
      
    }
    
  } ## End invididual test

  console <- printPush(msg ="                                                                                ", console)
  console <- printPop(console)
  console <- printClear(console)
  
  ## Return a list containing the statistic and its P-value
  ## bootstrapped In.vec for each variable...

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  sigtest(In=In,
          In.mat,
          P=P,
          bws = bws,
          ixvar = index,
          boot.method,
          pivot,
          joint,
          boot.type,
          boot.num)

}

npsigtest.default <- function(bws, xdat, ydat, ...){
  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat = if it
  ## was, we need to catch that and map it to dat = otherwise the call
  ## is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  xdat.named <- any(sc.names == "xdat")
  ydat.named <- any(sc.names == "ydat")

  no.bws <- missing(bws)
  no.xdat <- missing(xdat)
  no.ydat <- missing(ydat)

  ## if bws was passed in explicitly, do not compute bandwidths

  if(xdat.named)
    xdat <- toFrame(xdat)

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npregbw)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  tbw <- eval.parent(sc.bw)
  
  ## convention: drop 'bws' and up to two unnamed arguments (including bws)
  if(no.bws){
    tx.str <- ",xdat = xdat"
    ty.str <- ",ydat = ydat"
  } else {
    tx.str <- ifelse(xdat.named, ",xdat = xdat","")
    ty.str <- ifelse(ydat.named, ",ydat = ydat","")
    if((!bws.named) && (!xdat.named)){
      ty.str <- ifelse(ydat.named, ",ydat = ydat",
                       ifelse(no.ydat,"",",ydat"))
    }
  }

  ev <- eval(parse(text=paste("npsigtest(bws = tbw", tx.str, ty.str, ",...)")))

  ev$call <- match.call(expand.dots = FALSE)
  environment(ev$call) <- parent.frame()
  return(ev)
}

