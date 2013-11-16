# This function implements an individual test of significance for both
# discrete (Racine, hart, Li, 2006, ER) and continuous variables
# (Racine, 1997, JBES). It accepts a data frame for explanatory data
# (mixed datatypes allowed), a vector for y for a regression model, an
# npregbw object, and a set of indices for the columns of X for which
# the test is to be run (default = all).

# Note - this conducts _individual_ tests of significance only. It
# uses a wild bootstrap to handle potential heteroskedasticity (though
# it perhaps could be readily modified to resample (y.star, X) pairs
# and perhaps this is desirable).

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
                                 boot.num=399,
                                 boot.method=c("iid","wild","wild-rademacher","pairwise"),
                                 boot.type=c("I","II"),
                                 pivot=TRUE,
                                 index=seq(1,ncol(xdat)),
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

  In <- numeric(length(index))
  P <- numeric(length(index))

  ## Some constants for the wild bootstrap

  a <- -0.6180339887499  # (1-sqrt(5))/2
  b <- 1.6180339887499   # (1+sqrt(5))/2
  P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))

  ## A vector for storing the resampled statistics

  In.vec <- numeric(boot.num)

  ## ii is the counter for successive elements of In and P...

  In.mat = matrix(data = 0, ncol = length(index), nrow = boot.num)

  ii <- 0

  console <- newLineConsole()

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
                       bws=bws,
                       gradients=TRUE,
                       ...)

    In[ii] <- if(!pivot) {
      mean(npreg.out$grad[,i]^2)
    } else {
      ## Temporarily trap NaN XXX
      npreg.out$gerr[is.nan(npreg.out$gerr)] <- .Machine$double.xmax
      mean((npreg.out$grad[,i]/NZD(npreg.out$gerr[,i]))^2)
    }

    if(boot.method != "pairwise") {

      ## We now construct mhat.xi holding constant the variable whose
      ## significance is being tested at its median. First, make a copy
      ## of the data frame xdat
      
      xdat.eval <- xdat
      
      ## Impose the null by evaluating the conditional holding xdat[,i]
      ## constant at its median (numeric) or mode (factor/ordered) using
      ## uocquantile()
      
      xdat.eval[,i] <- uocquantile(xdat[,i], 0.5)
      
      mhat.xi <-  npreg(txdat = xdat,
                        tydat = ydat,
                        exdat=xdat.eval,
                        bws=bws,
                        ...)$mean
      
      ## Recenter the residuals
      
      delta.bar <- mean(ydat-mhat.xi)
      
      ei <- ydat - mhat.xi - delta.bar
      
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
                              bws=bws.original,
                              ...)

        } else {

          bws.boot <- npregbw(xdat = xdat,
                              ydat = ydat.star,
                              bws=bws.original,
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
                            bws=bws,
                            gradients=TRUE,
                            ...)

      } else {

        npreg.boot <- npreg(txdat = xdat,
                            tydat = ydat.star,
                            bws=bws,
                            gradients=TRUE,
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
          boot.type,
          boot.num)

}

npsigtest.default <- function(bws, xdat, ydat, ...){
  sc.names <- names(sys.call())

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

  mc <- match.call()

  tx.str <- ifelse(xdat.named, "xdat = xdat,",
                   ifelse(no.xdat, "", "xdat,"))
  ty.str <- ifelse(ydat.named, "ydat = ydat,",
                   ifelse(no.ydat, "", "ydat,"))

  tbw <- eval(parse(text = paste("npregbw(",
                      ifelse(bws.named,
                             paste(tx.str, ty.str,
                                   "bws = bws, bandwidth.compute = FALSE,"),
                             paste(ifelse(no.bws, "", "bws,"), tx.str, ty.str)),
                      "call = mc, ...",")",sep="")))

  repair.args <- c("data", "subset", "na.action")

  m.par <- match(repair.args, names(mc), nomatch = 0)
  m.child <- match(repair.args, names(tbw$call), nomatch = 0)

  if(any(m.child > 0)) {
    tbw$call[m.child] <- mc[m.par]
  }

  ## next we repair arguments portion of the call
  m.bws.par <- match(c("bws","xdat","ydat"), names(mc), nomatch = 0)
  m.bws.child <- match(c("bws","xdat","ydat"), as.character(tbw$call), nomatch = 0)
  m.bws.union <- (m.bws.par > 0) & (m.bws.child > 0)

  tbw$call[m.bws.child[m.bws.union]] <- mc[m.bws.par[m.bws.union]]

  environment(tbw$call) <- parent.frame()

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

