npcmstest <- function(formula,
                      data = NULL,
                      subset,
                      xdat,
                      ydat,
                      model = stop(paste(sQuote("model")," has not been provided")),
                      distribution = c("bootstrap", "asymptotic"),
                      boot.method=c("iid","wild","wild-rademacher"),
                      boot.num = 399,
                      pivot = TRUE,
                      density.weighted = TRUE,
                      random.seed = 42,
                      ...) {
  
  pcall = paste(deparse(model$call),collapse="")
  if(length(grep("x = TRUE", pcall)) == 0 | length(grep("y = TRUE", pcall)) == 0)
    stop(paste(sQuote("model")," is missing components ", sQuote("x"), " and ",
               sQuote("y"), ".\nTo fix this please invoke ", sQuote("lm"),
               " or ", sQuote("glm"),
               " with ", sQuote("x=TRUE"), " and ", sQuote("y=TRUE"),
               ".\nSee help for further info.", sep=""))

  if(boot.num < 9) stop("number of bootstrap replications must be >= 9")

  ## checking for consistent interface usage
  miss.xy = c(missing(xdat),missing(ydat))
  miss.f = missing(formula)
    
  if (any(miss.xy) && !all(miss.xy))
    stop("one of, but not both, xdat and ydat was specified")
  else if(all(miss.xy) & miss.f)
    stop("xdat, and ydat, are missing, and no formula is specified.")
  else if(all(miss.xy) & !miss.f){
    mf <- eval(parse(text=paste("model.frame(formula = formula, data = data,",
                       ifelse(missing(subset),"","subset = subset,"),
                           "na.action = na.omit)")))
    
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    na.index <- unclass(attr(xdat,"na.action"))
  } else if(!miss.f){
    stop(paste("A formula was specified along with xdat and ydat.\n",
               "Please see the documentation on proper interface usage."))
  } else {
    xdat = toFrame(xdat)

    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]

    na.index = which(goodrows==0)
  }

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  distribution = match.arg(distribution)
  boot.method = match.arg(boot.method)

  ## Here we go...

  model.resid <- residuals(model, type = "response")

  n = length(model.resid)

  ## ydat is model's residuals, xdat all regressors with types

  ##  bw <- npregbw(xdat=xdat,ydat=model.resid)

  console <- newLineConsole()
  console <- printPush("Bandwidth selection", console)

  bw <- npregbw(xdat=xdat, ydat=model$y, ...)

  console <- printPop(console)
  
  ## Now define the Jn test statistic that takes arguments xdat, the
  ## residual vector, the bandwidth object, and the number of bootstrap
  ## replications

  fhat <- 1

  prodh <- if (bw$ncon == 0) 1.0
  else
    prod(bw$bw[bw$icon])

  if (!density.weighted)
    fhat <- npksum(txdat = xdat,
                   bws = bw$bw, leave.one.out = TRUE,
                   bandwidth.divide = TRUE, ...)$ksum/n


  if(min(fhat) == 0)
  stop(paste(sep="","\nAttempt to divide by zero density.",
             "\nYou can try re-running the test with `density.weighted=TRUE'\n"))

  In <- function(xdat, model.resid, bw) {
    
    ## n is the number of observations

    n <- length(model.resid)

    ## Compute In (equation 2.10, Hsiao/Li/racine 2005)

    return( sum(model.resid*npksum(txdat=xdat,
                                   tydat=model.resid,
                                   bws=bw$bw,
                                   leave.one.out=TRUE,
                                   bandwidth.divide=TRUE, ...)$ksum/fhat)/n^2 )
  }

  Omega.hat <- function(xdat, model.resid, bw) {
  
    ## Variance of In (equation 2.11, Hsiao/Li/racine 2005)

    n <- length(model.resid)
    
    return( 2*prodh*
           sum(model.resid^2*
               npksum(txdat=xdat,
                      tydat=model.resid^2,
                      bws=bw$bw,
                      leave.one.out=TRUE,
                      kernel.pow=2,
                      bandwidth.divide=TRUE, ...)$ksum/fhat^2)/n^2 )
  }

  Jn <- function(xdat, model.resid, bw) {
    ## Compute the statistic, supposed to be N(0,1) asymptotically
    n <- length(model.resid)
    n*sqrt(prodh)*In(xdat, model.resid, bw)/sqrt(Omega.hat(xdat, model.resid, bw))
  }


  ## Now conduct a wild bootstrap.. yhat is the fitted model, and we have
  ## ols.resid above... these are external in scope to boot.wild

  yhat <- fitted(model)

  ## data is y,xdat for the OLS model...

  ## jracine March 8, 2006... not using boot() library (problematic I
  ## realized with [indices] hence unnecessary)

  boot.wild <- function(model.resid) {

    a <- -0.6180339887499 # (1-sqrt(5))/2
    P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))
    b <- 1.6180339887499 # (1+sqrt(5))/2

    ## Use the wild bootstrap to get a bootstrap vector for y under the
    ## null that the model is correct. Alternatively, we could pairwise
    ## resample Z={y,xdat}

    ## jracine removed [indices]

    y.star <- yhat + model.resid*ifelse(rbinom(length(model.resid),1,P.a)==1,a,b)
    resid <-
      if(is.null(model$family)) {
        residuals(glm(y.star~ model$x - 1), type = "response")
      } else {
        residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
      }
    
    return(if (pivot) Jn(xdat, resid, bw)
           else In(xdat, resid, bw))
  }

  boot.wild.rademacher <- function(model.resid) {

    a <- -1
    P.a <- 0.5
    b <- 1

    ## Use the wild bootstrap to get a bootstrap vector for y under
    ## the null that the model is correct, using Rademacher variables

    ## jracine removed [indices]

    y.star <- yhat + model.resid*ifelse(rbinom(length(model.resid),1,P.a)==1,a,b)
    resid <-
      if(is.null(model$family)) {
        residuals(glm(y.star~ model$x - 1), type = "response")
      } else {
        residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
      }
    
    return(if (pivot) Jn(xdat, resid, bw)
           else In(xdat, resid, bw))
  }

  boot.iid <- function(model.resid) {

    y.star <- yhat + sample(model.resid,replace=TRUE)
    resid <-
      if(is.null(model$family)) {
        residuals(glm(y.star~ model$x - 1), type = "response")
      } else {
        residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
      }
    
    return(if (pivot) Jn(xdat, resid, bw)
           else In(xdat, resid, bw))
  }

  if(distribution == "bootstrap"){
    Sn.bootstrap <- numeric(boot.num)

    for(ii in 1:boot.num) {
      console <- printPush(paste(sep="", "Bootstrap replication ",
                                 ii, "/", boot.num, "..."), console)
       if(boot.method == "iid"){
        Sn.bootstrap[ii] <- boot.iid(model.resid)
      } else if(boot.method == "wild"){
        Sn.bootstrap[ii] <- boot.wild(model.resid)
      } else if(boot.method == "wild-rademacher"){
        Sn.bootstrap[ii] <- boot.wild.rademacher(model.resid)
      }
      console <- printPop(console)
    }
    Sn.bootstrap <- sort(Sn.bootstrap)
    ##cat("\n")
  }


  ##  Return a list containing the test statistic etc.

  tIn = In(xdat, model.resid, bw)
  to.h = Omega.hat(xdat, model.resid, bw)

  s.d =
    if (pivot) 1.0
    else sqrt(to.h/prodh)/n

  if(distribution == "asymptotic") {
    
    tJn = list(
      Jn = n*sqrt(prodh)*tIn/sqrt(to.h),
      In = tIn,
      Omega.hat = to.h,
      q.90=qnorm(p = .90, sd = s.d),
      q.95=qnorm(p = .95, sd = s.d),
      q.99=qnorm(p = .99, sd = s.d),
      bw = bw,
      Jn.bootstrap = NA,
      In.bootstrap = NA,
      pivot = pivot)

    Sn = if (pivot) tJn$Jn else tIn

    tJn$P <- (1-pnorm(Sn, sd = s.d))

  } else {
    tJn = list(
      Jn = n*sqrt(prodh)*tIn/sqrt(to.h),
      In = tIn,
      Omega.hat = to.h,
      q.90=Sn.bootstrap[ceiling(0.90*boot.num)],
      q.95=Sn.bootstrap[ceiling(0.95*boot.num)],
      q.99=Sn.bootstrap[ceiling(0.99*boot.num)],
      bw=bw,
      Jn.bootstrap = if(pivot) Sn.bootstrap else NA,
      In.bootstrap = if(pivot) NA else Sn.bootstrap,
      pivot = pivot
      )

    Sn = if (pivot) tJn$Jn else tIn

    tJn$P <- mean(ifelse(Sn.bootstrap > Sn, 1, 0))

    
  }
  
  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
  
  cmstest(Jn = tJn$Jn,
          In = tJn$In,
          Omega.hat = tJn$Omega.hat,
          sd = s.d,
          q.90 = tJn$q.90,
          q.95 = tJn$q.95,
          q.99 = tJn$q.99,
          P = tJn$P,
          bws = bw,
          distribution = distribution,
          Jn.bootstrap = tJn$Jn.bootstrap,
          In.bootstrap = tJn$In.bootstrap,
          pivot = pivot,
          model = model,
          boot.method = boot.method,
          boot.num = boot.num,
          na.index = na.index)
}

