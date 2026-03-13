npqcmstest <- function(formula,
                       data = NULL,
                       subset,
                       xdat,
                       ydat,
                       model = stop(paste(sQuote("model")," has not been provided")),
                       tau = 0.5,
                       distribution = c("bootstrap", "asymptotic"),
                       bwydat = c("y","varepsilon"),
                       boot.method=c("iid","wild","wild-rademacher"),
                       boot.num = 399,
                       pivot = TRUE,
                       density.weighted = TRUE,
                       random.seed = 42,
                       ...) {
  .npRmpi_require_active_slave_pool(where = "npqcmstest()")
  if (.npRmpi_autodispatch_active())
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

  pcall = paste(deparse(model$call),collapse="")
  if(length(grep("model = TRUE", pcall)) == 0)
    stop(paste(sQuote("model")," is missing components ",
               sQuote("model"),
               ".\nTo fix this please invoke ",
               sQuote("rq"),
               " with ", sQuote("model=TRUE"),
               ".\nSee help for further info.",
               sep=""))

  if(tau <=0 || tau >=1) stop("tau must lie in (0,1)")

  if(boot.num < 9) stop("number of bootstrap replications must be >= 9")

  ## checking for consistent interface usage
  miss.xy = c(missing(xdat),missing(ydat))
  miss.f = missing(formula)
    
  if (any(miss.xy) && !all(miss.xy))
    stop("one of, but not both, xdat and ydat was specified")
  else if(all(miss.xy) & miss.f)
    stop("xdat, and ydat, are missing, and no formula is specified.")
  else if(all(miss.xy) & !miss.f){
    mf.args <- list(formula = formula, data = data, na.action = na.omit)
    if (!missing(subset))
      mf.args$subset <- subset
    mf <- do.call(model.frame, mf.args)
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    na.index <- unclass(attr(xdat,"na.action"))
  } else if(!miss.f){
    stop(paste("A formula was specified along with xdat and ydat.\n",
               "Please see the documentation on proper interface usage."))
  } else {
    xdat = toFrame(xdat)
    
    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]

    na.index <- which(!keep.rows)
  }

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)


  distribution = match.arg(distribution)
  boot.method = match.arg(boot.method)
  bwydat = match.arg(bwydat)  

  qresidual <- function(resid, tau) {
    n.obs <- length(resid)
    out <- rep.int(-tau, n.obs)
    nonmissing <- !is.na(resid)
    out[nonmissing & resid <= 0] <- 1 - tau
    out[!nonmissing] <- NA_real_
    out
  }

  ## Here we go...

  model.resid <- residuals(model, type = "response")

  n = length(model.resid)

  ## ydat is model's residuals, xdat all regressors with types

  .np_progress_note("Computing bandwidths")

  ## What are the optimal bandwidths? We could proceed to use those
  ## for the conditional expectation with raw y, or along the lines of
  ## Zheng (1998) does GCV where the dep var is varepsilon...

  if(bwydat == "y") {
    bw <- .np_progress_with_legacy_suppressed(npregbw(xdat=xdat, ydat=model$y, ...))
  } else if(bwydat == "varepsilon"){
    varepsilon <- qresidual(model.resid, tau)
    bw <- .np_progress_with_legacy_suppressed(npregbw(xdat=xdat, ydat=varepsilon, ...))
  }

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
                   bandwidth.divide = FALSE,...)$ksum/(n*prodh)

  if(min(fhat) == 0)
  stop(paste(sep="","\nAttempt to divide by zero density.",
             "\nYou can try re-running the test with `density.weighted=TRUE'\n"))

  In <- function(xdat, model.resid, bw) {
    
    ## n is the number of observations

    n <- length(model.resid)

    ## Compute In (equation 2.10, Hsiao/Li/racine 2005)
    ## Residuals in cms test replaced with varepsilon
    
    varepsilon <- qresidual(model.resid, tau)

    return( sum(varepsilon*npksum(txdat=xdat,
                                  tydat=varepsilon,
                                  bws=bw$bw,
                                  leave.one.out=TRUE,
                                  bandwidth.divide=TRUE,...)$ksum/fhat)/n^2 )
  }

  Omega.hat <- function(xdat, model.resid, bw) {
  
    ## Variance of In (equation 2.11, Hsiao/Li/racine 2005)
    n <- length(model.resid)

    ## Residuals in cms test replaced with varepsilon
    
    varepsilon <- qresidual(model.resid, tau)

    return( 2*prod(bw$bw[bw$icon])*
           sum(varepsilon^2*
               npksum(txdat=xdat,
                      tydat=varepsilon^2,
                      bws=bw$bw,
                      leave.one.out=TRUE,
                      kernel.pow=2,
                      bandwidth.divide=TRUE,...)$ksum/fhat^2)/n^2 )
  }

  Jn <- function(xdat, model.resid, bw) {
    ## Compute the statistic, supposed to be N(0,1) asymptotically
    n <- length(model.resid)
    n*sqrt(prodh)*In(xdat, model.resid, bw)/sqrt(Omega.hat(xdat, model.resid, bw))
  }


  ## Now conduct a wild bootstrap.. yhat is the fitted model, and we have
  ## rq.resid above... these are external in scope to boot.wild

  yhat <- fitted(model)

  ## data is y,xdat for the rq model...

  draw.wild.mult <- function(n.obs, a, b, p.a) {
    u <- stats::runif(n.obs)
    mult <- rep.int(b, n.obs)
    mult[u <= p.a] <- a
    mult
  }

  boot.wild <- function(model.resid) {

    a <- -0.6180339887499 # (1-sqrt(5))/2
    P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))
    b <- 1.6180339887499 #(1+sqrt(5))/2

    ## Use the wild bootstrap to get a bootstrap vector for y under
    ## the null that the model is correct. Alternatively, we could
    ## pairwise resample Z={y,xdat}. Since the errors will likely be
    ## non-zero, first render mean zero, apply the wild bootstrap,
    ## then add back in the mean

    resid.mean <- mean(model.resid)
    y.star <- yhat + (model.resid - resid.mean) *
      draw.wild.mult(length(model.resid), a, b, P.a) +
      resid.mean

    suppressWarnings(resid <- residuals(rq(y.star~ model$x - 1, tau=tau), type = "response"))

    return(if (pivot) Jn(xdat, resid, bw)
           else In(xdat, resid, bw))
  }

  boot.wild.rademacher <- function(model.resid) {

    a <- -1
    P.a <- 0.5
    b <- 1

    ## Use the wild bootstrap to get a bootstrap vector for y under
    ## the null that the model is correct, using Rademacher variables

    resid.mean <- mean(model.resid)
    y.star <- yhat + (model.resid - resid.mean) *
      draw.wild.mult(length(model.resid), a, b, P.a) +
      resid.mean

    suppressWarnings(resid <- residuals(rq(y.star~ model$x - 1, tau=tau), type = "response"))

    return(if (pivot) Jn(xdat, resid, bw)
           else In(xdat, resid, bw))
  }

  boot.iid <- function(model.resid) {

    ## Simple iid resampling

    y.star <- yhat + model.resid[sample.int(length(model.resid), replace = TRUE)]

    suppressWarnings(resid <- residuals(rq(y.star~ model$x - 1, tau=tau), type = "response"))

    return(if (pivot) Jn(xdat, resid, bw)
           else In(xdat, resid, bw))
  }

  if(distribution == "bootstrap"){
    Sn.bootstrap <- numeric(boot.num)
    progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")
    for (ii in seq_len(boot.num)) {
      if(boot.method == "iid"){
        Sn.bootstrap[ii] <- boot.iid(model.resid)
      } else if(boot.method == "wild"){
        Sn.bootstrap[ii] <- boot.wild(model.resid)
      } else if(boot.method == "wild-rademacher"){
        Sn.bootstrap[ii] <- boot.wild.rademacher(model.resid)
      }
      progress <- .np_progress_step(progress, done = ii)
    }
    progress <- .np_progress_end(progress)
    Sn.bootstrap <- sort(Sn.bootstrap)
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

    tJn$P <- mean(Sn.bootstrap > Sn)

    
  }
  
  ## Restore seed

  .np_seed_exit(seed.state)
  
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
