.np_cms_validate_model <- function(model, quantile = FALSE) {
  if (is.null(model[["x", exact = TRUE]]) || is.null(model[["y", exact = TRUE]]) ||
      (isTRUE(quantile) && is.null(model[["model", exact = TRUE]])))
    stop(if (quantile)
      "'model' must retain x, y and model components; fit rq with model = TRUE" else
      "'model' must retain x and y components; fit lm or glm with x = TRUE and y = TRUE",
      call. = FALSE)
  invisible(NULL)
}

# Model methods restore excluded rows for user display. Kernel statistics and
# bootstrap refits instead share the model's compact design/response sample.
.np_cms_compact_model <- function(model, xdat, ydat) {
  # Weighted rq retains a weighted x/y pair, but raw residuals/fitted values.
  # Its retained frame owns the raw sample (division by weights loses zero rows).
  if (inherits(model, "rq") && length(model[["weights", exact = TRUE]])) {
    frame <- model[["model", exact = TRUE]]
    model[["y"]] <- model.response(frame)
    model[["x"]] <- model.matrix(model[["terms", exact = TRUE]], frame,
                                contrasts.arg = model[["contrasts", exact = TRUE]])
  }
  n <- nrow(xdat)
  if (NROW(model[["x", exact = TRUE]]) != n ||
      length(model[["y", exact = TRUE]]) != n || length(ydat) != n ||
      !isTRUE(all.equal(as.numeric(model[["y", exact = TRUE]]),
                        as.numeric(ydat), tolerance = 0, check.attributes = FALSE)))
    stop("model and test data must contain the same complete observations in the same order",
         call. = FALSE)
  # Copy-on-modify: never change the caller's object or independently omit
  # elements from a residual vector. Keep the original na.action for reporting.
  model[["na.action"]] <- NULL
  if (length(residuals(model, type = "response")) != n || length(fitted(model)) != n)
    stop("model residuals and fitted values must match its retained compact sample",
         call. = FALSE)
  model
}

# One refit owner for serial and collective bootstrap callers. GLM working
# weights are not observation weights; only its retained prior weights apply.
.np_cms_refit_residuals <- function(model, y.star, tau = NULL) {
  if (!is.null(tau))
    return(residuals(rq(y.star ~ model$x - 1, tau = tau,
                        weights = model[["weights", exact = TRUE]]),
                     type = "response"))
  family <- model[["family", exact = TRUE]]
  weights <- if (is.null(family)) model[["weights", exact = TRUE]] else
    model[["prior.weights", exact = TRUE]]
  residuals(glm(y.star ~ model$x - 1,
                family = if (is.null(family)) gaussian() else family,
                weights = weights, offset = model[["offset", exact = TRUE]]),
            type = "response")
}

.np_cms_bootstrap_chunk_size <- function(n,
                                         boot.num,
                                         pivot,
                                         byte.budget = 16 * 1024^2,
                                         progress.cap = 16L) {
  n <- as.double(n)[1L]
  boot.num <- as.integer(boot.num)[1L]
  matrices <- if (isTRUE(pivot)) 6 else 4
  by.memory <- floor(byte.budget / (8 * n * matrices))
  as.integer(max(1, min(boot.num, progress.cap, by.memory)))
}

.np_cms_statistics_batch <- function(xdat,
                                     score,
                                     bw,
                                     fhat,
                                     prodh,
                                     pivot,
                                     kernel.args = list()) {
  score <- as.matrix(score)
  n <- nrow(score)
  fhat <- as.numeric(fhat)

  ksum <- do.call(
    npksum,
    c(list(txdat = xdat,
           tydat = score,
           bws = if (inherits(bw, "kbandwidth")) bw else bw[["bw", exact = TRUE]],
           leave.one.out = TRUE,
           bandwidth.divide = TRUE),
      kernel.args)
  )[["ksum", exact = TRUE]]
  dim(ksum) <- dim(score)
  In <- colSums(score * ksum / fhat) / n^2

  if (!isTRUE(pivot))
    return(list(In = In))

  score2 <- score^2
  ksum2 <- do.call(
    npksum,
    c(list(txdat = xdat,
           tydat = score2,
           bws = if (inherits(bw, "kbandwidth")) bw else bw[["bw", exact = TRUE]],
           leave.one.out = TRUE,
           kernel.pow = 2,
           bandwidth.divide = TRUE),
      kernel.args)
  )[["ksum", exact = TRUE]]
  dim(ksum2) <- dim(score2)
  Omega.hat <- 2 * prodh * colSums(score2 * ksum2 / fhat^2) / n^2

  list(In = In,
       Omega.hat = Omega.hat,
       Jn = n * sqrt(prodh) * In / sqrt(Omega.hat))
}

npcmstest <- function(formula,
                      data = NULL,
                      subset,
                      xdat,
                      ydat,
                      model = stop(paste(sQuote("model")," has not been provided")),
                      distribution = c("bootstrap", "asymptotic"),
                      boot.method=c("iid","wild","wild-rademacher"),
                      B = 399,
                      pivot = TRUE,
                      density.weighted = TRUE,
                      random.seed = 42,
                      ...) {

  if (...length())
    npRejectLegacyBootstrapCount(names(list(...)), "npcmstest")
  
  .np_cms_validate_model(model, quantile = FALSE)

  if(B < 9) stop("number of bootstrap replications must be >= 9")

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
      mf.args[c("data", "subset")] <- .np_formula_subset_inputs(
        data, substitute(subset), parent.frame())
    mf <- do.call(model.frame, mf.args)
    
    ydat <- model.response(mf)
    xdat <- mf[, .np_formula_term_names(attr(attr(mf, "terms"),"term.labels")), drop = FALSE]

    na.index <- unclass(attr(mf,"na.action"))
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

  model <- .np_cms_compact_model(model, xdat, ydat)

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)


  distribution = match.arg(distribution)
  boot.method = match.arg(boot.method)

  ## Here we go...

  model.resid <- residuals(model, type = "response")

  n = length(model.resid)

  ## ydat is model's residuals, xdat all regressors with types

  ##  bw <- npregbw(xdat=xdat,ydat=model.resid)

  bw <- .np_progress_select_bandwidth_enhanced(
    "Computing bandwidths",
    npregbw(xdat=xdat, ydat=model$y, ...))
  
  # Selection owns its dots. Every contraction uses its resolved kernel,
  # including manual bandwidths, scaling, nondefault kernels and bounds.
  kernel.bw <- kbandwidth(bw)

  ## Now define the Jn test statistic that takes arguments xdat, the
  ## residual vector, the bandwidth object, and the number of bootstrap
  ## replications

  fhat <- 1

  prodh <- if (bw$ncon == 0) 1.0
  else
    prod(bw$bw[bw$icon])

  if (!density.weighted)
    fhat <- npksum(txdat = xdat,
                   bws = kernel.bw, leave.one.out = TRUE,
                   bandwidth.divide = TRUE)$ksum/n


  if(min(fhat) == 0)
  stop(paste(sep="","\nAttempt to divide by zero density.",
             "\nYou can try re-running the test with `density.weighted=TRUE'\n"))

  In <- function(xdat, model.resid, bw) {
    
    ## n is the number of observations

    n <- length(model.resid)

    ## Compute In (equation 2.10, Hsiao/Li/racine 2005)

    return( sum(model.resid*npksum(txdat=xdat,
                                   tydat=model.resid,
                                   bws=kernel.bw,
                                   leave.one.out=TRUE,
                                   bandwidth.divide=TRUE)$ksum/fhat)/n^2 )
  }

  Omega.hat <- function(xdat, model.resid, bw) {
  
    ## Variance of In (equation 2.11, Hsiao/Li/racine 2005)

    n <- length(model.resid)
    
    return( 2*prodh*
           sum(model.resid^2*
               npksum(txdat=xdat,
                      tydat=model.resid^2,
                      bws=kernel.bw,
                      leave.one.out=TRUE,
                      kernel.pow=2,
                      bandwidth.divide=TRUE)$ksum/fhat^2)/n^2 )
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

  draw.wild.mult <- function(n.obs, a, b, p.a) {
    u <- stats::runif(n.obs)
    mult <- rep.int(b, n.obs)
    mult[u <= p.a] <- a
    mult
  }

  resid.wild <- function(model.resid) {

    a <- -0.6180339887499 # (1-sqrt(5))/2
    P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))
    b <- 1.6180339887499 # (1+sqrt(5))/2

    ## Use the wild bootstrap to get a bootstrap vector for y under the
    ## null that the model is correct. Alternatively, we could pairwise
    ## resample Z={y,xdat}

    ## jracine removed [indices]

    y.star <- yhat + model.resid * draw.wild.mult(length(model.resid), a, b, P.a)
    resid <-
      .np_cms_refit_residuals(model, y.star)
    
    resid
  }

  resid.wild.rademacher <- function(model.resid) {

    a <- -1
    P.a <- 0.5
    b <- 1

    ## Use the wild bootstrap to get a bootstrap vector for y under
    ## the null that the model is correct, using Rademacher variables

    ## jracine removed [indices]

    y.star <- yhat + model.resid * draw.wild.mult(length(model.resid), a, b, P.a)
    resid <-
      .np_cms_refit_residuals(model, y.star)
    
    resid
  }

  resid.iid <- function(model.resid) {

    y.star <- yhat + model.resid[sample.int(length(model.resid), replace = TRUE)]
    resid <-
      .np_cms_refit_residuals(model, y.star)
    
    resid
  }

  if(distribution == "bootstrap"){
    Sn.bootstrap <- numeric(B)
    progress <- .np_progress_begin("Bootstrap replications", total = B, surface = "bootstrap")

    chunk.size <- .np_cms_bootstrap_chunk_size(
      n = n,
      boot.num = B,
      pivot = pivot
    )
    for (start in seq.int(1L, B, by = chunk.size)) {
      stopi <- min(B, start + chunk.size - 1L)
      idx <- seq.int(start, stopi)
      residuals.chunk <- matrix(NA_real_, nrow = n, ncol = length(idx))

      for (jj in seq_along(idx)) {
        ii <- idx[[jj]]
        residuals.chunk[, jj] <- if(boot.method == "iid") {
          resid.iid(model.resid)
        } else if(boot.method == "wild") {
          resid.wild(model.resid)
        } else {
          resid.wild.rademacher(model.resid)
        }
        progress <- .np_progress_step(progress)
      }

      statistic <- .np_cms_statistics_batch(
        xdat = xdat,
        score = residuals.chunk,
        bw = kernel.bw,
        fhat = fhat,
        prodh = prodh,
        pivot = pivot,
        kernel.args = list()
      )
     Sn.bootstrap[idx] <- if (pivot) statistic[["Jn"]] else statistic[["In"]]
      progress <- .np_progress_step(progress, done = max(idx))
    }
    progress <- .np_progress_end(progress)
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
      q.90=Sn.bootstrap[ceiling(0.90*B)],
      q.95=Sn.bootstrap[ceiling(0.95*B)],
      q.99=Sn.bootstrap[ceiling(0.99*B)],
      bw=bw,
      Jn.bootstrap = if(pivot) Sn.bootstrap else NA,
      In.bootstrap = if(pivot) NA else Sn.bootstrap,
      pivot = pivot
      )

    Sn = if (pivot) tJn$Jn else tIn

    tJn$P <- .np_bootstrap_upper_tail_pvalue(Sn.bootstrap, Sn)

    
  }
  
  ## Restore seed

  .np_seed_exit(seed.state, remove_if_absent = TRUE)
  
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
          boot.num = B,
          na.index = na.index)
}
