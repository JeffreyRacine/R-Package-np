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
           bws = bw[["bw", exact = TRUE]],
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
           bws = bw[["bw", exact = TRUE]],
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

.npRmpi_cms_collective_context <- function() {
  isTRUE(.npRmpi_autodispatch_called_from_bcast())
}

.npRmpi_cms_iid_index_plan <- function(n, boot.num) {
  index <- matrix(NA_integer_, nrow = boot.num, ncol = n)

  for (b in seq_len(boot.num))
    index[b, ] <- sample.int(n, replace = TRUE)

  index
}

.npRmpi_cms_wild_multiplier_plan <- function(n, boot.num, a, b, p.a) {
  mult <- matrix(NA_real_, nrow = boot.num, ncol = n)

  for (ii in seq_len(boot.num)) {
    u <- stats::runif(n)
    row <- rep.int(b, n)
    row[u <= p.a] <- a
    mult[ii, ] <- row
  }

  mult
}

.npRmpi_cms_numeric_chunks <- function(gathered, size) {
  size <- as.integer(size)[1L]
  if (is.na(size) || size < 1L)
    stop("invalid MPI gather size")

  if (is.matrix(gathered)) {
    if (!identical(ncol(gathered), size))
      stop("npcmstest MPI gather returned malformed matrix output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.array(gathered)) {
    dims <- dim(gathered)
    if (length(dims) < 2L || !identical(dims[[length(dims)]], size))
      stop("npcmstest MPI gather returned malformed array output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.list(gathered)) {
    if (!identical(length(gathered), size))
      stop("npcmstest MPI gather returned malformed list output", call. = FALSE)
    return(gathered)
  }

  chunks <- as.list(gathered)
  if (!identical(length(chunks), size))
    stop("npcmstest MPI gather returned malformed atomic output", call. = FALSE)
  chunks
}

.npRmpi_cms_collective_multiplier_bootstrap <- function(plan,
                                                        model.resid,
                                                        yhat,
                                                        model,
                                                        pivot,
                                                        statistic.batch,
                                                        method,
                                                        progress = NULL,
                                                        comm = 1L) {
  boot.num <- nrow(plan)
  size <- mpi.comm.size(comm)
  rank <- mpi.comm.rank(comm)

  local.idx <- seq.int(rank + 1L, boot.num, by = size)
  .npRmpi_bootstrap_transport_trace(
    what = "npcmstest",
    event = "fanout.collective.start",
    fields = list(rank = rank, size = size, B = boot.num, local = length(local.idx), method = method)
  )

  local.Sn <- numeric(length(local.idx))
  chunk.size <- .np_cms_bootstrap_chunk_size(
    n = length(model.resid),
    boot.num = max(1L, length(local.idx)),
    pivot = pivot
  )

  if (length(local.idx)) {
    for (start in seq.int(1L, length(local.idx), by = chunk.size)) {
      stopi <- min(length(local.idx), start + chunk.size - 1L)
      pos <- seq.int(start, stopi)
      residuals.chunk <- matrix(NA_real_, nrow = length(model.resid),
                                ncol = length(pos))

      for (jj in seq_along(pos)) {
        ii <- local.idx[[pos[[jj]]]]
        y.star <- yhat + model.resid * plan[ii, ]
        residuals.chunk[, jj] <-
          if(is.null(model$family)) {
            residuals(glm(y.star~ model$x - 1), type = "response")
          } else {
            residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
          }
      }

      statistic <- .npRmpi_with_local_regression(
        statistic.batch(residuals.chunk)
      )
      local.Sn[pos] <- if (pivot) statistic[["Jn"]] else statistic[["In"]]
    }
  }

  invisible(gc(FALSE))

  payload <- c(as.numeric(local.idx), local.Sn)
  gathered <- mpi.gather.Robj(payload, root = 0L, comm = comm)

  if (rank == 0L) {
    chunks <- .npRmpi_cms_numeric_chunks(gathered, size = size)
    Sn.bootstrap <- numeric(boot.num)

    for (rr in seq_len(size)) {
      vals <- as.numeric(chunks[[rr]])
      if (!length(vals))
        next
      if ((length(vals) %% 2L) != 0L)
        stop("npcmstest MPI gather returned malformed bootstrap chunk", call. = FALSE)
      n.local <- length(vals) / 2L
      idx <- as.integer(vals[seq_len(n.local)])
      Sn.bootstrap[idx] <- vals[n.local + seq_len(n.local)]
    }

    if (!is.null(progress))
      progress <- .np_progress_step(progress, done = boot.num)
    .npRmpi_bootstrap_transport_trace(
      what = "npcmstest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num, method = method)
    )
    mpi.bcast.Robj(Sn.bootstrap, rank = 0L, comm = comm)
    Sn.bootstrap
  } else {
    .npRmpi_bootstrap_transport_trace(
      what = "npcmstest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num, method = method)
    )
    mpi.bcast.Robj(rank = 0L, comm = comm)
  }
}

.npRmpi_cms_collective_iid_bootstrap <- function(plan,
                                                 model.resid,
                                                 yhat,
                                                 model,
                                                 pivot,
                                                 statistic.batch,
                                                 progress = NULL,
                                                 comm = 1L) {
  boot.num <- nrow(plan)
  size <- mpi.comm.size(comm)
  rank <- mpi.comm.rank(comm)

  local.idx <- seq.int(rank + 1L, boot.num, by = size)
  .npRmpi_bootstrap_transport_trace(
    what = "npcmstest",
    event = "fanout.collective.start",
    fields = list(rank = rank, size = size, B = boot.num, local = length(local.idx), method = "iid")
  )

  local.Sn <- numeric(length(local.idx))
  chunk.size <- .np_cms_bootstrap_chunk_size(
    n = length(model.resid),
    boot.num = max(1L, length(local.idx)),
    pivot = pivot
  )

  if (length(local.idx)) {
    for (start in seq.int(1L, length(local.idx), by = chunk.size)) {
      stopi <- min(length(local.idx), start + chunk.size - 1L)
      pos <- seq.int(start, stopi)
      residuals.chunk <- matrix(NA_real_, nrow = length(model.resid),
                                ncol = length(pos))

      for (jj in seq_along(pos)) {
        ii <- local.idx[[pos[[jj]]]]
        y.star <- yhat + model.resid[plan[ii, ]]
        residuals.chunk[, jj] <-
          if(is.null(model$family)) {
            residuals(glm(y.star~ model$x - 1), type = "response")
          } else {
            residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
          }
      }

      statistic <- .npRmpi_with_local_regression(
        statistic.batch(residuals.chunk)
      )
      local.Sn[pos] <- if (pivot) statistic[["Jn"]] else statistic[["In"]]
    }
  }

  invisible(gc(FALSE))

  payload <- c(as.numeric(local.idx), local.Sn)
  gathered <- mpi.gather.Robj(payload, root = 0L, comm = comm)

  if (rank == 0L) {
    chunks <- .npRmpi_cms_numeric_chunks(gathered, size = size)
    Sn.bootstrap <- numeric(boot.num)

    for (rr in seq_len(size)) {
      vals <- as.numeric(chunks[[rr]])
      if (!length(vals))
        next
      if ((length(vals) %% 2L) != 0L)
        stop("npcmstest MPI gather returned malformed bootstrap chunk", call. = FALSE)
      n.local <- length(vals) / 2L
      idx <- as.integer(vals[seq_len(n.local)])
      Sn.bootstrap[idx] <- vals[n.local + seq_len(n.local)]
    }

    if (!is.null(progress))
      progress <- .np_progress_step(progress, done = boot.num)
    .npRmpi_bootstrap_transport_trace(
      what = "npcmstest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num, method = "iid")
    )
    mpi.bcast.Robj(Sn.bootstrap, rank = 0L, comm = comm)
    Sn.bootstrap
  } else {
    .npRmpi_bootstrap_transport_trace(
      what = "npcmstest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num, method = "iid")
    )
    mpi.bcast.Robj(rank = 0L, comm = comm)
  }
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
  .npRmpi_require_active_slave_pool(where = "npcmstest()")
  if (.npRmpi_autodispatch_active())
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))
  
  pcall = paste(deparse(model$call),collapse="")
  if(length(grep("x = (T|TRUE)[ ,)]", pcall)) == 0 || length(grep("y = (T|TRUE)[ ,)]", pcall)) == 0)
    stop(paste(sQuote("model")," is missing components ", sQuote("x"), " and ",
               sQuote("y"), ".\nTo fix this please invoke ", sQuote("lm"),
               " or ", sQuote("glm"),
               " with ", sQuote("x=TRUE"), " and ", sQuote("y=TRUE"),
               ".\nSee help for further info.", sep=""))

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
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)


  distribution = match.arg(distribution)
  boot.method = match.arg(boot.method)

  ## Here we go...

  model.resid <- residuals(model, type = "response")

  n = length(model.resid)

  ## ydat is model's residuals, xdat all regressors with types

  ##  bw <- npregbw(xdat=xdat,ydat=model.resid)

  .np_progress_note("Computing bandwidths")

  bw <- .np_progress_with_legacy_suppressed(npregbw(xdat=xdat, ydat=model$y, ...))
  
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

  statistic.batch <- function(model.resid) {
    .np_cms_statistics_batch(
      xdat = xdat,
      score = model.resid,
      bw = bw,
      fhat = fhat,
      prodh = prodh,
      pivot = pivot,
      kernel.args = list(...)
    )
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
      if(is.null(model$family)) {
        residuals(glm(y.star~ model$x - 1), type = "response")
      } else {
        residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
      }
    
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
      if(is.null(model$family)) {
        residuals(glm(y.star~ model$x - 1), type = "response")
      } else {
        residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
      }
    
    resid
  }

  resid.iid <- function(model.resid) {

    y.star <- yhat + model.resid[sample.int(length(model.resid), replace = TRUE)]
    resid <-
      if(is.null(model$family)) {
        residuals(glm(y.star~ model$x - 1), type = "response")
      } else {
        residuals(glm(y.star~ model$x - 1,family=model$family), type = "response")
      }
    
    resid
  }

  if(distribution == "bootstrap"){
    progress <- .np_progress_begin("Bootstrap replications", total = B, surface = "bootstrap")

    if (.npRmpi_cms_collective_context() && identical(boot.method, "iid")) {
      plan <- .npRmpi_cms_iid_index_plan(length(model.resid), B)
      post.boot.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      Sn.bootstrap <- .npRmpi_cms_collective_iid_bootstrap(
        plan = plan,
        model.resid = model.resid,
        yhat = yhat,
        model = model,
        pivot = pivot,
        statistic.batch = statistic.batch,
        progress = progress
      )
      assign(".Random.seed", post.boot.seed, envir = .GlobalEnv)
    } else if (.npRmpi_cms_collective_context() && identical(boot.method, "wild")) {
      plan <- .npRmpi_cms_wild_multiplier_plan(
        n = length(model.resid),
        boot.num = B,
        a = -0.6180339887499,
        b = 1.6180339887499,
        p.a = 0.72360679774998
      )
      post.boot.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      Sn.bootstrap <- .npRmpi_cms_collective_multiplier_bootstrap(
        plan = plan,
        model.resid = model.resid,
        yhat = yhat,
        model = model,
        pivot = pivot,
        statistic.batch = statistic.batch,
        method = "wild",
        progress = progress
      )
      assign(".Random.seed", post.boot.seed, envir = .GlobalEnv)
    } else if (.npRmpi_cms_collective_context() && identical(boot.method, "wild-rademacher")) {
      plan <- .npRmpi_cms_wild_multiplier_plan(
        n = length(model.resid),
        boot.num = B,
        a = -1,
        b = 1,
        p.a = 0.5
      )
      post.boot.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      Sn.bootstrap <- .npRmpi_cms_collective_multiplier_bootstrap(
        plan = plan,
        model.resid = model.resid,
        yhat = yhat,
        model = model,
        pivot = pivot,
        statistic.batch = statistic.batch,
        method = "wild-rademacher",
        progress = progress
      )
      assign(".Random.seed", post.boot.seed, envir = .GlobalEnv)
    } else {
      Sn.bootstrap <- numeric(B)
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
          progress <- .np_progress_step(progress, done = ii)
        }

        statistic <- statistic.batch(residuals.chunk)
        Sn.bootstrap[idx] <- if (pivot) statistic[["Jn"]] else statistic[["In"]]
      }
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

    tJn$P <- mean(Sn.bootstrap > Sn)

    
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
