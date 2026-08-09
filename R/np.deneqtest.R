## Function that implements the multivariate density equality test
## described in Li, Q., E. Maasoumi, and J.S. Racine (2009), "A
## Nonparametric Test for Equality of Distributions with Mixed
## Categorical and Continuous Data," Journal of Econometrics, Volume
## 148, pp 186-200.

.npdeneq_count_chunk_size <- function(pool.n,
                                      boot.num,
                                      byte.budget = 16 * 1024^2,
                                      progress.cap = 16L) {
  pool.n <- as.double(pool.n)[1L]
  boot.num <- as.integer(boot.num)[1L]
  by.memory <- floor(byte.budget / (8 * pool.n * 6))
  as.integer(max(1, min(boot.num, progress.cap, by.memory)))
}

.npdeneq_count_compression_eligible <- function(bw) {
  if (is.numeric(bw) || is.integer(bw))
    return(TRUE)

  kbw <- tryCatch(kbandwidth(bw), error = function(e) NULL)
  if (is.null(kbw) ||
      !identical(kbw[["type", exact = TRUE]], "fixed"))
    return(FALSE)

  lower <- kbw[["ckerlb", exact = TRUE]]
  upper <- kbw[["ckerub", exact = TRUE]]
  !any(is.finite(lower) | is.finite(upper))
}

.npRmpi_deneq_collective_context <- function() {
  isTRUE(.npRmpi_autodispatch_called_from_bcast())
}

.npRmpi_deneq_bootstrap_index_plan <- function(pool.n, n1, n2, boot.num) {
  x.index <- matrix(NA_integer_, nrow = boot.num, ncol = n1)
  y.index <- matrix(NA_integer_, nrow = boot.num, ncol = n2)

  for (i in seq_len(boot.num)) {
    x.index[i, ] <- sample.int(pool.n, size = n1, replace = TRUE)
    y.index[i, ] <- sample.int(pool.n, size = n2, replace = TRUE)
  }

  list(x = x.index, y = y.index)
}

.npRmpi_deneq_numeric_chunks <- function(gathered, size) {
  size <- as.integer(size)[1L]
  if (is.na(size) || size < 1L)
    stop("invalid MPI gather size")

  if (is.matrix(gathered)) {
    if (!identical(ncol(gathered), size))
      stop("npdeneqtest MPI gather returned malformed matrix output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.array(gathered)) {
    dims <- dim(gathered)
    if (length(dims) < 2L || !identical(dims[[length(dims)]], size))
      stop("npdeneqtest MPI gather returned malformed array output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.list(gathered)) {
    if (!identical(length(gathered), size))
      stop("npdeneqtest MPI gather returned malformed list output", call. = FALSE)
    return(gathered)
  }

  chunks <- as.list(gathered)
  if (!identical(length(chunks), size))
    stop("npdeneqtest MPI gather returned malformed atomic output", call. = FALSE)
  chunks
}

.npRmpi_deneq_collective_bootstrap <- function(plan,
                                               z,
                                               bw.x,
                                               bw.y,
                                               teststat,
                                               teststat.counted.batch = NULL,
                                               compress.bootstrap = FALSE,
                                               progress = NULL,
                                               comm = 1L) {
  boot.num <- nrow(plan$x)
  if (!identical(nrow(plan$y), boot.num))
    stop("npdeneqtest bootstrap index plan is malformed", call. = FALSE)

  size <- mpi.comm.size(comm)
  rank <- mpi.comm.rank(comm)

  local.idx <- seq.int(rank + 1L, boot.num, by = size)
  .npRmpi_bootstrap_transport_trace(
    what = "npdeneqtest",
    event = "fanout.collective.start",
    fields = list(rank = rank, size = size, B = boot.num, local = length(local.idx))
  )

  local.Tn <- numeric(length(local.idx))
  local.In <- numeric(length(local.idx))

  if (compress.bootstrap && length(local.idx)) {
    chunk.size <- .npdeneq_count_chunk_size(
      pool.n = nrow(z),
      boot.num = length(local.idx)
    )
    for (start in seq.int(1L, length(local.idx), by = chunk.size)) {
      stopi <- min(length(local.idx), start + chunk.size - 1L)
      pos <- seq.int(start, stopi)
      x.count <- matrix(0, nrow = nrow(z), ncol = length(pos))
      y.count <- matrix(0, nrow = nrow(z), ncol = length(pos))

      for (jj in seq_along(pos)) {
        i <- local.idx[[pos[[jj]]]]
        x.count[, jj] <- tabulate(plan$x[i, ], nbins = nrow(z))
        y.count[, jj] <- tabulate(plan$y[i, ], nbins = nrow(z))
      }

      output.boot <- .npRmpi_with_local_regression(
        teststat.counted.batch(
          z = z,
          x.count = x.count,
          y.count = y.count,
          bw.x = bw.x,
          bw.y = bw.y,
          n1 = ncol(plan$x),
          n2 = ncol(plan$y)
        )
      )
      local.Tn[pos] <- output.boot[["Tn"]]
      local.In[pos] <- output.boot[["In"]]
    }
  } else {
    for (jj in seq_along(local.idx)) {
      i <- local.idx[[jj]]
      x.bootstrap <- data.frame(z[plan$x[i, ], , drop = FALSE])
      y.bootstrap <- data.frame(z[plan$y[i, ], , drop = FALSE])
      output.boot <- .npRmpi_with_local_regression(
        teststat(x.bootstrap, y.bootstrap, bw.x, bw.y)
      )
      local.Tn[[jj]] <- output.boot$Tn
      local.In[[jj]] <- output.boot$In
    }
  }

  payload <- c(as.numeric(local.idx), local.Tn, local.In)
  gathered <- mpi.gather.Robj(payload, root = 0L, comm = comm)

  if (rank == 0L) {
    chunks <- .npRmpi_deneq_numeric_chunks(gathered, size = size)
    Tn.vector <- numeric(boot.num)
    In.vector <- numeric(boot.num)

    for (rr in seq_len(size)) {
      vals <- as.numeric(chunks[[rr]])
      if (!length(vals))
        next
      if ((length(vals) %% 3L) != 0L)
        stop("npdeneqtest MPI gather returned malformed bootstrap chunk", call. = FALSE)
      n.local <- length(vals) / 3L
      idx <- as.integer(vals[seq_len(n.local)])
      Tn.vector[idx] <- vals[n.local + seq_len(n.local)]
      In.vector[idx] <- vals[2L * n.local + seq_len(n.local)]
    }

    out <- list(Tn = Tn.vector, In = In.vector)
    if (!is.null(progress))
      progress <- .np_progress_step(progress, done = boot.num)
    .npRmpi_bootstrap_transport_trace(
      what = "npdeneqtest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num)
    )
    mpi.bcast.Robj(out, rank = 0L, comm = comm)
    out
  } else {
    .npRmpi_bootstrap_transport_trace(
      what = "npdeneqtest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num)
    )
    mpi.bcast.Robj(rank = 0L, comm = comm)
  }
}

npdeneqtest <- function(x = NULL,
                        y = NULL,
                        bw.x = NULL,
                        bw.y = NULL,
                        boot.num = 399,
                        random.seed = 42,
                        ...) {
  .npRmpi_require_active_slave_pool(where = "npdeneqtest()")
  if (.npRmpi_autodispatch_active())
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

  ## Some testing of input values

  if(is.null(x) || is.null(y)) stop(" you must provide x and y data")
  if(!is.data.frame(x) || !is.data.frame(y)) stop(" x and y must be data frames")
  if(!identical(names(data.frame(x)),names(data.frame(y)))) stop(" data frames x and y must have identical variable names")
  if(boot.num < 9) stop(" number of bootstrap replications must be >= 9")

  .np_progress_note("Computing bandwidths")

  if(is.null(bw.x) || is.null(bw.y)) {
    bw.x <- .np_progress_with_legacy_suppressed(npudensbw(dat=x,...))
    bw.y <- .np_progress_with_legacy_suppressed(npudensbw(dat=y,...))
  }

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)


  ## First, define test statistic function. This will return the
  ## standardized and unstandardized test statistic along with its
  ## estimated variance.

  teststat <- function(x,y,bw.x,bw.y) {

    ## Get n1 and n2, number of rows in x and y

    n1 <- nrow(x)
    n2 <- nrow(y)

    ## First, compute the In statistic

    ksum.1 <- .npksum_power12(txdat=x,
                              bws=bw.x,
                              leave.one.out=TRUE,
                              bandwidth.divide=TRUE)
    sum.1 <- sum(ksum.1$ksum)
    sum2.1 <- sum(ksum.1$ksum.power2)

    ksum.2 <- .npksum_power12(txdat=y,
                              bws=bw.y,
                              leave.one.out=TRUE,
                              bandwidth.divide=TRUE)
    sum.2 <- sum(ksum.2$ksum)
    sum2.2 <- sum(ksum.2$ksum.power2)

    ksum.3 <- .npksum_power12(txdat=x,
                              exdat=y,
                              bws=bw.x,
                              leave.one.out=FALSE,
                              bandwidth.divide=TRUE)
    sum.3 <- sum(ksum.3$ksum)
    sum2.3 <- sum(ksum.3$ksum.power2)

    ## sum.4 and sum.3 are identical...
    
    In <- sum.1/(n1*(n1-1))+sum.2/(n2*(n2-1))-2*sum.3/(n1*n2)

    ## Next, compute sigma^2_n

    ## sum.4 and sum.3 are identical

    sigma2.n<- 2*(sum2.1/(n1^2*(n1-1)^2)+sum2.2/(n2^2*(n2-1)^2)+2*sum2.3/(n1^2*n2^2))

    ## Finally, compute Tn, the standardized statistic

    Tn <- In/sqrt(sigma2.n)
    
    return(list(Tn=Tn,In=In))
    
  } ## End of test statistic

  teststat.counted.batch <- function(z, x.count, y.count,
                                     bw.x, bw.y, n1, n2) {
    x.count <- as.matrix(x.count)
    y.count <- as.matrix(y.count)

    ksum.x <- npksum(
      txdat = z,
      tydat = x.count,
      exdat = z,
      bws = bw.x,
      bandwidth.divide = TRUE
    )[["ksum", exact = TRUE]]
    dim(ksum.x) <- dim(x.count)
    sum.1 <- colSums(x.count *
                     (ksum.x - as.numeric(self.diagonal.x[["ksum", exact = TRUE]])))
    sum.3 <- colSums(y.count * ksum.x)

    ksum.x2 <- npksum(
      txdat = z,
      tydat = x.count,
      exdat = z,
      bws = bw.x,
      kernel.pow = 2,
      bandwidth.divide = TRUE
    )[["ksum", exact = TRUE]]
    dim(ksum.x2) <- dim(x.count)
    sum2.1 <- colSums(x.count *
                      (ksum.x2 - as.numeric(self.diagonal.x[["ksum.power2", exact = TRUE]])))
    sum2.3 <- colSums(y.count * ksum.x2)

    ksum.y <- npksum(
      txdat = z,
      tydat = y.count,
      exdat = z,
      bws = bw.y,
      bandwidth.divide = TRUE
    )[["ksum", exact = TRUE]]
    dim(ksum.y) <- dim(y.count)
    sum.2 <- colSums(y.count *
                     (ksum.y - as.numeric(self.diagonal.y[["ksum", exact = TRUE]])))

    ksum.y2 <- npksum(
      txdat = z,
      tydat = y.count,
      exdat = z,
      bws = bw.y,
      kernel.pow = 2,
      bandwidth.divide = TRUE
    )[["ksum", exact = TRUE]]
    dim(ksum.y2) <- dim(y.count)
    sum2.2 <- colSums(y.count *
                      (ksum.y2 - as.numeric(self.diagonal.y[["ksum.power2", exact = TRUE]])))

    In <- sum.1 / (n1 * (n1 - 1)) +
      sum.2 / (n2 * (n2 - 1)) -
      2 * sum.3 / (n1 * n2)
    sigma2.n <- 2 * (
      sum2.1 / (n1^2 * (n1 - 1)^2) +
      sum2.2 / (n2^2 * (n2 - 1)^2) +
      2 * sum2.3 / (n1^2 * n2^2)
    )

    list(Tn = In / sqrt(sigma2.n), In = In)
  }

  compress.bootstrap <-
    .npdeneq_count_compression_eligible(bw.x) &&
    .npdeneq_count_compression_eligible(bw.y)
  bootstrap.pool <- data.frame(rbind(x, y))
  self.diagonal.x <- self.diagonal.y <- NULL
  if (compress.bootstrap) {
    self.diagonal.x <- .npksum_power12(
      txdat = bootstrap.pool[1L, , drop = FALSE],
      bws = bw.x,
      bandwidth.divide = TRUE
    )
    self.diagonal.y <- .npksum_power12(
      txdat = bootstrap.pool[1L, , drop = FALSE],
      bws = bw.y,
      bandwidth.divide = TRUE
    )
  }

  ## Now write a bootstrap function for the test statistic
  
  teststat.boot <- function(x,y,bw.x,bw.y) {
    n1 <- nrow(x)
    n2 <- nrow(y)
    ## Resample from pooled data
    z <- bootstrap.pool
    x.index <- sample.int(nrow(z), size = n1, replace = TRUE)
    y.index <- sample.int(nrow(z), size = n2, replace = TRUE)
    x.bootstrap <- data.frame(z[x.index, , drop = FALSE])
    y.bootstrap <- data.frame(z[y.index, , drop = FALSE])
    output.boot <- teststat(x.bootstrap, y.bootstrap, bw.x, bw.y)
    return(list(Tn=output.boot$Tn,
                In=output.boot$In))
  }
  
  Tn.vector <- numeric(boot.num)
  In.vector <- numeric(boot.num)

  progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")

  if (.npRmpi_deneq_collective_context()) {
    z <- data.frame(rbind(x, y))
    plan <- .npRmpi_deneq_bootstrap_index_plan(
      pool.n = nrow(z),
      n1 = nrow(x),
      n2 = nrow(y),
      boot.num = boot.num
    )
    post.boot.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    boot.out <- .npRmpi_deneq_collective_bootstrap(
      plan = plan,
      z = z,
      bw.x = bw.x,
      bw.y = bw.y,
      teststat = teststat,
      teststat.counted.batch = teststat.counted.batch,
      compress.bootstrap = compress.bootstrap,
      progress = progress
    )
    Tn.vector <- boot.out$Tn
    In.vector <- boot.out$In
    assign(".Random.seed", post.boot.seed, envir = .GlobalEnv)
  } else if (compress.bootstrap) {
    pool.n <- nrow(bootstrap.pool)
    n1 <- nrow(x)
    n2 <- nrow(y)
    chunk.size <- .npdeneq_count_chunk_size(pool.n, boot.num)

    for (start in seq.int(1L, boot.num, by = chunk.size)) {
      stopi <- min(boot.num, start + chunk.size - 1L)
      idx <- seq.int(start, stopi)
      x.count <- matrix(0, nrow = pool.n, ncol = length(idx))
      y.count <- matrix(0, nrow = pool.n, ncol = length(idx))

      for (jj in seq_along(idx)) {
        x.count[, jj] <- tabulate(
          sample.int(pool.n, size = n1, replace = TRUE),
          nbins = pool.n
        )
        y.count[, jj] <- tabulate(
          sample.int(pool.n, size = n2, replace = TRUE),
          nbins = pool.n
        )
      }

      output.boot <- teststat.counted.batch(
        z = bootstrap.pool,
        x.count = x.count,
        y.count = y.count,
        bw.x = bw.x,
        bw.y = bw.y,
        n1 = n1,
        n2 = n2
      )
      Tn.vector[idx] <- output.boot[["Tn"]]
      In.vector[idx] <- output.boot[["In"]]
      for (i in idx)
        progress <- .np_progress_step(progress, done = i)
    }
  } else {
    for (i in seq_len(boot.num)) {
      output.boot <- teststat.boot(x,y,bw.x,bw.y)
      Tn.vector[i] <- output.boot$Tn
      In.vector[i] <- output.boot$In
      progress <- .np_progress_step(progress, done = i)
    }
  }

  progress <- .np_progress_end(progress)

  ## Compute the test statistic
  
  output <- if (.npRmpi_deneq_collective_context()) {
    .npRmpi_with_local_regression(teststat(x, y, bw.x, bw.y))
  } else {
    teststat(x,y,bw.x,bw.y)
  }
  
  ## Compute empirical P-values - the number of resampled statistics
  ## more extreme than the original statistic
  
  Tn.P <- mean(Tn.vector > output$Tn)
  In.P <- mean(In.vector > output$In)
  
  ## Restore seed

  .np_seed_exit(seed.state, remove_if_absent = TRUE)
  
  deneqtest(Tn=output$Tn,
            In=output$In,
            Tn.bootstrap=Tn.vector,
            In.bootstrap=In.vector,                
            Tn.P=Tn.P,
            In.P=In.P,
            boot.num=boot.num)
  
}
