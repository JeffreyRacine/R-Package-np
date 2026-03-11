align_singleindex_grad_index <- function(fit, target_index, beta) {
  fit_index <- as.vector(fit$index)
  ord <- match(target_index, fit_index)
  if (anyNA(ord))
    stop("failed to align single-index gradient fit to target index", call. = FALSE)
  as.matrix(fit$grad)[ord, 1L] / as.numeric(beta[1L])
}

explicit_singleindex_grad_index <- function(bw, tx, y, counts_vec) {
  idx <- rep.int(seq_len(nrow(tx)), counts_vec)
  target_index <- as.vector(as.matrix(tx) %*% bw$beta)
  fit <- npindex(
    bws = bw,
    txdat = tx[idx, , drop = FALSE],
    tydat = y[idx],
    exdat = tx,
    gradients = TRUE
  )
  align_singleindex_grad_index(fit = fit, target_index = target_index, beta = bw$beta)
}

test_that("fixed single-index gradient helper matches duplicate-sample refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  helper <- getFromNamespace(".np_inid_boot_from_index", "npRmpi")

  set.seed(20260310)
  n <- 36
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * tx$x1 - tx$x2) + rnorm(n, sd = 0.05)

  counts <- cbind(
    rep(1L, n),
    c(rep(2L, 6), rep(0L, 6), rep(1L, n - 12)),
    c(rep(0L, 3), rep(3L, 3), rep(1L, n - 6)),
    sample.int(3L, n, replace = TRUE)
  )

  normalize_counts <- function(v, n) {
    out <- as.integer(v)
    while (sum(out) < n)
      out[which.min(out)] <- out[which.min(out)] + 1L
    while (sum(out) > n) {
      hit <- which(out > 0L)[1L]
      out[hit] <- out[hit] - 1L
    }
    out
  }

  counts <- apply(counts, 2L, normalize_counts, n = n)
  if (!is.matrix(counts))
    counts <- matrix(counts, nrow = n)

  cfgs <- list(
    list(regtype = "lc", label = "lc"),
    list(regtype = "ll", label = "ll"),
    list(regtype = "lp", basis = "tensor", degree = 2L, label = "lp")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 0.25),
      bandwidth.compute = FALSE,
      regtype = cfg$regtype,
      bwtype = "fixed"
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- cfg$degree
    }
    bw <- do.call(npindexbw, bw.args)

    boot <- helper(
      xdat = tx,
      ydat = y,
      bws = bw,
      B = ncol(counts),
      counts = counts,
      gradients = TRUE
    )

    oracle_t0 <- explicit_singleindex_grad_index(
      bw = bw,
      tx = tx,
      y = y,
      counts_vec = rep(1L, n)
    )
    oracle_t <- vapply(
      seq_len(ncol(counts)),
      function(j) explicit_singleindex_grad_index(
        bw = bw,
        tx = tx,
        y = y,
        counts_vec = counts[, j]
      ),
      numeric(n)
    )

    expect_equal(
      as.vector(boot$t0),
      oracle_t0,
      tolerance = 1e-10,
      info = paste(cfg$label, "t0 parity")
    )
    expect_equal(
      boot$t,
      t(oracle_t),
      tolerance = 1e-8,
      info = paste(cfg$label, "counts parity")
    )
  }
})

test_that("fixed single-index gradient helper counts.drawer matches counts matrix", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  helper <- getFromNamespace(".np_inid_boot_from_index", "npRmpi")

  set.seed(20260310)
  n <- 30
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(tx$x1 + tx$x2) + rnorm(n, sd = 0.06)
  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 1, 0.25),
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwtype = "fixed"
  )

  counts <- cbind(
    rep(1L, n),
    c(rep(2L, 5), rep(0L, 5), rep(1L, n - 10)),
    c(rep(0L, 4), rep(2L, 4), rep(1L, n - 8))
  )
  if (!is.matrix(counts))
    counts <- matrix(counts, nrow = n)

  drawer <- function(start, stop) counts[, start:stop, drop = FALSE]

  boot.counts <- helper(
    xdat = tx,
    ydat = y,
    bws = bw,
    B = ncol(counts),
    counts = counts,
    gradients = TRUE
  )
  boot.drawer <- helper(
    xdat = tx,
    ydat = y,
    bws = bw,
    B = ncol(counts),
    counts.drawer = drawer,
    gradients = TRUE
  )

  expect_equal(boot.drawer$t0, boot.counts$t0, tolerance = 1e-12)
  expect_equal(boot.drawer$t, boot.counts$t, tolerance = 1e-12)
})

test_that("sibandwidth fixed gradient bootstrap now works for helper methods", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260310)
  n <- 30
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(tx$x1 + tx$x2) + rnorm(n, sd = 0.07)
  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 1, 0.25),
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwtype = "fixed"
  )

  for (boot.method in c("inid", "fixed", "geom")) {
    out <- suppressWarnings(plot(
      bw,
      xdat = tx,
      ydat = y,
      plot.behavior = "data",
      perspective = FALSE,
      gradients = TRUE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = boot.method,
      plot.errors.boot.num = 9L
    ))[[1]]

    expect_true(all(is.finite(out$grad)))
    expect_true(all(is.finite(out$glerr)))
    expect_true(all(is.finite(out$gherr)))
  }
})

test_that("single-index nonfixed gradient helper methods still fail fast", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260310)
  n <- 30
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(tx$x1 + tx$x2) + rnorm(n, sd = 0.08)

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      regtype = "ll",
      bwtype = bt
    )

    for (boot.method in c("inid", "fixed", "geom")) {
      expect_error(
        suppressWarnings(plot(
          bw,
          xdat = tx,
          ydat = y,
          plot.behavior = "data",
          perspective = FALSE,
          gradients = TRUE,
          plot.errors.method = "bootstrap",
          plot.errors.boot.method = boot.method,
          plot.errors.boot.num = 9L
        )),
        "single-index helper unavailable|requires helper mode with gradients=FALSE",
        info = paste(bt, boot.method, "gradient helper guard")
      )
    }
  }
})
