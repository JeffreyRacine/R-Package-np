test_that("npudenshat matches npudens and preserves matrix/apply parity across bwtypes", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudenshat <- getFromNamespace("npudenshat", "npRmpi")

  set.seed(20260310)
  n <- 60
  x <- sort(runif(n))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 21))
  iota <- rep.int(1, n)

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- npudensbw(
      dat = tx,
      bwtype = bwtype,
      bws = if (identical(bwtype, "fixed")) 0.18 else 9,
      bandwidth.compute = FALSE
    )

    fit.in <- npudens(bws = bw, tdat = tx)
    fit.ex <- npudens(bws = bw, tdat = tx, edat = ex)

    H.in <- npudenshat(bws = bw, tdat = tx, output = "matrix")
    H.ex <- npudenshat(bws = bw, tdat = tx, edat = ex, output = "matrix")

    apply.in <- npudenshat(bws = bw, tdat = tx, y = iota, output = "apply")
    apply.ex <- npudenshat(bws = bw, tdat = tx, edat = ex, y = iota, output = "apply")

    expect_s3_class(H.in, "npudenshat")
    expect_equal(as.vector(H.in %*% iota), as.vector(fit.in$dens), tolerance = 1e-10, info = bwtype)
    expect_equal(as.vector(H.ex %*% iota), as.vector(fit.ex$dens), tolerance = 1e-10, info = bwtype)
    expect_equal(as.vector(apply.in), as.vector(H.in %*% iota), tolerance = 1e-12, info = bwtype)
    expect_equal(as.vector(apply.ex), as.vector(H.ex %*% iota), tolerance = 1e-12, info = bwtype)
  }
})

test_that("npudenshat fixed-bandwidth count vectors reproduce resampled npudens fits", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudenshat <- getFromNamespace("npudenshat", "npRmpi")

  set.seed(20260311)
  n <- 45
  x <- sort(runif(n))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 17))
  bw <- npudensbw(
    dat = tx,
    bws = 0.16,
    bwtype = "fixed",
    bandwidth.compute = FALSE
  )
  H <- npudenshat(bws = bw, tdat = tx, edat = ex, output = "matrix")
  counts <- rmultinom(n = 3L, size = n, prob = rep.int(1 / n, n))

  for (b in seq_len(ncol(counts))) {
    idx <- rep.int(seq_len(n), counts[, b])
    fit <- npudens(bws = bw, tdat = tx[idx, , drop = FALSE], edat = ex)
    apply.out <- npudenshat(
      bws = bw,
      tdat = tx,
      edat = ex,
      y = counts[, b],
      output = "apply"
    )

    expect_equal(as.vector(H %*% counts[, b]), as.vector(apply.out), tolerance = 1e-12, info = b)
    expect_equal(as.vector(apply.out), as.vector(fit$dens), tolerance = 1e-10, info = b)
  }
})

test_that("npudenshat fixed apply mode matches matrix RHS multiplication", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudenshat <- getFromNamespace("npudenshat", "npRmpi")

  set.seed(20260311)
  n <- 48
  x <- sort(runif(n))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 19))
  bw <- npudensbw(
    dat = tx,
    bws = 0.16,
    bwtype = "fixed",
    bandwidth.compute = FALSE
  )
  rhs <- cbind(seq_len(n) / n, sin(seq_len(n) / 7))
  H <- npudenshat(bws = bw, tdat = tx, edat = ex, output = "matrix")
  apply.out <- npudenshat(bws = bw, tdat = tx, edat = ex, y = rhs, output = "apply")

  expect_equal(apply.out, H %*% rhs, tolerance = 1e-12)
})

test_that("npudenshat preserves bounded gaussian manual-bandwidth semantics", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudenshat <- getFromNamespace("npudenshat", "npRmpi")

  x <- sort(c(5, 11, 21, 31, 46, 75, 98, 122, 145, 165, 195, 224, 245, 293, 321, 330, 350, 420))
  tx <- data.frame(x = x)
  iota <- rep.int(1, nrow(tx))

  bw <- npudensbw(
    dat = tx,
    bws = 100,
    bandwidth.compute = FALSE,
    ckertype = "gaussian",
    ckerbound = "range"
  )

  fit <- npudens(bws = bw, tdat = tx)
  H <- npudenshat(bws = bw, tdat = tx, output = "matrix")
  apply.out <- npudenshat(bws = bw, tdat = tx, y = iota, output = "apply")

  expect_equal(as.vector(H %*% iota), as.vector(fit$dens), tolerance = 1e-7)
  expect_equal(as.vector(apply.out), as.vector(fit$dens), tolerance = 1e-7)
})
