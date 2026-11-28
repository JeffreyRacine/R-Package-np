test_that("npudisthat matches npudist and preserves matrix/apply parity across bwtypes", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudisthat <- getFromNamespace("npudisthat", "npRmpi")

  set.seed(20260310)
  n <- 60
  x <- sort(runif(n))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 21))
  iota <- rep.int(1, n)

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- npudistbw(
      dat = tx,
      bwtype = bwtype,
      bws = if (identical(bwtype, "fixed")) 0.18 else 9,
      bandwidth.compute = FALSE
    )

    fit.in <- npudist(bws = bw, tdat = tx)
    fit.ex <- npudist(bws = bw, tdat = tx, edat = ex)

    H.in <- npudisthat(bws = bw, tdat = tx, output = "matrix")
    H.ex <- npudisthat(bws = bw, tdat = tx, edat = ex, output = "matrix")

    apply.in <- npudisthat(bws = bw, tdat = tx, y = iota, output = "apply")
    apply.ex <- npudisthat(bws = bw, tdat = tx, edat = ex, y = iota, output = "apply")

    expect_s3_class(H.in, "npudisthat")
    expect_equal(as.vector(H.in %*% iota), as.vector(fit.in$dist), tolerance = 1e-10, info = bwtype)
    expect_equal(as.vector(H.ex %*% iota), as.vector(fit.ex$dist), tolerance = 1e-10, info = bwtype)
    expect_equal(as.vector(apply.in), as.vector(H.in %*% iota), tolerance = 1e-12, info = bwtype)
    expect_equal(as.vector(apply.ex), as.vector(H.ex %*% iota), tolerance = 1e-12, info = bwtype)
  }
})

test_that("npudisthat distinguishes training identity from equal external values", {
  npudisthat <- getFromNamespace("npudisthat", "npRmpi")

  set.seed(20260817)
  tx <- data.frame(x = sort(runif(48L)))
  iota <- rep.int(1.0, nrow(tx))
  bw <- npudistbw(
    dat = tx,
    bwtype = "generalized_nn",
    bws = 7,
    bandwidth.compute = FALSE
  )

  H.training <- npudisthat(bws = bw, tdat = tx)
  H.external <- npudisthat(bws = bw, tdat = tx, edat = tx)
  fit.training <- npudist(bws = bw, tdat = tx)
  fit.external <- npudist(bws = bw, tdat = tx, edat = tx)

  expect_equal(drop(H.training %*% iota), fit.training$dist,
               tolerance = 1e-12)
  expect_equal(drop(H.external %*% iota), fit.external$dist,
               tolerance = 1e-12)
  expect_gt(max(abs(H.training - H.external)), 1e-8)
})

test_that("npudisthat fixed-bandwidth count vectors reproduce resampled npudist fits", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudisthat <- getFromNamespace("npudisthat", "npRmpi")

  set.seed(20260311)
  n <- 45
  x <- sort(runif(n))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 17))
  bw <- npudistbw(
    dat = tx,
    bws = 0.16,
    bwtype = "fixed",
    bandwidth.compute = FALSE
  )
  H <- npudisthat(bws = bw, tdat = tx, edat = ex, output = "matrix")
  counts <- rmultinom(n = 3L, size = n, prob = rep.int(1 / n, n))

  for (b in seq_len(ncol(counts))) {
    idx <- rep.int(seq_len(n), counts[, b])
    fit <- npudist(bws = bw, tdat = tx[idx, , drop = FALSE], edat = ex)
    apply.out <- npudisthat(
      bws = bw,
      tdat = tx,
      edat = ex,
      y = counts[, b],
      output = "apply"
    )

    expect_equal(as.vector(H %*% counts[, b]), as.vector(apply.out), tolerance = 1e-12, info = b)
    expect_equal(as.vector(apply.out), as.vector(fit$dist), tolerance = 1e-10, info = b)
  }
})

test_that("npudisthat fixed apply mode matches matrix RHS multiplication", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudisthat <- getFromNamespace("npudisthat", "npRmpi")

  set.seed(20260311)
  n <- 48
  x <- sort(runif(n))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 19))
  bw <- npudistbw(
    dat = tx,
    bws = 0.16,
    bwtype = "fixed",
    bandwidth.compute = FALSE
  )
  rhs <- cbind(seq_len(n) / n, cos(seq_len(n) / 9))
  H <- npudisthat(bws = bw, tdat = tx, edat = ex, output = "matrix")
  apply.out <- npudisthat(bws = bw, tdat = tx, edat = ex, y = rhs, output = "apply")

  expect_equal(apply.out, H %*% rhs, tolerance = 1e-12)
})

test_that("npudisthat preserves bounded gaussian manual-bandwidth semantics", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npudisthat <- getFromNamespace("npudisthat", "npRmpi")

  x <- sort(c(5, 11, 21, 31, 46, 75, 98, 122, 145, 165, 195, 224, 245, 293, 321, 330, 350, 420))
  tx <- data.frame(x = x)
  iota <- rep.int(1, nrow(tx))

  bw <- npudistbw(
    dat = tx,
    bws = 100,
    bandwidth.compute = FALSE,
    ckertype = "gaussian",
    ckerbound = "range"
  )

  fit <- npudist(bws = bw, tdat = tx)
  H <- npudisthat(bws = bw, tdat = tx, output = "matrix")
  apply.out <- npudisthat(bws = bw, tdat = tx, y = iota, output = "apply")

  expect_equal(as.vector(H %*% iota), as.vector(fit$dist), tolerance = 1e-7)
  expect_equal(as.vector(apply.out), as.vector(fit$dist), tolerance = 1e-7)
})

test_that("multi-response hats reuse kernel weights and the first normalization", {
  .npRmpi_with_local_regression({
    x <- data.frame(x = seq(0, 1, length.out = 48L))
    y <- cbind(score = x$x, score = cos(3 * x$x))
    kernel.calls <- 0L
    operations <- character()
    kernel.weights <- .np_kernel_weights_direct
    exact.apply <- .np_exact_operator_apply
    testthat::local_mocked_bindings(
      .np_kernel_weights_direct = function(...) {
        kernel.calls <<- kernel.calls + 1L
        kernel.weights(...)
      },
      .np_exact_operator_apply = function(..., where) {
        operations <<- c(operations, where)
        exact.apply(..., where = where)
      }
    )

    for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
      bw <- npudistbw(dat = x, bwtype = type,
                      bws = if (type == "fixed") 0.2 else 7,
                      bandwidth.compute = FALSE)
      apply.hat <- function(response)
        npudisthat(bws = bw, tdat = x, y = response, output = "apply")
      reference <- cbind(apply.hat(y[, 1L]), apply.hat(y[, 2L]))
      kernel.calls <- 0L
      operations <- character()
      result <- apply.hat(y)
      expect_identical(unname(result), unname(reference))
      if (type == "fixed") {
        expect_identical(colnames(result), colnames(y))
        expect_identical(kernel.calls, 1L)
      } else {
        expect_identical(operations, paste("npudisthat exact operator",
                                          c("apply", "probe", "apply")))
      }
    }
    bw <- npudensbw(dat = x, bws = 0.2, bandwidth.compute = FALSE)
    apply.density <- function(response)
      npudenshat(bws = bw, tdat = x, y = response, output = "apply")
    reference <- cbind(apply.density(y[, 1L]), apply.density(y[, 2L]))
    kernel.calls <- 0L
    result <- apply.density(y)
    expect_identical(unname(result), unname(reference))
    expect_identical(colnames(result), colnames(y))
    expect_identical(kernel.calls, 1L)
  })
})
