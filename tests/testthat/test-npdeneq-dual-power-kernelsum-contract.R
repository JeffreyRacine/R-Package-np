test_that("internal dual-power kernel sums equal separate scalar calls", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")

  set.seed(20260717)
  x <- data.frame(
    c = rnorm(48),
    u = factor(sample(letters[1:3], 48, TRUE)),
    o = ordered(sample(1:4, 48, TRUE))
  )
  bw <- c(0.6, 0.2, 0.3)

  scalar1 <- local_eval(npksum(
    bws = bw, txdat = x, leave.one.out = TRUE,
    bandwidth.divide = TRUE))
  scalar2 <- local_eval(npksum(
    bws = bw, txdat = x, kernel.pow = 2,
    leave.one.out = TRUE, bandwidth.divide = TRUE))
  dual_sum <- getFromNamespace(".npksum_power12", "npRmpi")
  dual <- local_eval(dual_sum(
    bws = bw, txdat = x, leave.one.out = TRUE,
    bandwidth.divide = TRUE))

  expect_identical(dual$ksum, scalar1$ksum)
  expect_identical(dual$ksum.power2, scalar2$ksum)
  expect_false("ksum.power2" %in% names(scalar1))
})

test_that("dual-power route is exact for nearest-neighbor bandwidths", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)
  dual_sum <- getFromNamespace(".npksum_power12", "npRmpi")
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")

  set.seed(20260718)
  x <- data.frame(x = rnorm(52))
  ex <- data.frame(x = rnorm(31))

  for (type in c("generalized_nn", "adaptive_nn")) {
    args <- list(bws = 8, txdat = x, exdat = ex, bwtype = type,
                 bandwidth.divide = TRUE)
    scalar1 <- local_eval(do.call(npksum, args))
    scalar2 <- local_eval(do.call(npksum, c(args, list(kernel.pow = 2))))
    dual <- local_eval(do.call(dual_sum, args))

    expect_identical(dual$ksum, scalar1$ksum)
    expect_identical(dual$ksum.power2, scalar2$ksum)
  }
})

test_that("dual-power entry rejects unsupported kernel-sum configurations", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  default_sum <- getFromNamespace("npksum.default", "npRmpi")
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  set.seed(20260719)
  x <- data.frame(x = rnorm(20))

  expect_error(
    local_eval(default_sum(
      bws = 0.5, txdat = x, tydat = rnorm(20),
      .np.internal.power12 = TRUE)),
    "invalid use"
  )
  expect_error(
    local_eval(default_sum(
      bws = 0.5, txdat = x, operator = "integral",
      .np.internal.power12 = TRUE)),
    "invalid use"
  )
})

test_that("weighted dual-power sums preserve expanded-sample semantics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")

  set.seed(20260720)
  x <- data.frame(
    c = rnorm(24),
    u = factor(rep(c("a", "b", "c", "a"), 6),
               levels = c("a", "b", "c")),
    o = ordered(rep(c(1, 2, 3, 2), 6), levels = 1:3)
  )
  counts <- c(1, 3, 2, 4, rep(1, 20))
  bw <- c(0.6, 0.2, 0.3)
  weighted_sum <- getFromNamespace(".npksum_power12_weighted", "npRmpi")
  dual_sum <- getFromNamespace(".npksum_power12", "npRmpi")

  weighted <- local_eval(weighted_sum(
    bws = bw, txdat = x, counts = counts,
    bandwidth.divide = TRUE
  ))
  scalar1 <- local_eval(npksum(
    bws = bw, txdat = x, tydat = counts,
    bandwidth.divide = TRUE
  ))
  scalar2 <- local_eval(npksum(
    bws = bw, txdat = x, tydat = counts, kernel.pow = 2,
    bandwidth.divide = TRUE
  ))
  expect_identical(weighted$ksum, scalar1$ksum)
  expect_identical(weighted$ksum.power2, scalar2$ksum)

  expanded <- x[rep(seq_len(nrow(x)), counts), , drop = FALSE]
  expanded_sum <- local_eval(dual_sum(
    bws = bw, txdat = expanded, leave.one.out = TRUE,
    bandwidth.divide = TRUE
  ))
  diagonal <- local_eval(dual_sum(
    bws = bw, txdat = x[1L, , drop = FALSE],
    bandwidth.divide = TRUE
  ))
  compressed1 <- sum(counts *
    (as.numeric(weighted$ksum) - as.numeric(diagonal$ksum)))
  compressed2 <- sum(counts *
    (as.numeric(weighted$ksum.power2) -
     as.numeric(diagonal$ksum.power2)))
  expect_equal(compressed1, sum(expanded_sum$ksum), tolerance = 1e-12)
  expect_equal(compressed2, sum(expanded_sum$ksum.power2), tolerance = 1e-12)

  expect_error(
    weighted_sum(bws = bw, txdat = x, counts = c(counts[-1L], 0)),
    "positive finite"
  )
  expect_error(
    local_eval(weighted_sum(
      bws = 6, txdat = x["c"], counts = counts,
      bwtype = "generalized_nn"
    )),
    "invalid use"
  )
})

test_that("count compression accepts numeric fixed bandwidths only by default", {
  eligible <- getFromNamespace(
    ".npdeneq_count_compression_eligible",
    "npRmpi"
  )
  expect_true(eligible(0.2))
  expect_false(eligible(structure(list(), class = "not-a-bandwidth")))
})
