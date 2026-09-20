test_that("density-equality observed sums use the collective owner", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(193)
  x <- data.frame(x = rnorm(43))
  y <- data.frame(x = rnorm(27))
  original <- getFromNamespace(".npksum_power12", "npRmpi")
  observed.local <- logical()
  with_nprmpi_progress_bindings(list(.npksum_power12 = function(...) {
    args <- list(...)
    if (nrow(args$txdat) > 1L)
      observed.local <<- c(observed.local,
        isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    original(...)
  }), {
  result <- npdeneqtest(x, y, bw.x = .6, bw.y = .6, B = 9)
  })
  expect_identical(observed.local, rep(FALSE, 3L))
  sums <- function(a, b = a, h, loo = FALSE, power = 1L) {
    args <- list(txdat = a, bws = h, leave.one.out = loo,
                 bandwidth.divide = TRUE, kernel.pow = power)
    if (!loo) args$exdat <- b
    sum(do.call(npksum, args)$ksum)
  }
  n1 <- nrow(x); n2 <- nrow(y)
  In <- sums(x, h = .6, loo = TRUE) / (n1 * (n1 - 1)) +
    sums(y, h = .6, loo = TRUE) / (n2 * (n2 - 1)) -
    2 * sums(x, y, h = .6) / (n1 * n2)
  variance <- 2 * (sums(x, h = .6, loo = TRUE, power = 2L) / (n1^2 * (n1-1)^2) +
    sums(y, h = .6, loo = TRUE, power = 2L) / (n2^2 * (n2-1)^2) +
    2 * sums(x, y, h = .6, power = 2L) / (n1^2 * n2^2))
  expect_equal(result$In, In, tolerance = 2e-12)
  expect_equal(result$Tn, In / sqrt(variance), tolerance = 2e-12)
  expect_null(getFromNamespace(".np_progress_runtime", "npRmpi")$fit_state)
  expect_null(getFromNamespace(".np_progress_runtime", "npRmpi")$fit_forward)
})
