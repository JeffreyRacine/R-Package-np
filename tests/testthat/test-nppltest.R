test_that("partial-linearity test is exported and both result methods are registered", {
  ns <- asNamespace("npRmpi")
  expect_identical(getExportedValue("npRmpi", "nppltest"), ns$nppltest)
  for (generic in c("print", "summary"))
    expect_identical(getS3method(generic, "pltest", envir = ns),
                     get(paste0(generic, ".pltest"), envir = ns, inherits = FALSE))
  expect_error(nppltest(NULL, B = TRUE), "B must")
  expect_error(nppltest(NULL), "object must")
})

test_that("mixed training operator uses the fitted categorical weights", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(481)
  n <- 30L
  x <- data.frame(x = rnorm(n), w = factor(rep(1:2, length.out = n)))
  z <- data.frame(group = factor(rep(1:2, each = n/2)), v = rnorm(n))
  y <- rnorm(n)
  responses <- cbind(y, cos(y), seq_len(n)/n)
  for (reg in c("lc", "ll")) {
    bw <- npplregbw(xdat = x, ydat = y, zdat = z, bandwidth.compute = FALSE,
                    regtype = reg, bws = rbind(c(.496, .8), c(.2, .55), c(.496, .6)))
    H <- .pltest_operator(bw, x, y, z)$H
    oracle <- apply(responses, 2L, function(v)
      as.numeric(fitted(npplreg(bws = bw, txdat = x, tydat = v, tzdat = z))))
    expect_equal(unname(H %*% responses), unname(oracle), tolerance = 1e-10)
    expect_equal(.pltest_residual_norms(H), rowSums((diag(n) - H)^2))
  }
})

test_that("selected angular columns reproduce every serial entry in anchor order", {
  V <- cbind(c(0, 0, 1, -1, 2, 1e-14), c(1, 1, 0, 2, -1, 1))
  full <- .Call(C_np_pltest_angles, V)
  for (ids in list(c(6L, 2L, 4L), c(3L, 5L, 1L), 1L)) {
    columns <- .Call(C_np_pltest_angle_columns, V, ids)
    for (b in seq_along(ids)) {
      rows <- seq.int(ids[b], nrow(V))
      expect_identical(columns[rows, b], full[rows, ids[b]])
    }
  }
  expect_error(.Call(C_np_pltest_angle_columns, V, c(0L, 1L)), "column index")
})
