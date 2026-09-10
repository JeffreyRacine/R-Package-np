independent_lp_delete_one_weights <- function(xdat,
                                              bandwidth,
                                              degree,
                                              row_index) {
  x <- as.numeric(xdat[[1L]])
  powers <- 0L:as.integer(degree)
  basis <- outer(x, powers, `^`)
  weights <- dnorm((x[[row_index]] - x) / bandwidth)
  weights[[row_index]] <- 0

  gram <- crossprod(basis, basis * weights)
  projection <- solve(gram, basis[row_index, ])
  as.numeric(weights * (basis %*% projection))
}

bounded_gaussian_kernel_nprmpi <- function(x0, X, h, lower, upper) {
  denom <- h * (pnorm((upper - x0) / h) - pnorm((lower - x0) / h))
  dnorm((x0 - X) / h) / denom
}

test_that("npcdensbw cv.ml LP degree-0 bounded objective matches delete-one reconstruction", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(1)
  n <- 80L
  xdat <- data.frame(x = runif(n))
  ydat <- rbeta(n, 1, 1)

  bw <- npRmpi::npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ml",
    regtype = "lp",
    degree = 0L,
    bws = c(0.15, 0.12),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )

  manual_rows <- vapply(seq_len(n), function(i) {
    weights <- independent_lp_delete_one_weights(
      xdat,
      bandwidth = bw$xbw[[1L]],
      degree = bw$degree.engine[[1L]],
      row_index = i
    )
    ky <- bounded_gaussian_kernel_nprmpi(
      ydat[i],
      ydat,
      bw$ybw[1L],
      bw$cykerlb[bw$iycon][1L],
      bw$cykerub[bw$iycon][1L]
    )
    sum(weights * ky)
  }, numeric(1))

  expect_true(all(manual_rows > .Machine$double.xmin))

  manual_objective <- sum(log(manual_rows))
  np_objective <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)$objective

  expect_equal(np_objective, manual_objective, tolerance = 1e-5)
})

test_that("npcdensbw cv.ml LP rejects negative delete-one fits as raw invalid", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(1)
  n <- 80L
  xdat <- data.frame(x = runif(n))
  ydat <- rbeta(n, 1, 1)

  bw <- npRmpi::npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ml",
    regtype = "lp",
    degree = 3L,
    bws = c(0.15, 0.12),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )

  manual_rows <- vapply(seq_len(n), function(i) {
    weights <- independent_lp_delete_one_weights(
      xdat,
      bandwidth = bw$xbw[[1L]],
      degree = bw$degree.engine[[1L]],
      row_index = i
    )
    ky <- bounded_gaussian_kernel_nprmpi(
      ydat[i],
      ydat,
      bw$ybw[1L],
      bw$cykerlb[bw$iycon][1L],
      bw$cykerub[bw$iycon][1L]
    )
    sum(weights * ky)
  }, numeric(1))

  expect_gt(sum(manual_rows < 0), 0L)

  raw <- npRmpi:::.npcdensbw_eval_only(
    xdat, ydat, bw, invalid.penalty = "dbmax")
  guided <- npRmpi:::.npcdensbw_eval_only(
    xdat, ydat, bw, invalid.penalty = "baseline")
  expect_identical(as.numeric(raw$objective), -.Machine$double.xmax)
  expect_true(is.finite(guided$objective))
  expect_lt(abs(guided$objective), .Machine$double.xmax)
})
