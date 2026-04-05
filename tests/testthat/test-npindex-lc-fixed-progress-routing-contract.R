make_npindex_lc_fixed_progress_fixture <- function() {
  set.seed(20260405)
  n <- 24L
  dat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- dat$x1 + 0.6 * dat$x2
  dat$y <- sin(index) + 0.2 * index^2 + rnorm(n, sd = 0.05)

  list(
    dat = dat,
    tx = dat[c("x1", "x2")],
    y = dat$y,
    newdata = dat[c(2L, 7L), c("x1", "x2")],
    bw_fixed = npindexbw(
      xdat = dat[c("x1", "x2")],
      ydat = dat$y,
      bws = c(1, 0.6, 0.35),
      method = "ichimura",
      regtype = "lc",
      bwtype = "fixed",
      bandwidth.compute = FALSE
    ),
    bw_largeh = npindexbw(
      xdat = dat[c("x1", "x2")],
      ydat = dat$y,
      bws = c(1, 0.6, 5),
      method = "ichimura",
      regtype = "lc",
      bwtype = "fixed",
      bandwidth.compute = FALSE
    )
  )
}

npindex_lc_fixed_scalar_index <- function(txdat, beta) {
  data.frame(index = as.vector(as.matrix(txdat) %*% beta))
}

npindex_lc_fixed_oracle <- function(index.tx, y, bws, index.ex = index.tx) {
  spec <- getFromNamespace(".npindex_resolve_spec", "npRmpi")(bws, where = "npindex")
  regbw <- getFromNamespace(".np_index_regression_bw_state", "npRmpi")(
    index.df = index.tx,
    ydat = y,
    bws = bws
  )
  npreg(
    exdat = index.ex,
    bws = regbw,
    txdat = index.tx,
    tydat = y
  )
}

test_that("npRmpi fixed lc direct and predict routes preserve the scalar-index oracle", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_lc_fixed_progress_fixture()

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  fit <- npindex(
    bws = fixture$bw_fixed,
    txdat = fixture$tx,
    tydat = fixture$y
  )

  index.tx <- npindex_lc_fixed_scalar_index(fixture$tx, fit$beta)
  oracle <- npindex_lc_fixed_oracle(index.tx, fixture$y, fit$bws)

  pred <- predict(fit, newdata = fixture$newdata)
  index.ex <- npindex_lc_fixed_scalar_index(fixture$newdata, fit$beta)
  pred.oracle <- npindex_lc_fixed_oracle(index.tx, fixture$y, fit$bws, index.ex = index.ex)

  expect_equal(as.numeric(fit$mean), as.numeric(oracle$mean), tolerance = 1e-12)
  expect_equal(as.numeric(pred), as.numeric(pred.oracle$mean), tolerance = 1e-12)
})

test_that("npRmpi fixed lc large-h route preserves the scalar-index oracle", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_lc_fixed_progress_fixture()

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  fit <- npindex(
    bws = fixture$bw_largeh,
    txdat = fixture$tx,
    tydat = fixture$y
  )

  index.tx <- npindex_lc_fixed_scalar_index(fixture$tx, fit$beta)
  oracle <- npindex_lc_fixed_oracle(index.tx, fixture$y, fit$bws)

  expect_equal(as.numeric(fit$mean), as.numeric(oracle$mean), tolerance = 1e-12)
})

test_that("npRmpi fixed lc bw-to-fit route preserves the scalar-index oracle", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_lc_fixed_progress_fixture()

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260405)
  fit <- npindex(
    y ~ x1 + x2,
    data = fixture$dat,
    method = "ichimura",
    regtype = "lc",
    bwtype = "fixed",
    nmulti = 1L
  )

  index.tx <- npindex_lc_fixed_scalar_index(fixture$tx, fit$beta)
  oracle <- npindex_lc_fixed_oracle(index.tx, fixture$y, fit$bws)

  expect_equal(as.numeric(fit$mean), as.numeric(oracle$mean), tolerance = 1e-12)
})
