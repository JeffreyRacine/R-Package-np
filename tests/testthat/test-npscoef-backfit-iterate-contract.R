npscoef_backfit_iterate_fixture <- function(n = 24L) {
  set.seed(20260702)
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = sort(runif(n)))
  ydat <- (1 + zdat$z) * xdat$x + rnorm(n, sd = 0.05)

  list(xdat = xdat, zdat = zdat, ydat = ydat)
}

test_that("npscoefbw retains fitted bandwidth state for iterated backfit CV", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  fixture <- npscoef_backfit_iterate_fixture()

  bw <- suppressWarnings(npscoefbw(
    xdat = fixture$xdat,
    zdat = fixture$zdat,
    ydat = fixture$ydat,
    bws = 0.25,
    bandwidth.compute = TRUE,
    cv.iterate = TRUE,
    backfit.iterate = TRUE,
    cv.num.iterations = 1L,
    backfit.maxiter = 2L,
    optim.maxit = 2L,
    nmulti = 1L,
    random.seed = 11L
  ))

  expect_s3_class(bw, "scbandwidth")
  expect_equal(dim(bw$bw.fitted), c(length(bw$bw), ncol(fixture$xdat) + 1L))
  expect_true(all(is.finite(bw$bw.fitted)))
  expect_true(is.finite(as.numeric(bw$fval.fitted)))
})

test_that("npscoef iterated fits preserve residual identity when errors are requested", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  fixture <- npscoef_backfit_iterate_fixture()
  bw <- suppressWarnings(npscoefbw(
    xdat = fixture$xdat,
    zdat = fixture$zdat,
    ydat = fixture$ydat,
    bws = 0.25,
    bandwidth.compute = TRUE,
    cv.iterate = TRUE,
    backfit.iterate = TRUE,
    cv.num.iterations = 1L,
    backfit.maxiter = 2L,
    optim.maxit = 2L,
    nmulti = 1L,
    random.seed = 11L
  ))

  fit <- NULL
  expect_warning(
    fit <- npscoef(
      bws = bw,
      txdat = fixture$xdat,
      tzdat = fixture$zdat,
      tydat = fixture$ydat,
      iterate = TRUE,
      residuals = TRUE,
      errors = TRUE,
      maxiter = 50L,
      tol = 1e-8
    ),
    "standard errors are not available for iterated npscoef fits"
  )

  expect_s3_class(fit, "smoothcoefficient")
  expect_lt(max(abs(as.numeric(fit$mean) + as.numeric(fit$resid) - fixture$ydat)), 1e-10)
  expect_true(all(is.na(fit$merr)))
  expect_true(all(is.na(fit$gerr)))
})

test_that("npscoef nearest-neighbor candidate bandwidths reject invalid k values", {
  old_opts <- options(np.extendednn = FALSE)
  on.exit(options(old_opts), add = TRUE)

  disabled <- npRmpi:::.npscoef_nn_candidate_bandwidth(
    param = c(0, 1, 2, 100),
    bwtype = "adaptive_nn",
    nobs = 10L
  )
  expect_equal(disabled, c(NA_real_, NA_real_, 2, NA_real_))

  options(np.extendednn = TRUE)
  enabled <- npRmpi:::.npscoef_nn_candidate_bandwidth(
    param = c(0, 1, 2, 100),
    bwtype = "adaptive_nn",
    nobs = 10L
  )
  expect_equal(enabled, c(NA_real_, NA_real_, 2, 100))
})

test_that("npscoef fast large-h gate rejects unsafe kernel metadata and non-finite evals", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(
    np.messages = FALSE,
    np.largeh.rel.tol = 0.05,
    np.disc.upper.rel.tol = 0.05
  )
  on.exit(options(old_opts), add = TRUE)

  fixture <- npscoef_backfit_iterate_fixture()
  make_bw <- function(...) {
    npscoefbw(
      xdat = fixture$xdat,
      zdat = fixture$zdat,
      ydat = fixture$ydat,
      bws = 1e6,
      bandwidth.compute = FALSE,
      ...
    )
  }

  safe <- make_bw(ckertype = "gaussian", ckerorder = 2L)
  expect_true(npRmpi:::.npscoefbw_fast_eligible(safe, fixture$zdat))

  zbad <- fixture$zdat
  zbad$z[1L] <- NA_real_
  expect_false(npRmpi:::.npscoefbw_fast_eligible(safe, zbad))
  expect_false(npRmpi:::.npscoefbw_fast_eligible(
    make_bw(ckertype = "epanechnikov", ckerorder = 4L),
    fixture$zdat
  ))
  expect_false(npRmpi:::.npscoefbw_fast_eligible(
    make_bw(ckertype = "truncated gaussian", ckerorder = 2L),
    fixture$zdat
  ))
  expect_false(npRmpi:::.npscoefbw_fast_eligible(
    make_bw(ckertype = "gaussian", ckerorder = 2L, ckerbound = "range"),
    fixture$zdat
  ))
})
