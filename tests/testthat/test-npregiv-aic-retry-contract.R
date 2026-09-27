test_that("IV AIC optimization retries retain the polynomial degree", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(42)
  n <- 64L
  w <- rnorm(n)
  z <- w + rnorm(n)
  y <- z + rnorm(n)
  # This bounded search reaches the retry that formerly dropped degree.
  fit <- withCallingHandlers(
    npregiv(y = y, z = z, w = w, p = 2L, bwmethod = "cv.aic",
            nmulti = 1L, iterate.max = 2L, optim.maxit = 2L,
            optim.maxattempts = 1L),
    warning = function(w) {
      if (grepl("optim failed to converge|iterate.max reached|Stopping rule increases",
                conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(fit, "npregiv")
  expect_length(fit$phi, n)
  expect_true(all(is.finite(fit$phi)))
})
