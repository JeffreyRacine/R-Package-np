test_that("npregbw preserves deprecated remin compatibility", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260515)
  n <- 40L
  x <- rnorm(n)
  y <- x + rnorm(n)

  set.seed(20260516)
  expect_warning(
    bw.legacy <- npregbw(y ~ x, nmulti = 1, remin = FALSE),
    "deprecated"
  )
  set.seed(20260516)
  bw.modern <- npregbw(y ~ x, nmulti = 1, powell.remin = FALSE)

  expect_equal(as.numeric(bw.legacy$bw), as.numeric(bw.modern$bw),
               tolerance = 1e-12)
  expect_warning(
    fit.legacy <- npreg(y ~ x, nmulti = 1, remin = FALSE),
    "deprecated"
  )
  expect_s3_class(fit.legacy, "npregression")
})
