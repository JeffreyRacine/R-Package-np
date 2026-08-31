test_that("npsigtest upper-tail owner counts equality only", {
  upper_tail <- getFromNamespace(".np_npsig_upper_tail_p", "npRmpi")

  expect_identical(upper_tail(rep(0, 9), 0), 1)
  expect_identical(upper_tail(c(0, 0, 1), 0), 1)

  bootstrap <- c(-2, -1, 0, 1, 2)
  observed <- 0.5
  expect_identical(
    upper_tail(bootstrap, observed),
    mean(bootstrap > observed)
  )
})

test_that("npsigtest reports a smoothed-out categorical point mass as nonrejection under MPI", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_options <- options(
    np.messages = FALSE,
    np.largelambda = TRUE,
    npRmpi.autodispatch = TRUE
  )
  on.exit(options(old_options), add = TRUE)

  set.seed(280828)
  n <- 40L
  x <- seq(-1, 1, length.out = n)
  z <- factor(rep(c("a", "b"), length.out = n))
  y <- x + rnorm(n, sd = 0.1)
  bw <- npregbw(
    xdat = data.frame(x, z),
    ydat = y,
    bws = c(0.35, 0.5),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  result <- npsigtest(
    bw,
    index = 2L,
    pivot = FALSE,
    boot.num = 9L,
    random.seed = 42L
  )

  expect_identical(unname(result$In), 0)
  expect_true(all(result$In.bootstrap == 0))
  expect_identical(unname(result$P), 1)
  expect_identical(unname(result$reject), "")
  expect_true(is.na(unname(result$rejectNum)))
})
