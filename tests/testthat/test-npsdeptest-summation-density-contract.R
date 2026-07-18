test_that("npsdeptest summation density helper matches npudens", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  density_fast <- getFromNamespace(".np_sdeptest_fixed_density", "npRmpi")

  cases <- list(
    list(
      data = c(-2.5, -1, -1, -0.2, 0.3, 1.8, 4.5),
      bw = 0.73
    ),
    list(
      data = cbind(
        x = c(-2.5, -1, -1, -0.2, 0.3, 1.8, 4.5),
        y = c(3.2, 0.1, 0.1, -1.7, 0.8, 2.4, 8.0)
      ),
      bw = c(0.73, 1.21)
    )
  )

  for (case in cases) {
    reference <- fitted(npudens(tdat = case$data, bws = case$bw))
    expect_identical(density_fast(case$data, case$bw), reference)
  }
})

test_that("npsdeptest summation route retains its public result shape", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(1927)
  data <- c(rexp(22) - 1, -2, -2, 7)
  out <- npsdeptest(
    data = data,
    lag.num = 2,
    method = "summation",
    bootstrap = TRUE,
    boot.num = 9,
    random.seed = 4817
  )

  expect_s3_class(out, "sdeptest")
  expect_identical(dim(out$Srho.bootstrap.mat), c(9L, 2L))
  expect_identical(dim(out$Srho.cumulant.bootstrap.mat), c(9L, 2L))
  expect_length(out$P, 2L)
  expect_length(out$P.cumulant, 2L)
})
