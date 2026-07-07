test_that("npdeptest basic functionality works", {
  skip_on_cran()
  run_slow <- tolower(Sys.getenv("NP_RUN_SLOW_NPTESTS_MPI", ""))
  skip_if_not(
    run_slow %in% c("1", "true", "yes"),
    "slow npRmpi dependency-test bootstrap smoke; set NP_RUN_SLOW_NPTESTS_MPI=true to run"
  )
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 20
  x <- rnorm(n)
  y <- x + rnorm(n, sd=0.1)
  # Basic smoke only: keep n small and use the minimum bootstrap count.
  test <- npdeptest(x, y, method="summation", boot.num=9)
  
  expect_s3_class(test, "deptest")
  expect_output(summary(test))
})

test_that("npsdeptest basic functionality works", {
  skip_on_cran()
  run_slow <- tolower(Sys.getenv("NP_RUN_SLOW_NPTESTS_MPI", ""))
  skip_if_not(
    run_slow %in% c("1", "true", "yes"),
    "slow npRmpi serial-dependency bootstrap smoke; set NP_RUN_SLOW_NPTESTS_MPI=true to run"
  )
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 20
  y <- arima.sim(n=n, list(ar=0.5))
  test <- npsdeptest(y, lag.num=1, method="summation", boot.num=9)
  
  expect_s3_class(test, "sdeptest")
  expect_output(summary(test))
})
