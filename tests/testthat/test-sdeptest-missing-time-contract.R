test_that("serial-dependence input never compresses missing time points", {
  ns <- asNamespace("npRmpi")
  fun <- get("npsdeptest", ns)
  env <- new.env(parent = ns)
  # Isolate input validation: any attempted dispatch, seed entry or bandwidth
  # selection would fail this focused test instead of masking the contract.
  env$.npRmpi_require_active_slave_pool <- function(...) invisible(TRUE)
  env$.npRmpi_autodispatch_active <- function(...) FALSE
  env$.np_seed_enter <- function(...) stop("seed entry reached")
  env$npudensbw <- function(...) stop("bandwidth search reached")
  environment(fun) <- env
  set.seed(54)
  seed <- .Random.seed
  for (pos in c(1L, 5L, 12L)) for (missing in c(NA_real_, NaN)) {
    data <- seq_len(12L) / 10
    data[pos] <- missing
    for (series in list(data, ts(data))) for (method in c("summation", "integration"))
      for (bootstrap in c(FALSE, TRUE)) {
        expect_error(fun(series, lag.num = 2, method = method, B = 9,
                         bootstrap = bootstrap),
                     "complete time series.*contiguous complete segment")
        expect_identical(.Random.seed, seed)
      }
  }
  expect_error(fun(rep(NA_real_, 12L), B = 9), "complete time series")
  expect_error(fun(seq_len(12L), B = 9), "seed entry reached")
})
