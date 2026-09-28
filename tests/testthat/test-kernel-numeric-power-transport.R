test_that("numeric kernel sums separate paired computation from construction", {
  numeric_sum <- getFromNamespace("npksum.numeric", "npRmpi")
  probe <- new.env(parent = environment(numeric_sum))
  environment(numeric_sum) <- probe
  probe$kbandwidth <- function(...) {
    probe$constructor <- list(...)
    structure(list(), class = "kbandwidth")
  }
  probe$npksum.default <- function(...) list(...)
  # This seam test examines only argument partitioning, not MPI lifecycle.
  probe$.npRmpi_require_active_slave_pool <- function(...) invisible(NULL)
  probe$.npRmpi_master_local_entry_needed <- function(...) FALSE
  probe$.npRmpi_npksum_should_localize <- function(...) FALSE
  probe$.npRmpi_autodispatch_active <- function(...) FALSE

  x <- data.frame(x = c(0, 0.5, 1))
  for (value in list(TRUE, FALSE, NULL)) {
    args <- list(bws = 0.3, txdat = x, bwtype = "fixed",
                 .np.internal.power12 = value,
                 .np.internal.power12.weighted = value)
    got <- do.call(numeric_sum, args)
    expect_true(all(c(".np.internal.power12",
                      ".np.internal.power12.weighted") %in% names(got)))
    expect_identical(got[[".np.internal.power12"]], value)
    expect_identical(got[[".np.internal.power12.weighted"]], value)
    expect_false(any(c(".np.internal.power12",
                       ".np.internal.power12.weighted") %in%
                     names(probe$constructor)))
    expect_identical(probe$constructor$bwtype, "fixed")
  }
  ordinary <- numeric_sum(bws = 0.3, txdat = x)
  expect_identical(names(ordinary), c("txdat", "bws"))
  derivative <- numeric_sum(bws = 0.3, txdat = x,
                            return.derivative.kernel.weights = TRUE)
  expect_true(derivative$return.derivative.kernel.weights)
  expect_false("return.derivative.kernel.weights" %in% names(probe$constructor))
})
