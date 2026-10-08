test_that("npksum materializes master-only extra arguments before pooled work", {
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)),
          "This transport contract requires actual pooled autodispatch")
  if (!spawn_mpi_slaves()) skip("An MPI pool is required")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  expect_true(.mpi_pool_active())
  set.seed(26135)
  x <- data.frame(x = rnorm(80))
  local.call <- function() {
    master.only <- 17L
    count <- 0L
    value <- npksum(txdat = x, bws = .4,
                   foo = { count <- count + 1L; master.only })
    reference <- npksum(txdat = x, bws = .4)
    expect_identical(count, 1L)
    expect_equal(value$ksum, reference$ksum, tolerance = 0)
    expect_error(npksum(txdat = x, bws = .4,
                       foo = stop("master argument sentinel")),
                 "master argument sentinel")
    expect_equal(npksum(txdat = x, bws = .4)$ksum, reference$ksum,
                 tolerance = 0)
  }
  local.call()
})
