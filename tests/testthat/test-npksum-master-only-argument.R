test_that("npksum materializes master-only extra arguments before pooled work", {
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)),
          "This transport contract requires actual pooled autodispatch")
  if (!spawn_mpi_slaves()) skip("An MPI pool is required")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  expect_true(.mpi_pool_active())
  set.seed(26135)
  x <- data.frame(x = rnorm(80))
  local.call <- function() {
    master.only <- 4L
    count <- 0L
    value <- npksum(txdat = x, bws = .4,
                   ckerorder = { count <- count + 1L; master.only })
    reference <- npksum(txdat = x, bws = .4)
    fourth <- npksum(txdat = x, bws = .4, ckerorder = 4L)
    expect_identical(count, 1L)
    expect_equal(value$ksum, fourth$ksum, tolerance = 0)
    expect_false(isTRUE(all.equal(fourth$ksum, reference$ksum)))
    expect_error(npksum(txdat = x, bws = .4,
                       foo = { count <- count + 1L; stop("wrong forcing") }),
                 "unused npksum argument")
    expect_identical(count, 1L)
    expect_error(npksum(txdat = x, bws = .4,
                       ckerorder = stop("master argument sentinel")),
                 "master argument sentinel")
    expect_equal(npksum(txdat = x, bws = .4)$ksum, reference$ksum,
                 tolerance = 0)
  }
  local.call()
})
