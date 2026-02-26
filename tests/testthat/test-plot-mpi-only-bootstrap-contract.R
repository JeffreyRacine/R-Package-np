test_that("bootstrap fanout is rejected inside mpi.bcast.cmd context", {
  fanout_enabled <- getFromNamespace(".npRmpi_bootstrap_fanout_enabled", "npRmpi")

  mpi.bcast.cmd <- function(expr) force(expr)

  expect_error(
    mpi.bcast.cmd(
      fanout_enabled(
        comm = 1L,
        n = 10L,
        B = 10L,
        chunk.size = 2L,
        what = "unit-test"
      )
    ),
    "cannot run inside mpi\\.bcast\\.cmd context"
  )
})

test_that("wild bootstrap helper requires MPI fanout and never falls back locally", {
  wild_boot_t <- getFromNamespace(".np_wild_boot_t", "npRmpi")

  withr::local_options(npRmpi.mpi.initialized = FALSE)

  H <- diag(2)
  fit.mean <- c(0.1, -0.2)
  residuals <- c(0.3, 0.4)

  expect_error(
    wild_boot_t(H = H, fit.mean = fit.mean, residuals = residuals, B = 3L, wild = "rademacher"),
    "requires an active MPI slave pool"
  )
})
