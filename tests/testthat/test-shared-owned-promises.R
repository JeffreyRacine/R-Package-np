test_that("pooled constructors use the master's first argument realization", {
  skip_on_cran()
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
  if (!spawn_mpi_slaves()) skip("An MPI pool is required")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(np.messages = FALSE)
  set.seed(5460)
  x <- data.frame(x = rnorm(90), z = runif(90, -1, 1))
  y <- sin(x$x) + .5*x$z + rnorm(90, sd = .2)
  set.seed(2024); first <- runif(1, -1, 1); set.seed(2024)
  count <- 0L
  b <- npindexbw(xdat = x, ydat = y,
                bws = { count <- count + 1L; c(1, runif(1, -1, 1), .3) },
                bandwidth.compute = FALSE)
  expect_identical(count, 1L)
  expect_identical(unname(b$beta[2L]), first)
  local({
    metadata <- c("A", "B")
    count <- 0L
    b <- npregbw(xdat = x, ydat = y, bws = c(.4, .5),
      bandwidth.compute = FALSE, xnames = { count <- count + 1L; metadata })
    ref <- npregbw(xdat = x, ydat = y, bws = c(.4, .5),
                  bandwidth.compute = FALSE, xnames = c("A", "B"))
    expect_identical(count, 1L)
    expect_identical(b$bw, ref$bw)
    expect_identical(b$xnames, ref$xnames)
    expect_error(npregbw(xdat = x, ydat = y, bws = c(.4, .5),
      bandwidth.compute = FALSE, xnames = stop("metadata failure")), "metadata failure")
    expect_equal(fitted(npreg(ref)), fitted(npreg(b)), tolerance = 0)
  })
})
