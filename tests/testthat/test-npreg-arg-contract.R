test_that("npreg rejects non-logical gradients and residuals under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260223)
  dat <- data.frame(y = rnorm(30), x = runif(30))
  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  expect_error(npRmpi::npreg(bws = bw, data = dat, gradients = "foo"),
               "'gradients' must be TRUE or FALSE")
  expect_error(npRmpi::npreg(bws = bw, data = dat, residuals = "bar"),
               "'residuals' must be TRUE or FALSE")
})

test_that("npreg.rbandwidth validates gradient flags before autodispatch branch", {
  fn.body <- paste(deparse(body(npreg.rbandwidth), width.cutoff = 500L), collapse = " ")
  pos.grad <- regexpr("gradients <- as\\.logical\\(gradients\\)", fn.body)[1]
  pos.resid <- regexpr("residuals <- as\\.logical\\(residuals\\)", fn.body)[1]
  pos.auto <- regexpr("if \\(\\.npRmpi_autodispatch_active\\(\\)\\)", fn.body)[1]

  expect_true(pos.grad > 0L && pos.resid > 0L && pos.auto > 0L)
  expect_true(pos.grad < pos.auto)
  expect_true(pos.resid < pos.auto)
})
