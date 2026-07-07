test_that("npreg rejects non-logical gradients and residuals under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
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
  npreg_rbandwidth <- getFromNamespace("npreg.rbandwidth", "npRmpi")
  fn.body <- paste(deparse(body(npreg_rbandwidth), width.cutoff = 500L), collapse = " ")
  pos.grad <- regexpr("gradients <- npValidateScalarLogical\\(gradients, \"gradients\"\\)", fn.body)[1]
  pos.resid <- regexpr("residuals <- npValidateScalarLogical\\(residuals, \"residuals\"\\)", fn.body)[1]
  pos.auto <- regexpr("if \\(\\.npRmpi_autodispatch_active\\(\\)\\)", fn.body)[1]

  expect_true(pos.grad > 0L && pos.resid > 0L && pos.auto > 0L)
  expect_true(pos.grad < pos.auto)
  expect_true(pos.resid < pos.auto)
})

test_that("lp regtype remains lp internally for degree 0/1", {
  npRegtypeToC <- getFromNamespace("npRegtypeToC", "npRmpi")
  REGTYPE_LP <- getFromNamespace("REGTYPE_LP", "npRmpi")

  expect_identical(npRegtypeToC(regtype = "lp", degree = 0L, ncon = 1L)$code, REGTYPE_LP)
  expect_identical(npRegtypeToC(regtype = "lp", degree = 1L, ncon = 1L)$code, REGTYPE_LP)
})

test_that("lp degree-0 gradients allow first derivatives and reject higher orders under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260305)
  dat <- data.frame(y = rnorm(30), x = runif(30))
  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 0
  )

  expect_identical(bw$regtype, "lp")
  expect_no_error(
    fit <- npRmpi::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 1L)
  )
  expect_error(
    npRmpi::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 2L),
    "no available derivative components"
  )
  expect_no_error(
    value.fit <- npRmpi::npreg(bws = bw, data = dat, gradients = FALSE)
  )
  expect_identical(fit$bws$regtype, "lp")
  expect_identical(value.fit$bws$regtype, "lp")
})

test_that("lc gradients reject higher-order derivative requests under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260617)
  dat <- data.frame(y = rnorm(40), x = runif(40))
  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  expect_no_error(
    fit <- npRmpi::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 1L)
  )
  expect_error(
    npRmpi::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 2L),
    "supports only first derivatives for regtype='lc'"
  )
  expect_error(
    npRmpi::gradients(fit, gradient.order = 2L),
    "supports only first derivatives for regtype='lc'"
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot(fit, gradients = TRUE, gradient_order = 2L, output = "data"),
    "supports only first derivatives for regtype='lc'"
  )
})

test_that("npreg.rbandwidth no longer contains bernstein OOS direct fallback", {
  npreg_rbandwidth <- getFromNamespace("npreg.rbandwidth", "npRmpi")
  fn.body <- paste(deparse(body(npreg_rbandwidth), width.cutoff = 500L), collapse = " ")
  expect_no_match(fn.body, "bernstein\\.oos")
  expect_no_match(fn.body, "\\.np_kernel_weights_direct")
  expect_no_match(fn.body, "\\.npreghat_solve_eval")
})

test_that("lp bernstein OOS warns but keeps canonical npreg path under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260305)
  dat <- data.frame(y = rnorm(40), x = runif(40))
  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    bws = 0.35,
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2,
    bernstein.basis = TRUE
  )
  ex <- data.frame(x = c(-0.2, 1.2))
  expect_warning(
    fit <- npRmpi::npreg(bws = bw, txdat = dat["x"], tydat = dat$y, exdat = ex),
    "outside training support"
  )
  expect_identical(fit$bws$regtype, "lp")
})
