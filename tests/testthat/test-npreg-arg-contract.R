test_that("npreg rejects non-logical gradients and residuals", {
  set.seed(20260223)
  dat <- data.frame(y = rnorm(30), x = runif(30))
  bw <- np::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  expect_error(np::npreg(bws = bw, data = dat, gradients = "foo"),
               "'gradients' must be TRUE or FALSE")
  expect_error(np::npreg(bws = bw, data = dat, residuals = "bar"),
               "'residuals' must be TRUE or FALSE")
})

test_that("npreg.rbandwidth validates gradient flags as logical scalars", {
  npreg_rbandwidth <- getFromNamespace("npreg.rbandwidth", "np")
  fn.body <- paste(deparse(body(npreg_rbandwidth), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "gradients <- npValidateScalarLogical\\(gradients, \"gradients\"\\)")
  expect_match(fn.body, "residuals <- npValidateScalarLogical\\(residuals, \"residuals\"\\)")
})

test_that("lp regtype remains lp internally for degree 0/1", {
  npRegtypeToC <- getFromNamespace("npRegtypeToC", "np")
  REGTYPE_LP <- getFromNamespace("REGTYPE_LP", "np")

  expect_identical(npRegtypeToC(regtype = "lp", degree = 0L, ncon = 1L)$code, REGTYPE_LP)
  expect_identical(npRegtypeToC(regtype = "lp", degree = 1L, ncon = 1L)$code, REGTYPE_LP)
  expect_identical(npRegtypeToC(regtype = "ll", degree = NULL, ncon = 2L),
                   list(code = REGTYPE_LP, degree = c(1L, 1L)))
})

test_that("lp degree-0 gradients allow first derivatives and reject higher orders", {
  set.seed(20260305)
  dat <- data.frame(y = rnorm(30), x = runif(30))
  bw <- np::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 0
  )

  expect_identical(bw$regtype, "lp")
  expect_no_error(
    fit <- np::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 1L)
  )
  expect_error(
    np::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 2L),
    "no available derivative components"
  )
  expect_no_error(
    value.fit <- np::npreg(bws = bw, data = dat, gradients = FALSE)
  )
  expect_identical(fit$bws$regtype, "lp")
  expect_identical(value.fit$bws$regtype, "lp")
})

test_that("lc gradients reject higher-order derivative requests", {
  set.seed(20260617)
  dat <- data.frame(y = rnorm(40), x = runif(40))
  bw <- np::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  expect_no_error(
    fit <- np::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 1L)
  )
  expect_error(
    np::npreg(bws = bw, data = dat, gradients = TRUE, gradient.order = 2L),
    "supports only first derivatives for regtype='lc'"
  )
  expect_error(
    np::gradients(fit, gradient.order = 2L),
    "supports only first derivatives for regtype='lc'"
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot(fit, gradients = TRUE, gradient_order = 2L, output = "data"),
    "supports only first derivatives for regtype='lc'"
  )
})

test_that("lp bernstein OOS warns but keeps canonical npreg path", {
  set.seed(20260305)
  dat <- data.frame(y = rnorm(40), x = runif(40))
  bw <- np::npregbw(
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
    fit <- np::npreg(bws = bw, txdat = dat["x"], tydat = dat$y, exdat = ex),
    "outside training support"
  )
  expect_identical(fit$bws$regtype, "lp")
})

test_that("npreg.rbandwidth no longer contains bernstein OOS direct fallback", {
  npreg_rbandwidth <- getFromNamespace("npreg.rbandwidth", "np")
  fn.body <- paste(deparse(body(npreg_rbandwidth), width.cutoff = 500L), collapse = " ")
  expect_no_match(fn.body, "bernstein\\.oos")
  expect_no_match(fn.body, "\\.np_kernel_weights_direct")
  expect_no_match(fn.body, "\\.npreghat_solve_eval")
})
