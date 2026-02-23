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
  fn.body <- paste(deparse(body(npreg.rbandwidth), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "gradients <- npValidateScalarLogical\\(gradients, \"gradients\"\\)")
  expect_match(fn.body, "residuals <- npValidateScalarLogical\\(residuals, \"residuals\"\\)")
})
