test_that("public npindex formula NOMAD shortcut completes on the small Ichimura smoke", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 50L
  dat <- data.frame(
    x1 = runif(n, min = -1, max = 1),
    x2 = runif(n, min = -1, max = 1)
  )
  dat$y <- dat$x1 - dat$x2 + rnorm(n)

  fit <- npindex(
    y ~ x1 + x2,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )

  expect_s3_class(fit, "singleindex")
  expect_s3_class(fit$bws, "sibandwidth")
  expect_true(is.finite(as.numeric(fit$bws$fval)))
  expect_true(is.finite(as.numeric(fit$bws$bw)))
})
