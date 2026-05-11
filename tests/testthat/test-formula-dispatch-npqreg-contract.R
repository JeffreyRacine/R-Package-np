test_that("named bws formula dispatch matches explicit npqreg bandwidth route", {
  set.seed(20260323)
  d <- data.frame(
    x = seq(0.1, 0.9, length.out = 9),
    y = seq(0.2, 1.0, length.out = 9)
  )
  nd <- data.frame(x = c(0.2, 0.5, 0.8))

  bw <- npcdistbw(y ~ x, data = d, newdata = nd)
  fit.pos <- npqreg(bws = bw, data = d, newdata = nd, tau = 0.4)
  fit.named <- npqreg(bws = y ~ x, data = d, newdata = nd, tau = 0.4)

  expect_s3_class(fit.named, "qregression")
  expect_equal(length(fit.named$quantile), nrow(nd), tolerance = 0)
  expect_lt(max(abs(as.numeric(fit.named$quantile) - as.numeric(fit.pos$quantile))), 1e-2)
  expect_lt(max(abs(as.numeric(fit.named$quanterr) - as.numeric(fit.pos$quanterr))), 1e-2)
})

test_that("named bws formula dispatch forwards NOMAD shortcut controls", {
  skip_if_not_installed("crs")

  set.seed(20260511)
  d <- data.frame(
    x = seq(0.05, 0.95, length.out = 30),
    y = sin(seq(0.05, 0.95, length.out = 30)) + rnorm(30, sd = 0.05)
  )

  capture.output(
    fit.named <- npqreg(
      bws = y ~ x,
      data = d,
      tau = 0.5,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L
    )
  )

  expect_s3_class(fit.named, "qregression")
  expect_true(isTRUE(fit.named$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit.named$bws$regtype.engine), "lp")
})
