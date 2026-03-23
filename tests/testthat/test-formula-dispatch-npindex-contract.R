test_that("named bws formula dispatch matches positional npindex route", {
  set.seed(20260323)
  dat <- data.frame(
    y = rnorm(24),
    x1 = runif(24),
    x2 = runif(24)
  )
  ex <- data.frame(
    x1 = seq(0.1, 0.9, length.out = 8),
    x2 = seq(0.9, 0.1, length.out = 8)
  )

  fit.pos <- npindex(
    y ~ x1 + x2,
    data = dat,
    newdata = ex,
    method = "ichimura"
  )
  fit.named <- npindex(
    bws = y ~ x1 + x2,
    data = dat,
    newdata = ex,
    method = "ichimura"
  )

  expect_equal(as.numeric(fit.named$mean), as.numeric(fit.pos$mean), tolerance = 0)
  expect_equal(as.numeric(fit.named$merr), as.numeric(fit.pos$merr), tolerance = 0)
})
