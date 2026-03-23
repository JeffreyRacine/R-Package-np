test_that("named bws formula dispatch matches positional semiparametric routes", {
  set.seed(20260323)
  dat <- data.frame(
    y = rnorm(24),
    x1 = runif(24),
    z1 = runif(24)
  )
  ex <- data.frame(
    x1 = seq(0.1, 0.8, length.out = 7),
    z1 = seq(0.8, 0.1, length.out = 7)
  )

  sc.pos <- npscoef(y ~ x1 | z1, data = dat, newdata = ex, errors = FALSE, iterate = FALSE)
  sc.named <- npscoef(bws = y ~ x1 | z1, data = dat, newdata = ex, errors = FALSE, iterate = FALSE)
  expect_equal(as.numeric(sc.named$mean), as.numeric(sc.pos$mean), tolerance = 0)
  expect_equal(as.numeric(sc.named$merr), as.numeric(sc.pos$merr), tolerance = 0)

  pl.pos <- npplreg(y ~ x1 | z1, data = dat, newdata = ex)
  pl.named <- npplreg(bws = y ~ x1 | z1, data = dat, newdata = ex)
  expect_equal(as.numeric(pl.named$mean), as.numeric(pl.pos$mean), tolerance = 0)
  expect_equal(as.numeric(pl.named$merr), as.numeric(pl.pos$merr), tolerance = 0)
})
