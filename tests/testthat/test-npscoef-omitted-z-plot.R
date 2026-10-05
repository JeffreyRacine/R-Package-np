test_that("smooth-coefficient plots and hats use x metadata when z is omitted", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(2403)
  x <- data.frame(x = runif(60))
  y <- sin(3 * x$x) + rnorm(60, sd = 0.1)
  fit <- npscoef(txdat = x, tydat = y, bws = 0.4, se = FALSE)
  expect_null(fit$bws[["zdati", exact = TRUE]])
  applied <- npscoefhat(bws = fit$bws, txdat = x, y = y, output = "apply")
  expect_equal(as.vector(applied), fitted(fit), tolerance = 1e-8)
  panel <- plot(fit, output = "data", neval = 9L)[[1L]]
  expect_equal(nrow(panel$eval), 9L)
  expected <- npscoef(bws = fit$bws, txdat = x, tydat = y,
                      exdat = panel$eval, se = FALSE)
  expect_equal(as.vector(panel$mean), fitted(expected), tolerance = 1e-8)

  formula.fit <- npscoef(y ~ x, data = data.frame(y = y, x = x$x), bws = 0.4)
  formula.panel <- plot(formula.fit, output = "data", neval = 9L)[[1L]]
  expect_equal(as.vector(formula.panel$mean),
               as.vector(predict(formula.fit, newdata = formula.panel$eval)),
               tolerance = 1e-8)
})
