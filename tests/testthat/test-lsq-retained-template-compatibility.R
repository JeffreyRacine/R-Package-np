test_that("LSQ conversion preserves legacy formula readers and kernel templates", {
  set.seed(727)
  d <- data.frame(x = runif(65), z = runif(65), y = rnorm(65))
  b <- npregbw(y ~ x + z, data = d, bws = c(.6, .7), regtype = "lp",
    degree = c(1L, 1L), ckertype = "epanechnikov", ckerbound = "range",
    bandwidth.compute = FALSE)
  legacy <- b
  legacy[[".np.formula.training"]] <- NULL
  args <- list(bws = legacy, bandwidth.compute = FALSE,
               scale = rep(1, 65), delta = .5)
  old <- do.call(nplsqreg, args)
  args$bws <- b
  fit <- do.call(nplsqreg, args)
  expect_equal(fitted(old), fitted(fit), tolerance = 0)
  for (field in c("bw", "regtype", "basis", "degree", "bernstein.basis",
                  "type", "ckertype", "ckerorder", "ckerbound", "ckerlb", "ckerub"))
    expect_equal(fit$reg.bws[[field]], b[[field]], tolerance = 0)
  expect_equal(predict(fit, newdata = d[1:7, ]),
               predict(old, newdata = d[1:7, ]), tolerance = 0)
})
