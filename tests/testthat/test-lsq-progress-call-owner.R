test_that("vector LSQ search re-enters its private method with owned controls", {
  set.seed(291)
  n <- 61L
  x <- data.frame(x = rnorm(n), z = runif(n))
  y <- sin(x$x) + x$z + rnorm(n, sd = .2)
  start <- npregbw(xdat = x, ydat = y, bws = c(.6, .6),
                   bandwidth.compute = FALSE)
  count <- 0L
  fit <- nplsqreg(start, txdat = x, tydat = y, tau = .5,
    scale = { count <- count + 1L; rep(.3, n) }, bandwidth.compute = FALSE)
  expect_identical(count, 1L)
  expect_length(fitted(fit), n)
  count <- 0L
  bw <- nplsqregbw(start, xdat = { count <- count + 1L; x },
                   ydat = y, tau = c(.2, .5, .8),
                   scale = rep(.3, n), nmulti = 1L, itmax = 2L)
  expect_identical(count, 1L)
  expect_s3_class(bw, "lsqregressionbandwidth")
  expect_equal(bw$tau, c(.2, .5, .8), tolerance = 0)
  expect_identical(bw$call[[1L]], as.name("nplsqregbw"))
  fit <- nplsqreg(bw, se = TRUE)
  expect_identical(dim(fitted(fit)), c(n, 3L))
  expect_identical(dim(predict(fit, exdat = x[1:9, ], se.fit = TRUE)$fit), c(9L, 3L))
})
