test_that("smooth-coefficient prediction and lazy residuals replay fitting controls", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(190907)
  x <- data.frame(x = runif(30))
  z <- data.frame(z = sort(runif(30)))
  y <- (1 + z$z)*x$x + rnorm(30, sd = .05)
  bw <- suppressWarnings(npscoefbw(xdat = x, zdat = z, ydat = y, bws = .3,
    bandwidth.compute = TRUE, cv.iterate = TRUE, backfit.iterate = TRUE,
    cv.num.iterations = 1L, backfit.maxiter = 2L, optim.maxit = 2L,
    nmulti = 1L, random.seed = 11L))
  expect_false(is.null(bw$bw.fitted))
  for (it in c(FALSE, TRUE)) {
    fit <- suppressWarnings(npscoef(bw, iterate = it, maxiter = 2L, tol = .001))
    expect_identical(fit$fit.controls,
      list(iterate = it, maxiter = 2L, tol = .001, leave.one.out = FALSE))
    expect_equal(suppressWarnings(predict(fit)), fitted(fit), tolerance = 1e-12)
    expect_equal(as.vector(suppressWarnings(residuals(fit))),
                 as.vector(y - fitted(fit)), tolerance = 1e-12)
    direct <- suppressWarnings(npscoef(bw, iterate = !it, maxiter = 1L, tol = .1))
    expect_equal(suppressWarnings(predict(fit, iterate = !it, maxiter = 1L, tol = .1)),
                 fitted(direct), tolerance = 1e-12)
  }
  fit <- npscoef(bw, iterate = FALSE, leave.one.out = TRUE)
  expect_equal(predict(fit), fitted(fit), tolerance = 1e-12)
  expect_equal(as.vector(residuals(fit)), as.vector(y - fitted(fit)), tolerance = 1e-12)
  fit$fit.controls <- NULL
  expect_equal(suppressWarnings(predict(fit)),
               suppressWarnings(fitted(npscoef(bw))), tolerance = 1e-12)
})
