library(np)

test_that("npcdens formula fit keeps lc ll and lp estimator paths live", {
  set.seed(20260307)
  n <- 60L
  x <- rnorm(n)
  y <- x + rnorm(n)

  fit.lc <- npcdens(y ~ x, regtype = "lc", nmulti = 1)
  fit.ll <- npcdens(y ~ x, regtype = "ll", nmulti = 1)
  fit.lp0 <- npcdens(y ~ x, regtype = "lp", degree = 0L, nmulti = 1)
  fit.lp2 <- npcdens(y ~ x, regtype = "lp", degree = 2L, nmulti = 1)

  expect_true(all(is.finite(fitted(fit.lc))))
  expect_true(all(is.finite(fitted(fit.ll))))
  expect_true(all(is.finite(fitted(fit.lp0))))
  expect_true(all(is.finite(fitted(fit.lp2))))

  expect_equal(fitted(fit.lc), fitted(fit.lp0), tolerance = 1e-8)

  for (fit in list(fit.lc, fit.ll, fit.lp0, fit.lp2)) {
    expect_true(length(capture.output(summary(fit))) > 0L)

    pdf(file = tempfile(fileext = ".pdf"))
    on.exit(dev.off(), add = TRUE)
    expect_error(plot(fit, view = "fixed"), NA)
  }
})
