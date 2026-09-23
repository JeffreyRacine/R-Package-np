test_that("quantile predictions retain resolved extraction controls", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(-1, 1, length.out = 19L))
  d$y <- d$x + .2*cos(seq_len(nrow(d)))
  controls <- list(tol = 1e-9, small = 1e-10, itmax = 80L)
  for (route in c("formula", "native")) for (tau in list(.5, c(.25, .65))) {
    b <- if (route == "formula")
      npcdistbw(y ~ x, data = d, bws = c(.4, .6), bandwidth.compute = FALSE) else
      npcdistbw(xdat = d["x"], ydat = d["y"], bws = c(.4, .6), bandwidth.compute = FALSE)
    fit <- do.call(npqreg, c(list(bws = b, tau = tau), controls))
    expect_identical(fit[["fit.controls", exact = TRUE]], controls)
    expect_identical(predict(fit), fitted(fit))
    expect_identical(predict(fit, newdata = d["x"]), fitted(fit))
    expect_identical(predict(unserialize(serialize(fit, NULL)), newdata = d["x"]), fitted(fit))
    expect_identical(predict(fit, tol = .02),
      fitted(npqreg(b, tau = tau, tol = .02, small = controls$small, itmax = controls$itmax)))
    expect_identical(predict(fit, small = .02),
      fitted(npqreg(b, tau = tau, tol = controls$tol, small = .02, itmax = controls$itmax)))
    expect_error(predict(fit, tol = NULL), "'tol'")
    expect_error(predict(fit, itmax = 1L), "failed to converge")
    legacy <- fit; legacy$fit.controls <- NULL
    expect_identical(predict(legacy), fitted(npqreg(b, tau = tau)))
    expect_identical(predict(fit, se.fit = TRUE)$fit, fitted(fit))
  }
})

test_that("quantile plot fit construction records its actual default controls", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(-1, 1, length.out = 15L))
  d$y <- d$x + .2*sin(seq_len(nrow(d)))
  b <- npcdistbw(y ~ x, data = d, bws = c(.4, .6), bandwidth.compute = FALSE)
  fit <- npqreg(b, tau = c(.3, .6))
  out <- plot(fit, output = "data", perspective = FALSE, errors = "none", neval = 3L)
  expect_true(length(out) > 0L)
  for (panel in out) {
    expect_identical(panel[["fit.controls", exact = TRUE]],
      list(tol = 1.490116e-4, small = 1.490116e-5, itmax = 10000L))
    expect_identical(predict(panel, exdat = panel$xeval), fitted(panel))
  }
})
