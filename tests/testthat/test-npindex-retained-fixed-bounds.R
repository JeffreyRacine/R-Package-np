test_that("single-index fitting carries retained scalar kernel bounds", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.04, .96, length.out = 29L))
  y <- sin(5*x$x) + .05*cos(seq_len(29L))
  e <- data.frame(x = c(.1, .4, .9))
  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn"))
    for (regtype in c("lc", "ll", "lp")) {
      h <- if (bwtype == "fixed") .25 else 15
      controls <- list(bwtype = bwtype, regtype = regtype,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1, bandwidth.compute = FALSE)
      if (regtype == "lp") controls$degree <- 2L
      expect_warning(b <- do.call(npindexbw,
        c(list(xdat = x, ydat = y, bws = c(1, h)), controls)), "one dimension")
      rb <- do.call(npregbw, c(list(xdat = x, ydat = y, bws = h), controls))
      actual <- npindex(b, txdat = x, tydat = y, exdat = e, se = FALSE, gradients = TRUE)
      plain <- npindex(b, txdat = x, tydat = y, exdat = e, se = FALSE)
      oracle <- npreg(rb, txdat = x, tydat = y, exdat = e, gradients = TRUE, se = FALSE)
      expect_equal(fitted(actual), fitted(oracle), tolerance = 1e-11)
      expect_equal(as.numeric(gradients(actual)), as.numeric(gradients(oracle)), tolerance = 1e-10)
      expect_equal(fitted(plain), fitted(oracle), tolerance = 1e-11)
      expect_equal(as.numeric(npindexhat(b, txdat = x, exdat = e, y = y, output = "apply")),
        fitted(oracle), tolerance = 1e-11)
      expect_equal(as.numeric(predict(actual, exdat = e, se = FALSE)), fitted(oracle), tolerance = 1e-11)
      expect_error(npindex(b, txdat = x, tydat = y, exdat = data.frame(x = 1.1), se = FALSE),
        "bounds")
    }
})

test_that("single-index inference and resampling retain the fixed domain", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.04, .96, length.out = 29L))
  y <- sin(5*x$x) + .05*cos(seq_len(29L))
  e <- data.frame(x = c(.1, .4, .9))
  for (kernel in c("gaussian", "epanechnikov")) for (regtype in c("lc", "ll", "lp")) {
    controls <- list(regtype = regtype, ckertype = kernel, ckerbound = "fixed",
      ckerlb = 0, ckerub = 1, bandwidth.compute = FALSE)
    if (regtype == "lp") controls$degree <- 2L
    expect_warning(b <- do.call(npindexbw,
      c(list(xdat = x, ydat = y, bws = c(1, .5)), controls)), "one dimension")
    rb <- do.call(npregbw, c(list(xdat = x, ydat = y, bws = .5), controls))
    actual <- npindex(b, txdat = x, tydat = y, exdat = e, se = TRUE, gradients = TRUE)
    oracle <- npreg(rb, txdat = x, tydat = y, exdat = e, se = TRUE, gradients = TRUE)
    expect_equal(se(actual), se(oracle), tolerance = 1e-11)
    for (gradient in c(FALSE, TRUE)) {
      set.seed(211)
      boot <- npindex(b, txdat = x, tydat = y, exdat = e, se = TRUE,
        gradients = gradient, se.type = "bootstrap", B = 3L)
      expect_true(all(is.finite(se(boot))))
      expect_equal(fitted(boot), fitted(oracle), tolerance = 1e-11)
    }
  }
})
