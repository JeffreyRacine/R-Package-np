test_that("implicit hat responses use the retained regression coding", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  for (labels in list(c("10", "30", "50"), c("alpha", "beta", "gamma"))) {
    y <- ordered(rep(labels, 8), levels = labels)
    numeric.y <- if (labels[1L] == "10") rep(c(10, 30, 50), 8) else as.double(y)
    d <- data.frame(x = seq(.05, .95, length.out = 24), y = y)
    for (type in c("fixed", "generalized_nn", "adaptive_nn")) for (rt in c("lc", "ll")) {
      settings <- list(bws = if (type == "fixed") .25 else 12,
                       bwtype = type, regtype = rt, bandwidth.compute = FALSE)
      b <- do.call(npregbw, c(list(xdat = d["x"], ydat = y), settings))
      bf <- do.call(npregbw, c(list(formula = y ~ x, data = d), settings))
      model <- npreg(b, se = FALSE)
      H <- npreghat(b, txdat = d["x"])
      expected <- as.vector(H %*% numeric.y)
      expect_equal(expected, as.numeric(fitted(model)), tolerance = 1e-11)
      for (input in list(b, bf, model)) {
        actual <- do.call(npreghat, list(bws = input, output = "apply"))
        expect_equal(as.numeric(actual), expected, tolerance = 1e-11)
        constrained <- do.call(npreghat, list(bws = input, output = "constraint"))
        expect_equal(as.numeric(constrained), as.numeric(t(H) * numeric.y), tolerance = 1e-11)
      }
      expect_equal(as.numeric(npreghat(b, y = as.double(y), output = "apply")),
                   as.numeric(H %*% as.double(y)), tolerance = 1e-11)
      payload <- cbind(numeric.y, as.double(y))
      expect_equal(as.numeric(npreghat(b, y = payload, output = "apply")),
                   as.numeric(H %*% payload), tolerance = 1e-11)
    }
  }
})

test_that("implicit factor response mapping preserves retained omission rows", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.1, .9, length.out = 24),
                  y = ordered(rep(c(10, 30, 50), 8)))
  d$x[3L] <- NA_real_; d$y[9L] <- NA
  b <- npregbw(y ~ x, data = d, bws = .3, bandwidth.compute = FALSE)
  keep <- complete.cases(d)
  H <- npreghat(b)
  expected <- as.numeric(H %*% as.numeric(as.character(d$y[keep])))
  expect_equal(as.numeric(npreghat(b, output = "apply")), expected, tolerance = 1e-11)
  expect_equal(ncol(H), sum(keep))
})
