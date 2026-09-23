test_that("quantile prediction keeps native coordinates out of formula transforms", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = exp(seq(-.8, .8, length.out = 21)))
  d$y <- 2 * log(d$x) + .2 * cos(seq_len(nrow(d)))
  for (shape in c("ordinary", "log")) {
    f <- if (shape == "log") y ~ log(x) else y ~ x
    nd <- data.frame(x = exp(c(-.5, -.1, .4)))
    native <- data.frame(if (shape == "log") log(nd$x) else nd$x)
    names(native) <- if (shape == "log") "log(x)" else "x"
    positional <- native; names(positional) <- "x"
    for (type in c("fixed", "generalized_nn")) {
      b <- npcdistbw(f, data = d, bws = if (type == "fixed") c(.25, .3) else c(8, 8),
                    bwtype = type, bandwidth.compute = FALSE)
      for (tau in list(.5, c(.3, .7))) {
        fit <- npqreg(b, tau = tau)
        oracle <- fitted(npqreg(b, exdat = native, tau = tau))
        expect_identical(predict(fit, newdata = nd), oracle)
        expect_identical(predict(fit, exdat = native), oracle)
        expect_identical(predict(fit, exdat = positional), oracle)
        expect_identical(predict(fit, exdat = native, newdata = data.frame(wrong = 1:3)), oracle)
        # Explicit NULL is an invalid native argument, not a missing argument.
        expect_error(predict(fit, exdat = NULL, newdata = nd), "exdat")
        expect_error(predict(fit, newdata = data.frame(wrong = 1:3)), "columns")
        missing <- native; missing[2, 1] <- NA_real_
        reference <- npqreg(b, exdat = missing, tau = tau, se = TRUE)
        pred <- predict(fit, exdat = missing, se.fit = TRUE)
        expect_identical(pred$fit, fitted(reference))
        expect_identical(pred$se.fit, se(reference))
        expect_true(all(is.na(if (is.null(dim(pred$fit))) pred$fit[2] else pred$fit[2, ])))
        expect_identical(predict(unserialize(serialize(fit, NULL)), exdat = native), oracle)
      }
    }
  }
})
