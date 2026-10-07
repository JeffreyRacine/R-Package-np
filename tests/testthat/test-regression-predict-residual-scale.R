# These tests own prediction metadata; numerical estimator tests stay in their owners.
test_that("prediction residual scale follows known training MSE", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(664)
  d <- data.frame(x = rnorm(120), z = runif(120))
  d$y <- sin(d$x + .3*d$z) + rnorm(120, sd = .3)
  e <- d[91:115, ]; e$y <- e$y + 3
  w <- dnorm(outer(d$x, d$x, `-`)/.7) * dnorm(outer(d$z, d$z, `-`)/.8)
  training.mse <- mean((d$y - drop(w %*% d$y/rowSums(w)))^2)
  for (formula in c(FALSE, TRUE)) {
    b <- if (formula) npregbw(y ~ x + z, data = d, bws = c(.7, .8),
                             bandwidth.compute = FALSE) else
      npregbw(xdat = d[c("x", "z")], ydat = d$y, bws = c(.7, .8),
              bandwidth.compute = FALSE)
    f <- npreg(b, se = TRUE)
    p <- predict(f, newdata = e[c("x", "z")], se.fit = TRUE)
    scored <- predict(f, exdat = e[c("x", "z")], eydat = e$y, se.fit = TRUE)
    expect_equal(f$MSE, training.mse, tolerance = 1e-10)
    expect_equal(p$residual.scale, training.mse, tolerance = 1e-10)
    expect_equal(scored$residual.scale, training.mse, tolerance = 1e-10)
    expect_identical(p$fit, scored$fit)
    expect_identical(p$se.fit, scored$se.fit)
    expect_identical(predict(f, newdata = e[c("x", "z")]), p$fit)
    expect_equal(predict(f, eydat = d$y + 3, se.fit = TRUE)$residual.scale,
                 training.mse, tolerance = 1e-10)

    external <- npreg(b, exdat = e[c("x", "z")], eydat = e$y)
    expect_true(is.na(predict(external, newdata = e[c("x", "z")],
                             se.fit = TRUE)$residual.scale))
    alternate <- npreg(b, eydat = d$y + 3)
    expect_true(alternate$trainiseval)
    expect_true(is.na(predict(alternate, newdata = e[c("x", "z")],
                             se.fit = TRUE)$residual.scale))
    expect_equal(predict(alternate, se.fit = TRUE)$residual.scale,
                 training.mse, tolerance = 1e-10)

    tr <- d[1:87, ]; tr$y <- tr$y + cos(tr$x)
    replace.args <- if (formula) list(data = tr) else
      list(txdat = tr[c("x", "z")], tydat = tr$y)
    replacement <- do.call(npreg, c(list(bws = b), replace.args))
    at.training <- do.call(predict, c(list(object = f, se.fit = TRUE), replace.args))
    at.external <- do.call(predict, c(list(object = f, se.fit = TRUE,
                                          newdata = e[c("x", "z")]), replace.args))
    expect_equal(at.training$residual.scale, replacement$MSE)
    expect_true(is.na(at.external$residual.scale))
    expect_equal(predict(replacement, newdata = e[c("x", "z")], se.fit = TRUE)$residual.scale,
                 replacement$MSE)
    expect_true(is.na(predict(f, newdata = e[c("x", "z")],
                             na.action = na.exclude, se.fit = TRUE)$residual.scale))
  }
})
