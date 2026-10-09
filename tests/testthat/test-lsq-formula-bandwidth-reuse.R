test_that("LSQ formula regression bandwidths retain data and specification", {
  set.seed(952)
  n <- 80L
  d <- data.frame(x = rnorm(n), z = runif(n))
  d$y <- sin(d$x) + d$z + rnorm(n, sd = .2)
  b <- npregbw(y ~ x + z, data = d, bws = c(.55, .55),
               bandwidth.compute = FALSE)
  ref <- nplsqreg(b, txdat = d[c("x", "z")], tydat = d$y,
                  bandwidth.compute = FALSE, scale = rep(1, n), delta = .5)
  for (with.data in c(FALSE, TRUE)) {
    args <- list(bws = b, bandwidth.compute = FALSE,
                 scale = rep(1, n), delta = .5)
    if (with.data) args$data <- d
    fit <- do.call(nplsqreg, args)
    expect_equal(fit$reg.bws$bw, b$bw, tolerance = 0)
    expect_equal(fitted(fit), fitted(ref), tolerance = 0)
    expect_s3_class(fit$bws$formula, "formula")
    expect_equal(predict(fit, newdata = d[1:8, ]),
                 predict(ref, exdat = d[1:8, c("x", "z")]), tolerance = 0)
    expect_equal(fitted(nplsqreg(fit$bws)), fitted(fit), tolerance = 0)
  }
  saved <- d
  d <- NULL
  fit <- nplsqreg(b, bandwidth.compute = FALSE, scale = rep(1, n), delta = .5)
  expect_equal(fitted(fit), fitted(ref), tolerance = 0)
  expect_equal(fit$bws$xdat, saved[c("x", "z")], tolerance = 0)
})

test_that("LSQ formula bandwidth reuse owns transformed and selected rows", {
  set.seed(753)
  n <- 67L
  d <- data.frame(x = rnorm(n), z = runif(n), y = rnorm(n))
  d$y[8L] <- NA_real_
  b <- npregbw(y ~ I(x^2) + z, data = d, subset = 3:64,
               na.action = na.exclude, bws = c(.65, .45),
               regtype = "ll", bandwidth.compute = FALSE)
  frame <- b[[".np.formula.training"]]$frame
  x <- frame[c("I(x^2)", "z")]
  y <- frame$y
  fit <- nplsqreg(b, bandwidth.compute = FALSE,
                  scale = rep(1, nrow(frame)), delta = .5)
  ref <- nplsqreg(b, txdat = x, tydat = y, bandwidth.compute = FALSE,
                  scale = rep(1, nrow(frame)), delta = .5)
  expect_equal(fit$reg.bws$bw, b$bw, tolerance = 0)
  expect_equal(fit$reg.bws$regtype, b$regtype)
  expect_equal(fitted(fit)[-as.integer(attr(frame, "na.action"))],
               fitted(ref), tolerance = 0)
  expect_identical(fit$rows.omit, as.integer(attr(frame, "na.action")))
  replacement <- d[1:55, ]
  replacement$y <- seq_len(55) / 55
  fresh <- nplsqreg(b, data = replacement, subset = 2:50,
                    bandwidth.compute = FALSE, scale = rep(1, 55), delta = .5)
  chosen <- model.frame(y ~ I(x^2) + z, data = replacement, subset = 2:50,
                        na.action = na.exclude)
  chosen[["I(x^2)"]] <- as.numeric(chosen[["I(x^2)"]])
  expected <- nplsqreg(b, txdat = chosen[c("I(x^2)", "z")], tydat = chosen$y,
                       bandwidth.compute = FALSE, scale = rep(1, nrow(chosen)), delta = .5)
  expect_equal(fitted(fresh), fitted(expected), tolerance = 0)
  expect_equal(fresh$bws$scale, rep(1, nrow(chosen)), tolerance = 0)
})

test_that("LSQ explicit subset precedence does not change later saved-sample reuse", {
  d <- data.frame(x = seq(.01, .99, length.out = 64L))
  d$y <- sin(4 * d$x) + d$x^2
  b <- npregbw(y ~ x, data = d, subset = x < .5, bws = .2,
               bandwidth.compute = FALSE)
  original <- b
  retained <- nplsqreg(b, bandwidth.compute = FALSE,
                       scale = rep(1, 32L), delta = .5)
  selections <- list(NULL, TRUE, 29:64, d$x > .5)
  rows <- list(seq_len(64L), seq_len(64L), 29:64, which(d$x > .5))
  for (i in seq_along(selections)) {
    fresh <- do.call(nplsqreg, list(bws = b, data = d,
      subset = selections[[i]], bandwidth.compute = FALSE,
      scale = rep(1, 64L), delta = .5))
    chosen <- rows[[i]]
    expected <- nplsqreg(b, txdat = d[chosen, "x", drop = FALSE],
                         tydat = d$y[chosen], bandwidth.compute = FALSE,
                         scale = rep(1, length(chosen)), delta = .5)
    expect_equal(fitted(fresh), fitted(expected), tolerance = 0)
    expect_equal(fresh$bws$xdat, d[chosen, "x", drop = FALSE], tolerance = 0)
    expect_equal(fresh$bws$ydat, expected$bws$ydat, tolerance = 0)
    expect_identical(b, original)
  }
  replay <- nplsqreg(b, data = d, bandwidth.compute = FALSE,
                     scale = rep(1, 64L), delta = .5)
  expect_equal(fitted(replay), fitted(retained), tolerance = 0)
  expect_equal(fitted(nplsqreg(b, bandwidth.compute = FALSE,
                              scale = rep(1, 32L), delta = .5)),
               fitted(retained), tolerance = 0)
})
