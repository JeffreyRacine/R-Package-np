test_that("named bws formula dispatch matches explicit npqreg bandwidth route", {
  set.seed(20260323)
  d <- data.frame(
    x = seq(0.1, 0.9, length.out = 9),
    y = seq(0.2, 1.0, length.out = 9)
  )
  nd <- data.frame(x = c(0.2, 0.5, 0.8))

  bw <- npcdistbw(y ~ x, data = d)
  fit.pos <- npqreg(bws = bw, data = d, newdata = nd, tau = 0.4)
  fit.named <- npqreg(bws = y ~ x, data = d, newdata = nd, tau = 0.4)

  expect_s3_class(fit.named, "qregression")
  expect_equal(length(fit.named$quantile), nrow(nd), tolerance = 0)
  expect_lt(max(abs(as.numeric(fit.named$quantile) - as.numeric(fit.pos$quantile))), 1e-2)
  expect_lt(max(abs(as.numeric(fit.named$quanterr) - as.numeric(fit.pos$quanterr))), 1e-2)
})

test_that("formula npqreg keeps native exdat out of bandwidth selection", {
  set.seed(20260514)
  d <- data.frame(
    x = seq(0.1, 0.9, length.out = 30)
  )
  d$y <- sin(2 * pi * d$x) + rnorm(nrow(d), sd = 0.08)
  nd <- data.frame(x = c(0.2, 0.5, 0.8))

  fit <- npqreg(
    y ~ x,
    data = d,
    exdat = nd,
    tau = 0.5,
    nmulti = 1
  )

  expect_s3_class(fit, "qregression")
  expect_equal(length(fitted(fit)), nrow(nd))
})

test_that("named bws formula dispatch forwards NOMAD shortcut controls", {
  skip_if_not_installed("crs")

  set.seed(20260511)
  d <- data.frame(
    x = seq(0.05, 0.95, length.out = 30),
    y = sin(seq(0.05, 0.95, length.out = 30)) + rnorm(30, sd = 0.05)
  )

  capture.output(
    fit.named <- npqreg(
      bws = y ~ x,
      data = d,
      tau = 0.5,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L
    )
  )

  expect_s3_class(fit.named, "qregression")
  expect_true(isTRUE(fit.named$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit.named$bws$regtype.engine), "lp")
})

test_that("npqreg NOMAD fit metadata propagates through prediction, gradients, and plot", {
  skip_if_not_installed("crs")

  count_namespace_calls <- function(fun, expr) {
    ns <- asNamespace("np")
    ctr <- new.env(parent = emptyenv())
    ctr$n <- 0L
    trace(
      what = fun,
      where = ns,
      tracer = bquote(.(ctr)$n <- .(ctr)$n + 1L),
      print = FALSE
    )
    on.exit(try(untrace(fun, where = ns), silent = TRUE), add = TRUE)
    force(expr)
    ctr$n
  }

  set.seed(20260511)
  n <- 36L
  d <- data.frame(x = seq(0.05, 0.95, length.out = n))
  d$y <- sin(2 * pi * d$x) + rnorm(n, sd = 0.08)
  nd <- data.frame(x = c(0.2, 0.4, 0.7))

  capture.output(
    fit <- npqreg(
      y ~ x,
      data = d,
      newdata = nd,
      tau = 0.5,
      gradients = TRUE,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L
    )
  )

  expect_s3_class(fit, "qregression")
  expect_true(isTRUE(fit$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit$bws$regtype.engine), "lp")
  expect_length(fitted(fit), nrow(nd))
  expect_equal(nrow(gradients(fit)), nrow(nd))
  expect_equal(nrow(gradients(fit, errors = TRUE)), nrow(nd))

  bw.calls <- count_namespace_calls("npcdistbw", {
    pred <- predict(fit, newdata = nd, tau = c(0.25, 0.5, 0.75))
    expect_equal(nrow(pred), nrow(nd))
    expect_equal(ncol(pred), 3L)
  })
  expect_identical(bw.calls, 0L)

  qreg.calls <- count_namespace_calls("npqreg", {
    out <- plot(
      fit,
      output = "data",
      perspective = FALSE,
      errors = "none",
      tau = c(0.25, 0.5, 0.75),
      neval = 6L
    )
    expect_type(out, "list")
  })
  expect_identical(qreg.calls, 0L)
})
