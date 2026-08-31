h8_direct <- function(...) {
  package <- environmentName(environment(npreg))
  getFromNamespace(".np_regression_direct", package)(...)
}

h8_bw <- function(x, y, degree) {
  npregbw(
    xdat = x, ydat = y,
    bws = c(0.42, 0.39, 0.25, 0.3),
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = "fixed", regtype = "lp", degree = degree,
    degree.select = "manual"
  )
}

h8_fixture <- function(n = 19L) {
  x <- data.frame(
    x1 = seq(-1, 1, length.out = n),
    x2 = (((seq_len(n) * 5L) %% n) + 0.5) / n,
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high"), length.out = n),
                levels = c("low", "mid", "high"))
  )
  y <- sin(1.4 * x$x1) + 0.25 * x$x2 + 0.35 * (x$u == "b") -
    0.2 * (x$o == "high") + seq_len(n) / 200
  list(x = x, y = y, ex = x[c(1L, 5L, 10L, 15L, n), , drop = FALSE])
}

test_that("direct regression helper consumes the complete native HC0 payload", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)
  fixture <- h8_fixture()

  for (degree in list(c(2L, 2L), c(2L, 0L), c(0L, 0L))) {
    bw <- h8_bw(fixture$x, fixture$y, degree)
    public <- suppressWarnings(npreg(
      bws = bw, txdat = fixture$x, tydat = fixture$y,
      exdat = fixture$ex, se = TRUE, gradients = TRUE,
      warn.glp.gradient = FALSE
    ))
    direct <- suppressWarnings(h8_direct(
      bws = bw, txdat = fixture$x, tydat = fixture$y,
      exdat = fixture$ex, se = TRUE, gradients = TRUE,
      local.mode = TRUE
    ))

    expect_identical(direct$mean, public$mean)
    expect_identical(direct$merr, public$merr)
    expect_identical(direct$grad, public$grad)
    expect_identical(direct$gerr, public$gerr)
    expect_true(all(is.finite(direct$gerr[, 3:4, drop = FALSE])))
  }

  direct.body <- paste(deparse(body(h8_direct), width.cutoff = 500L),
                       collapse = " ")
  implementation <- paste(deparse(body(getFromNamespace(
    ".np_regression_direct", environmentName(environment(npreg))
  )), width.cutoff = 500L), collapse = " ")
  expect_false(grepl("npreghat", implementation, fixed = TRUE))
  expect_match(implementation,
               "do.compiled.gradients <- isTRUE(gradients)", fixed = TRUE)
  expect_true(nzchar(direct.body))
})

test_that("regression accessors expose the same cumulative HC0 result", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  fixture <- h8_fixture()
  bw <- h8_bw(fixture$x, fixture$y, c(2L, 2L))
  fit <- npreg(
    bws = bw, txdat = fixture$x, tydat = fixture$y,
    se = TRUE, gradients = TRUE, residuals = TRUE
  )

  expect_identical(fitted(fit), fit$mean)
  expect_identical(se(fit), fit$merr)
  expect_identical(gradients(fit), fit$grad)
  expect_identical(gradients(fit, se = TRUE), fit$gerr)
  expect_equal(residuals(fit), fixture$y - fit$mean, tolerance = 0)

  predicted <- predict(fit, exdat = fixture$ex, se.fit = TRUE)
  explicit <- npreg(
    bws = bw, txdat = fixture$x, tydat = fixture$y,
    exdat = fixture$ex, se = TRUE
  )
  expect_identical(predicted$fit, explicit$mean)
  expect_identical(predicted$se.fit, explicit$merr)
  expect_output(print(fit), "Regression Data", fixed = TRUE)
  expect_output(summary(fit), "Regression Data", fixed = TRUE)
})

test_that("categorical asymptotic plot panels consume native HC0 errors", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  fixture <- h8_fixture()
  bw <- h8_bw(fixture$x, fixture$y, c(2L, 2L))
  fit <- npreg(
    bws = bw, txdat = fixture$x, tydat = fixture$y,
    se = TRUE, gradients = TRUE
  )
  out <- suppressWarnings(plot(
    fit, xdat = fixture$x, ydat = fixture$y,
    output = "data", perspective = FALSE, gradients = TRUE,
    errors = "asymptotic", data_overlay = FALSE, neval = 7L
  ))

  expect_type(out, "list")
  expect_true(all(vapply(out, inherits, logical(1), "npregression")))
  expect_true(all(is.finite(out[[3L]]$grad)))
  expect_true(all(is.finite(out[[3L]]$gerr)))
  expect_true(all(is.finite(out[[4L]]$grad)))
  expect_true(all(is.finite(out[[4L]]$gerr)))
})

test_that("formula NA restoration preserves mean and every gradient SE column", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  fixture <- h8_fixture(23L)
  dat <- cbind(y = fixture$y, fixture$x)
  dat$y[3L] <- NA_real_
  dat$x1[7L] <- NA_real_
  newdata <- fixture$x[c(1L, 4L, 8L, 12L, 17L, 23L), , drop = FALSE]
  newdata$x2[2L] <- NA_real_
  newdata$u[5L] <- NA

  bw <- npregbw(
    y ~ x1 + x2 + u + o, data = dat,
    bws = c(0.42, 0.39, 0.25, 0.3),
    bandwidth.compute = FALSE, bwscaling = FALSE,
    regtype = "lp", degree = c(2L, 2L), degree.select = "manual"
  )
  fit <- npreg(
    bws = bw, data = dat, newdata = newdata,
    se = TRUE, gradients = TRUE
  )

  omitted <- c(2L, 5L)
  kept <- setdiff(seq_len(nrow(newdata)), omitted)
  expect_identical(nrow(fit$grad), nrow(newdata))
  expect_identical(nrow(fit$gerr), nrow(newdata))
  expect_true(all(is.na(fit$mean[omitted])))
  expect_true(all(is.na(fit$merr[omitted])))
  expect_true(all(is.na(fit$grad[omitted, , drop = FALSE])))
  expect_true(all(is.na(fit$gerr[omitted, , drop = FALSE])))
  expect_true(all(is.finite(fit$mean[kept])))
  expect_true(all(is.finite(fit$merr[kept])))
  expect_true(all(is.finite(fit$grad[kept, , drop = FALSE])))
  expect_true(all(is.finite(fit$gerr[kept, , drop = FALSE])))
})
