test_that("shared NN errors name the variable without rejecting all ties", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(group = factor(rep(c("a", "b"), 6)),
                  good = seq_len(12) / 12, tied = c(rep(0, 4), 1:8))
  y <- 1 + sin(seq_len(12)) / 4 + x$good
  bw <- npregbw(xdat = x, ydat = y, bws = c(.2, 8, 4),
                 bwtype = "generalized_nn", regtype = "lc",
                 bandwidth.compute = FALSE)
  fit <- npreg(bws = bw, txdat = x, tydat = y, se = FALSE)
  expect_true(all(is.finite(fitted(fit))))
  check <- function(expr, excluded = 0L, k = 4L) {
    e <- tryCatch(force(expr), np_nn_zero_radius = identity)
    expect_s3_class(e, "np_nn_zero_radius")
    expect_s3_class(e, "simpleError")
    expect_identical(e$variable, "tied")
    expect_identical(e$lookup.k, k)
    expect_identical(e$excluded, excluded)
    expect_match(conditionMessage(e), "repeated values", fixed = TRUE)
    expect_match(conditionMessage(e), "zero literal radius", fixed = TRUE)
    expect_match(conditionMessage(e), 'bwtype="fixed", the default', fixed = TRUE)
  }
  check(npreg(bws = bw, txdat = x, tydat = y, gradients = TRUE, se = FALSE))
  check(predict(fit, newdata = x))
  check(np:::.np_regression_direct(bw, x, y, exdat = x, gradients = TRUE))
  check(npreghat(bws = bw, txdat = x, exdat = x, y = y, output = "apply"))
  check(npksum(txdat = x, exdat = x, bws = c(.2, 8, 4),
               bwtype = "generalized_nn"))

  for (joint in c(FALSE, TRUE)) {
    set.seed(527)
    seed <- .Random.seed
    e <- tryCatch(npsigtest(fit, B = 9L, joint = joint), np_nn_zero_radius = identity)
    expect_identical(e$variable, "tied")
    expect_match(e$stage, "unrestricted gradient evaluation", fixed = TRUE)
    expect_identical(.Random.seed, seed)
  }

  for (type in c("generalized_nn", "adaptive_nn", "fixed")) {
    safe <- npregbw(xdat = x, ydat = y,
      bws = if (type == "fixed") c(.2, .8, 2) else c(.2, 8, 5),
      bwtype = type, regtype = "lc", bandwidth.compute = FALSE)
    out <- npreg(bws = safe, txdat = x, tydat = y, exdat = x,
                  gradients = TRUE, se = FALSE)
    expect_true(all(is.finite(out$mean)))
    expect_true(all(is.finite(out$grad)))
  }
  bad <- npregbw(xdat = x, ydat = y, bws = c(.2, 8, 3),
    bwtype = "adaptive_nn", regtype = "lc", bandwidth.compute = FALSE)
  check(npreg(bws = bad, txdat = x, tydat = y), excluded = 1L, k = 3L)
  expect_identical(npreg(bws = bw, txdat = x, tydat = y, se = FALSE)$mean, fit$mean)
})

test_that("density and conditional owners use the same named condition", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(tied = c(0, 0, 0, 1, 2, 3) / 4)
  z <- data.frame(smooth = seq_len(6) / 6)
  for (type in c("generalized_nn", "adaptive_nn")) {
    for (distribution in c(FALSE, TRUE)) {
      bw.fun <- if (distribution) npudistbw else npudensbw
      fit.fun <- if (distribution) npudist else npudens
      bw <- bw.fun(dat = x, bwtype = type, bws = 2, bandwidth.compute = FALSE)
      e <- tryCatch(fit.fun(bws = bw, tdat = x, edat = x), np_nn_zero_radius = identity)
      expect_s3_class(e, "np_nn_zero_radius")
      expect_identical(e$variable, "tied")

      cbw.fun <- if (distribution) npcdistbw else npcdensbw
      cfit.fun <- if (distribution) npcdist else npcdens
      for (response in c(FALSE, TRUE)) for (regtype in c("lc", "ll", "lp"))
        for (kernel in c("gaussian", "beta")) {
        tx <- if (response) z else x
        ty <- if (response) x else z
        args <- list(xdat = tx, ydat = ty, bwtype = type,
          bws = if (response) c(2, 4) else c(4, 2), regtype = regtype,
          bandwidth.compute = FALSE)
        if (regtype == "lp") args$degree <- 3L
        if (kernel == "beta") args <- c(args, list(
          cxkertype = "beta", cykertype = "beta", cxkerbound = "fixed",
          cykerbound = "fixed", cxkerlb = 0, cxkerub = 1, cykerlb = 0, cykerub = 1))
        cbw <- do.call(cbw.fun, args)
        e <- tryCatch(cfit.fun(bws = cbw, txdat = tx, tydat = ty,
                               exdat = tx, eydat = ty), np_nn_zero_radius = identity)
        expect_s3_class(e, "np_nn_zero_radius")
        expect_identical(e$variable, "tied")
      }
    }
  }
})

test_that("NN context remains lazy and unrelated conditions propagate unchanged", {
  expect_identical(np:::.np_with_nn_radius_context(17L, stop("names forced"),
                                                   stop("context forced")), 17L)
  original <- structure(list(message = "unrelated", call = NULL),
                         class = c("unrelated_error", "error", "condition"))
  actual <- tryCatch(np:::.np_with_nn_radius_context(stop(original), "x"),
                    error = identity)
  expect_identical(actual, original)
})
