test_that("single-index refits report the current complete training sample", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(7121)
  x <- data.frame(x = rnorm(120), z = rnorm(120))
  responses <- list(ichimura = sin(x$x + .3*x$z),
                    kleinspady = rbinom(120, 1, plogis(x$x + .3*x$z)))
  for (method in names(responses)) {
    y <- responses[[method]]
    b <- npindexbw(xdat = x, ydat = y, method = method,
                   bws = c(1, .3, .8), bandwidth.compute = FALSE)
    train <- x[seq_len(95), ]; train$x[3] <- NA_real_
    for (compute in c(FALSE, TRUE)) {
      set.seed(22)
      refit <- npindexbw(xdat = train, ydat = y[seq_len(95)], bws = b,
        bandwidth.compute = compute, only.optimize.beta = TRUE,
        nmulti = 1L, optim.maxit = 20L)
      set.seed(22)
      direct <- npindexbw(xdat = train, ydat = y[seq_len(95)],
        method = method, bws = c(1, .3, .8), bandwidth.compute = compute,
        only.optimize.beta = TRUE, nmulti = 1L, optim.maxit = 20L)
      expect_equal(refit$nobs, 94L)
      expect_equal(refit$nobs.omit, 1L)
      expect_equal(as.integer(refit$rows.omit), 3L)
      expect_equal(coef(refit), coef(direct), tolerance = 1e-12)
      expect_equal(refit$bw, direct$bw, tolerance = 1e-12)
      expect_equal(refit$fval, direct$fval, tolerance = 1e-12)
      expect_match(paste(capture.output(print(refit)), collapse = " "), "94 observations")
      expect_match(paste(capture.output(summary(refit)), collapse = " "), "94 observations")
    }
    expect_equal(b$nobs, 120L)
  }
})

test_that("beta-only refits preserve held h without re-anchoring the search floor", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(6241)
  x <- data.frame(x = rnorm(90), z = rnorm(90))
  responses <- list(ichimura = sin(x$x + .3*x$z),
                    kleinspady = rep(c(0, 1), 45))
  for (method in names(responses)) {
    y <- responses[[method]]
    b <- npindexbw(xdat = x, ydat = y, method = method, regtype = "lc",
      bws = c(1, .3, .2), bandwidth.compute = FALSE,
      scale.factor.search.lower = 1)
    tx <- x[1:75, ]; ty <- y[1:75]
    floor <- .npindex_start_bandwidth_scale(as.matrix(tx) %*% b$beta, nrow(tx))
    expect_gt(floor, b$bw) # The witness must fail the old re-anchored guard.
    reference <- npindexbw(xdat = tx, ydat = ty, bws = b,
      only.optimize.beta = TRUE, nmulti = 1L, scale.factor.search.lower = 0)
    held <- npindexbw(xdat = tx, ydat = ty, bws = b,
      only.optimize.beta = TRUE, nmulti = 1L)
    expect_identical(held$bw, b$bw)
    expect_identical(held$beta, reference$beta)
    expect_identical(held$fval, reference$fval)
    expect_true(is.finite(held$fval))
    expect_identical(npGetScaleFactorSearchLower(held), 1)
    expect_error(npindexbw(xdat = tx, ydat = ty, bws = held,
      only.optimize.beta = FALSE, nmulti = 1L), "below the continuous")
    # An internal refinement owner may still impose its original physical floor.
    expect_error(npindexbw(xdat = tx, ydat = ty, bws = b,
      only.optimize.beta = TRUE, nmulti = 1L, .fixed.h.lower = .3),
      "below the continuous")
  }
})
