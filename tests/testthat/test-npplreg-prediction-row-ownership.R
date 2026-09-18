plreg_rows_fixture <- function() {
  set.seed(91617)
  d <- data.frame(y = rnorm(36), x = runif(36), z = runif(36))
  d$y <- .7 * d$x + sin(3 * d$z) + .3 * d$y
  rownames(d) <- paste0("training-", seq_len(nrow(d)))
  d
}

plreg_rows_bw <- function(d, action = na.exclude, regtype = "lc",
                          bwtype = "fixed", subset = seq_len(nrow(d))) {
  # Formula is positional; bws supplies manual, unsearched widths.
  args <- list(y ~ x | z, data = d, na.action = action, subset = subset,
    bandwidth.compute = FALSE, regtype = regtype, bwtype = bwtype)
  args$bws <- matrix(if (bwtype == "fixed") .6 else 12, 2L, 1L)
  if (regtype == "lp") args$degree <- 2L
  do.call(npplregbw, args)
}

plreg_rows_oracle <- function(bw, d, nd, rows = which(complete.cases(d))) {
  ev <- which(complete.cases(nd[c("x", "z")]))
  x <- d[rows, "x", drop = FALSE]
  z <- d[rows, "z", drop = FALSE]
  native <- npplreg(bw, txdat = x, tydat = d$y[rows], tzdat = z,
    exdat = nd[ev, "x", drop = FALSE], ezdat = nd[ev, "z", drop = FALSE],
    residuals = TRUE, se = TRUE)
  H <- npplreghat(bw, txdat = x, tzdat = z,
    exdat = nd[ev, "x", drop = FALSE], ezdat = nd[ev, "z", drop = FALSE],
    output = "matrix")
  out <- list(fit = rep(NA_real_, nrow(nd)), se.fit = rep(NA_real_, nrow(nd)))
  out$fit[ev] <- fitted(native)
  out$se.fit[ev] <- sqrt(pmax(drop((H^2) %*% (native$resid^2)), 0))
  list(pred = out, native = native)
}

test_that("plreg prediction SE preserves separate training and evaluation omissions", {
  old <- options(np.messages = FALSE, na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  d <- plreg_rows_fixture()
  d$y[3L] <- NA; d$x[8L] <- NA; d$z[14L] <- NA
  nd <- data.frame(x = c(.2, NA, .5, .7, .8), z = c(.3, .4, .6, NA, .5))
  rownames(nd) <- paste0("evaluation-", seq_len(nrow(nd)))
  for (regtype in c("lc", "ll", "lp")) for (bwtype in c("fixed", "generalized_nn")) {
    bw <- plreg_rows_bw(d, regtype = regtype, bwtype = bwtype)
    fit <- npplreg(bw, residuals = TRUE, se = TRUE)
    before <- fit[c("mean", "resid", "xcoef", "xcoeferr", "xcoefvcov", "omit")]
    oracle <- plreg_rows_oracle(bw, d, nd)
    pred <- predict(fit, newdata = nd, se.fit = TRUE)
    expect_equal(unname(pred$fit), oracle$pred$fit, tolerance = 0)
    expect_equal(unname(pred$se.fit), oracle$pred$se.fit, tolerance = 0)
    expect_identical(which(is.na(pred$se.fit)), c(2L, 4L))
    evaluation <- npplreg(bw, newdata = nd)
    expect_s3_class(evaluation$omit, "exclude")
    expect_identical(names(evaluation$omit), rownames(nd)[c(2L, 4L)])
    expect_identical(names(pred$se.fit), names(pred$fit))
    expect_identical(which(is.na(residuals(fit))), c(3L, 8L, 14L))
    expect_equal(unname(coef(fit)), unname(coef(oracle$native)), tolerance = 0)
    expect_equal(vcov(fit), vcov(oracle$native), tolerance = 0)
    expect_identical(fit[c("mean", "resid", "xcoef", "xcoeferr", "xcoefvcov", "omit")], before)
    expect_identical(predict(fit, newdata = nd), pred$fit)
    training <- predict(fit, se.fit = TRUE)
    # Training evaluation also excludes response-only missing rows.
    expected <- rep(NA_real_, nrow(d))
    keep <- which(complete.cases(d))
    expected[keep] <- plreg_rows_oracle(bw, d, d[keep, c("x", "z")])$pred$se.fit
    expect_equal(unname(training$se.fit), expected, tolerance = 0)
    expect_identical(which(is.na(training$fit)), which(is.na(training$se.fit)))
  }
})

test_that("plreg prediction SE uses subset row identities with omit and exclude", {
  old <- options(np.messages = FALSE, na.action = "na.omit")
  on.exit(options(old), add = TRUE)
  d <- plreg_rows_fixture()
  d$y[3L] <- NA; d$x[8L] <- NA; d$z[14L] <- NA
  selected <- c(2L, 3L, 5L, 8L, 10:32)
  nd <- data.frame(x = c(.2, NA, .7), z = c(.3, .4, .6))
  rows <- selected[complete.cases(d[selected, ])]
  for (action in list(na.omit, na.exclude)) {
    options(na.action = action)
    bw <- plreg_rows_bw(d, action = action, subset = selected)
    fit <- npplreg(bw, residuals = TRUE, se = TRUE)
    ans <- predict(fit, newdata = nd, se.fit = TRUE)
    oracle <- plreg_rows_oracle(bw, d, nd, rows)
    expected <- oracle$pred
    if (identical(action, na.omit)) expected <- lapply(expected, function(x) x[c(1L, 3L)])
    expect_equal(lapply(ans, unname), expected, tolerance = 0)
    expect_s3_class(fit$omit, if (identical(action, na.omit)) "omit" else "exclude")
    expect_identical(names(fit$omit), rownames(d)[selected[!complete.cases(d[selected, ])]])
    expect_equal(vcov(fit), vcov(oracle$native), tolerance = 0)
    expect_identical(predict(fit, newdata = nd), ans$fit)
  }
})

test_that("plreg prediction SE honors fresh training data and native precedence", {
  old <- options(np.messages = FALSE, na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  d <- plreg_rows_fixture()
  bw <- plreg_rows_bw(d)
  fit <- npplreg(bw, residuals = TRUE, se = TRUE)
  changed <- d
  changed$y <- d$y + .6 * d$x^2
  changed$x[3L] <- NA; changed$y[8L] <- NA; changed$z[14L] <- NA
  nd <- data.frame(x = c(.2, NA, .5, .7), z = c(.3, .4, .6, NA))
  oracle <- plreg_rows_oracle(bw, changed, nd)$pred
  for (object in list(fit, unserialize(serialize(fit, NULL)))) {
    ans <- predict(object, data = changed, newdata = nd, se.fit = TRUE)
    expect_equal(lapply(ans, unname), oracle, tolerance = 0)
    native <- predict(object, txdat = changed["x"], tydat = changed$y,
      tzdat = changed["z"], exdat = nd["x"], ezdat = nd["z"],
      newdata = data.frame(invalid = 1), se.fit = TRUE)
    expect_equal(lapply(native, unname), oracle, tolerance = 0)
    formula.native <- predict(object, data = changed, exdat = nd["x"],
      ezdat = nd["z"], newdata = data.frame(invalid = 1), se.fit = TRUE)
    expect_equal(lapply(formula.native, unname), oracle, tolerance = 0)
  }
})

test_that("plreg prediction SE consumes each prepared formula role once", {
  old <- options(np.messages = FALSE, na.action = "na.omit")
  on.exit(options(old), add = TRUE)
  d <- plreg_rows_fixture()
  hits <- 0L
  bump <- function(x) { hits <<- hits + 1L; x }
  bw <- npplregbw(y ~ bump(x) | z, data = d,
    bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  fit <- npplreg(bw)
  hits <- 0L
  ans <- predict(fit, se.fit = TRUE)
  expect_identical(hits, 0L)
  expect_true(all(is.finite(ans$se.fit)))
})

test_that("plreg prediction SE rejects malformed maps without residual value filtering", {
  old <- options(np.messages = FALSE, na.action = "na.omit")
  on.exit(options(old), add = TRUE)
  d <- plreg_rows_fixture()
  bw <- plreg_rows_bw(d)
  fit <- npplreg(bw)
  ns <- asNamespace(getNamespaceName(environment(npplreg)))
  prepare <- get(".np_plreg_predict_se_data", ns)
  infer <- get(".np_plreg_predict_se", ns)
  broken <- fit
  broken$call$txdat <- NULL
  expect_error(prepare(bw, broken), "call does not contain 'txdat'", fixed = TRUE)
  broken <- fit
  broken$eval.keep <- TRUE
  expect_error(infer(bw, broken), "evaluation row map does not match evaluation data", fixed = TRUE)
  empty <- fit
  empty$eval.keep <- rep(FALSE, nrow(d))
  expect_identical(infer(bw, empty), rep(NA_real_, nrow(d)))
  options(na.action = "na.exclude")
  nd <- data.frame(x = c(.2, NA, .7), z = c(.3, .4, .6))
  padded <- npplreg(bw, newdata = nd)
  padded$eval.keep <- rep(FALSE, nrow(padded$evalx))
  expect_identical(infer(bw, padded), rep(NA_real_, nrow(nd)))
  original <- get("npplreg", ns)
  contaminated <- function(...) {
    out <- original(...)
    if (isTRUE(out$residuals)) out$resid[2L] <- NA_real_
    out
  }
  local_mocked_bindings(npplreg = contaminated, .package = getNamespaceName(ns))
  expect_identical(infer(bw, fit), rep(NA_real_, nrow(d)))
})

test_that("plreg native-bandwidth prediction retains complete-case and no-NA controls", {
  old <- options(np.messages = FALSE, na.action = "na.omit")
  on.exit(options(old), add = TRUE)
  d <- plreg_rows_fixture()
  bw <- npplregbw(xdat = d["x"], ydat = d$y, zdat = d["z"],
    bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  fit <- npplreg(bw, residuals = TRUE, se = TRUE)
  rng <- .Random.seed
  ans <- predict(fit, se.fit = TRUE)
  oracle <- plreg_rows_oracle(bw, d, d[c("x", "z")])$pred
  expect_identical(.Random.seed, rng)
  expect_equal(lapply(ans, unname), oracle, tolerance = 0)
  expect_identical(predict(fit), ans$fit)
  d$y[3L] <- NA; d$x[8L] <- NA; d$z[14L] <- NA
  nd <- data.frame(x = c(.2, NA, .7), z = c(.3, .4, .6))
  ans <- predict(fit, txdat = d["x"], tydat = d$y, tzdat = d["z"],
    exdat = nd["x"], ezdat = nd["z"], se.fit = TRUE)
  expect_equal(lapply(ans, unname), plreg_rows_oracle(bw, d, nd)$pred, tolerance = 0)
})

test_that("plreg prediction SE composes response and predictor evaluation masks once", {
  old <- options(np.messages = FALSE, na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  d <- plreg_rows_fixture()
  d$y[3L] <- NA
  bw <- plreg_rows_bw(d)
  fit <- npplreg(bw)
  nd <- data.frame(y = c(0, 0, NA, 0, 0),
    x = c(.2, NA, .5, .7, .8), z = c(.3, .4, .6, NA, .5))
  expected <- lapply(plreg_rows_oracle(bw, d, nd[c(1L, 5L), ])$pred, function(x) {
    out <- rep(NA_real_, nrow(nd)); out[c(1L, 5L)] <- x; out
  })
  for (extra in list(
    list(newdata = nd, y.eval = TRUE),
    list(newdata = nd, eydat = c(0, NA, 0)),
    list(exdat = nd["x"], ezdat = nd["z"], eydat = nd$y,
      newdata = data.frame(invalid = 1)))) {
    ans <- do.call(predict, c(list(object = fit, se.fit = TRUE), extra))
    expect_equal(lapply(ans, unname), expected, tolerance = 0)
    expect_identical(which(is.na(ans$se.fit)), c(2L, 3L, 4L))
  }
})
