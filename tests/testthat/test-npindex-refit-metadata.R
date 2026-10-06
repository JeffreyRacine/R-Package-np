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

test_that("beta-only calls retain degree and use the supplied zero-tail start", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-2, 2, length.out = 80), z = sin(seq_len(80)))
  y <- sin(x$x + .3*x$z)
  seen <- new.env(parent = emptyenv()); seen$starts <- list()
  base.optim <- get("optim", envir = asNamespace("npRmpi"))
  testthat::local_mocked_bindings(
    .npindexbw_nomad_search = function(...) stop("unexpected NOMAD degree search"),
    .np_degree_search = function(...) stop("unexpected cell degree search"),
    lm = function(...) stop("unexpected OLS initialization"),
    optim = function(par, ...) {
      seen$starts[[length(seen$starts) + 1L]] <- par
      base.optim(par, ...)
    }, .package = "npRmpi")
  for (engine in c("nomad", "nomad+powell", "cell")) {
    b <- npindexbw(xdat = x, ydat = y, bws = c(1, 0, .8),
      regtype = "lp", degree = 1L, bernstein.basis = FALSE,
      only.optimize.beta = TRUE, nmulti = 1L, degree.select = "exhaustive",
      search.engine = engine, degree.min = 0L, degree.max = 2L)
    expect_identical(b$degree, 1L)
    expect_identical(b$bw, .8)
    expect_false(b$bernstein.basis)
    expect_null(b$degree.search)
    expect_true(is.finite(b$fval))
  }
  expect_length(seen$starts, 3L)
  expect_true(all(vapply(seen$starts, function(p) identical(as.double(p), 0), logical(1))))
  shortcut <- npindexbw(xdat = x, ydat = y, bws = c(1, 0, .8),
    regtype = "lp", degree = 1L, only.optimize.beta = TRUE, nmulti = 1L,
    nomad = TRUE, nomad.nmulti = 1L)
  expect_identical(shortcut$degree, 1L)
  expect_null(shortcut$nomad.shortcut)
  expect_error(npindexbw(xdat = x, ydat = y, bws = c(1, 0, .8),
    regtype = "lp", only.optimize.beta = TRUE, degree.select = "exhaustive"),
    "degree must be supplied explicitly")
})
