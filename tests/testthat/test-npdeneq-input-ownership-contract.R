di_payload <- function(x) x[c("Tn", "In", "Tn.bootstrap", "In.bootstrap", "Tn.P", "In.P", "boot.num")]

di_fixture <- function() {
  set.seed(91925)
  make <- function(n) data.frame(
    v = runif(n, .1, .9),
    u = factor(rep(c("a", "b"), length.out = n), levels = c("a", "b", "c")),
    o = ordered(rep(2:4, length.out = n), levels = 2:5))
  list(x = make(24L), y = make(21L))
}

test_that("density equality rejects nonsymmetric kernels without disabling density fits", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- di_fixture()
  x <- d$x["v"]; y <- d$y["v"]
  cases <- list(
    list(bwtype = "generalized_nn", bws = 7),
    list(bwtype = "adaptive_nn", bws = 7),
    list(ckerbound = "range", bws = .3),
    list(ckerbound = "fixed", ckerlb = 0, ckerub = 1, bws = .3),
    list(ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1, bws = .3))
  validate <- getFromNamespace(".npdeneq_validate_bandwidth", "np")
  for (spec in cases) {
    bw <- do.call(npudensbw, c(list(dat = x, bandwidth.compute = FALSE), spec))
    expect_true(all(is.finite(fitted(npudens(bws = bw, tdat = x)))))
    kbw <- getFromNamespace("kbandwidth", "np")(bw)
    for (object in list(bw, kbw)) {
      set.seed(123); seed <- .Random.seed
      expect_error(npdeneqtest(x, y, bw.x = object, B = 9),
                   "symmetric kernels without boundary normalization")
      expect_error(npdeneqtest(x, y, bw.y = object, B = 9),
                   "symmetric kernels without boundary normalization")
      expect_identical(.Random.seed, seed)
    }
    set.seed(123); seed <- .Random.seed
    expect_error(do.call(npdeneqtest, c(list(x = x, y = y, B = 9), spec)),
                 "symmetric kernels without boundary normalization")
    expect_identical(.Random.seed, seed)
  }
  # Retain compact symmetric kernels: a finite support is not a domain bound.
  uniform.call <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_gt(length(warnings), 0L)
    expect_true(all(grepl("ignoring kernel order specified with uniform kernel type",
                          warnings, fixed = TRUE)))
    value
  }
  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    bw <- npudensbw(dat = x, bws = .3, bandwidth.compute = FALSE,
                   ckertype = kernel)
    expect_identical(validate(bw), bw)
    if (kernel == "uniform") {
      # Existing kernel-sum conversion reports an ignored uniform order.
      fit <- uniform.call(npdeneqtest(x, y, bw.x = bw, B = 9))
      reverse <- uniform.call(npdeneqtest(y, x, bw.x = bw, B = 9))
    } else {
      fit <- npdeneqtest(x, y, bw.x = bw, B = 9)
      reverse <- npdeneqtest(y, x, bw.x = bw, B = 9)
    }
    expect_equal(fit$In, reverse$In, tolerance = 2e-12)
    expect_true(all(is.finite(c(fit$In, fit$Tn.bootstrap))))
  }
  expect_identical(validate(.3), .3)
  # The preselection seam must reject before invoking even a mocked selector.
  local_mocked_bindings(npudensbw = function(...) stop("SEARCH_ENTERED"),
                        .package = "np")
  select <- getFromNamespace(".npdeneq_select_bandwidth", "np")
  for (spec in cases)
    expect_error(do.call(select, c(list(x = x), spec)),
                 "symmetric kernels without boundary normalization")
  expect_error(select(x, bwtype = "g"), "symmetric kernels")
  expect_error(select(x, ckerbound = "r"), "symmetric kernels")
  expect_error(select(x, bwtype = "f", ckertype = "e", ckerbound = "n"),
               "SEARCH_ENTERED")
})


di_literal <- function(x, y, bx, by, B = 9L, seed = 42L) {
  statistic <- function(a, b) {
    sums <- function(dat, bw, exdat = NULL, power = 1L) {
      args <- list(txdat = dat, bws = bw, leave.one.out = is.null(exdat),
        bandwidth.divide = TRUE, kernel.pow = power)
      if (!is.null(exdat)) args$exdat <- exdat
      sum(do.call(npksum, args)$ksum)
    }
    n <- nrow(a); m <- nrow(b)
    In <- sums(a, bx)/(n*(n-1)) + sums(b, by)/(m*(m-1)) -
      2*sums(a, bx, b)/(n*m)
    variance <- 2*(sums(a, bx, power = 2L)/(n^2*(n-1)^2) +
      sums(b, by, power = 2L)/(m^2*(m-1)^2) +
      2*sums(a, bx, b, 2L)/(n^2*m^2))
    c(Tn = In/sqrt(variance), In = In)
  }
  set.seed(seed)
  pool <- rbind(x, y)
  draws <- vapply(seq_len(B), function(i) {
    a <- pool[sample.int(nrow(pool), nrow(x), TRUE), , drop = FALSE]
    b <- pool[sample.int(nrow(pool), nrow(y), TRUE), , drop = FALSE]
    statistic(a, b)
  }, numeric(2L))
  observed <- statistic(x, y)
  list(Tn = unname(observed[1L]), In = unname(observed[2L]),
    Tn.bootstrap = unname(draws[1L, ]), In.bootstrap = unname(draws[2L, ]),
    Tn.P = mean(draws[1L, ] >= observed[1L]),
    In.P = mean(draws[2L, ] >= observed[2L]))
}

test_that("density equality selects one first-sample common bandwidth", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- di_fixture()
  bx <- npudensbw(dat = d$x, bws = c(.3,.2,.2), bandwidth.compute = FALSE)
  by <- npudensbw(dat = d$y, bws = c(.3,.2,.2), bandwidth.compute = FALSE)
  selector <- getFromNamespace("npudensbw", "np")
  state <- new.env(parent = emptyenv())
  state$calls <- list()
  local_mocked_bindings(npudensbw = function(...) {
    args <- list(...)
    state$calls <- c(state$calls, list(args))
    selector(...)
  }, .package = "np")
  for (object in c(FALSE, TRUE)) {
    xbw <- if (object) bx else bx$bw
    ybw <- if (object) by else by$bw
    state$calls <- list()
    expected <- npdeneqtest(d$x, d$y, xbw, ybw, B = 9)
    for (side in c("x", "y")) {
      actual <- if (side == "x")
        npdeneqtest(d$x, d$y, bw.x = xbw, B = 9) else
        npdeneqtest(d$x, d$y, bw.y = ybw, B = 9)
      expect_identical(di_payload(actual), di_payload(expected))
    }
    expect_length(state$calls, 0L)
  }
  state$calls <- list()
  actual <- npdeneqtest(d$x, d$y, B = 9,
    bws = bx$bw, bandwidth.compute = FALSE)
  expect_length(state$calls, 1L)
  expect_identical(state$calls[[1L]]$dat, d$x)
  expect_identical(state$calls[[1L]]$bwmethod, "cv.ls")
  expected <- npdeneqtest(d$x, d$y, bx, bx, B = 9)
  expect_equal(di_payload(actual), di_payload(expected), tolerance = 2e-12)
  state$calls <- list()
  npdeneqtest(d$x, d$y, B = 9, bws = bx$bw,
    bandwidth.compute = FALSE, bwmethod = "cv.ml")
  expect_identical(state$calls[[1L]]$bwmethod, "cv.ml")
  expect_length(state$calls, 1L)
  set.seed(123); seed <- .Random.seed
  expect_error(npdeneqtest(d$x, d$y, .3, .4, B = 9), "one common bandwidth")
  altered <- bx; altered$ukertype <- "liracine"
  expect_error(npdeneqtest(d$x, d$y, bx, altered, B = 9), "one common bandwidth")
  expect_identical(.Random.seed, seed)
})

test_that("density equality owns complete rows within its kernel domain", {
  
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- di_fixture()
  d$x$v[1L] <- NA_real_; d$x$o[4L] <- NA
  d$y$u[2L] <- NA
  xc <- stats::na.omit(d$x); yc <- stats::na.omit(d$y)
  expect_identical(levels(xc$u), levels(d$x$u))
  for (bounded in c(FALSE, TRUE)) {
    args <- if (bounded) list(ckerbound = "fixed", ckerlb = 0, ckerub = 1) else list()
    bx <- do.call(npudensbw, c(list(dat = xc, bws = c(.3,.2,.2),
      bandwidth.compute = FALSE), args))
    by <- do.call(npudensbw, c(list(dat = yc, bws = c(.3,.2,.2),
      bandwidth.compute = FALSE), args))
    set.seed(71); before <- .Random.seed
    if (bounded) {
      expect_error(npdeneqtest(d$x, d$y, bx, by, B = 9),
                   "symmetric kernels without boundary normalization")
      expect_identical(.Random.seed, before)
      next
    }
    actual <- npdeneqtest(d$x, d$y, bx, by, B = 9)
    expect_identical(.Random.seed, before)
    expected <- npdeneqtest(xc, yc, bx, by, B = 9)
    expect_identical(di_payload(actual), di_payload(expected))
    oracle <- di_literal(xc, yc, bx, by)
    expect_equal(di_payload(actual)[names(oracle)], oracle, tolerance = 2e-12)
    # Row ordering may change summation roundoff, never the observed functional.
    permuted <- npdeneqtest(xc[nrow(xc):1, ], yc[nrow(yc):1, ], bx, by, B = 9)
    expect_equal(c(actual$Tn, actual$In),
      c(permuted$Tn, permuted$In), tolerance = 2e-12)
    # A single supplied bandwidth is reused on the already-clean samples.
    partial <- do.call(npdeneqtest, c(list(x = d$x, y = d$y,
      bw.x = bx, bws = by$bw, bandwidth.compute = FALSE, B = 9), args))
    expect_equal(di_payload(partial), di_payload(expected), tolerance = 2e-12)
  }
})

test_that("density equality rejects insufficient samples and remains usable", {
  
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- di_fixture()
  for (side in c("x", "y")) for (n in 0:1) {
    bad <- d
    bad[[side]]$v[seq_len(nrow(bad[[side]]) - n)] <- NA_real_
    set.seed(123); seed <- .Random.seed
    expect_error(npdeneqtest(bad$x, bad$y, c(.3,.2,.2), c(.3,.2,.2), B = 9),
      "at least two complete observations")
    expect_identical(.Random.seed, seed)
  }
  expect_s3_class(npdeneqtest(d$x, d$y, c(.3,.2,.2), c(.3,.2,.2), B = 9),
    "deneqtest")
})
