di_payload <- function(x) x[c("Tn", "In", "Tn.bootstrap", "In.bootstrap", "Tn.P", "In.P", "boot.num")]

di_fixture <- function() {
  set.seed(91925)
  make <- function(n) data.frame(
    v = runif(n, .1, .9),
    u = factor(rep(c("a", "b"), length.out = n), levels = c("a", "b", "c")),
    o = ordered(rep(2:4, length.out = n), levels = 2:5))
  list(x = make(24L), y = make(21L))
}

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
    Tn.P = mean(draws[1L, ] > observed[1L]),
    In.P = mean(draws[2L, ] > observed[2L]))
}

test_that("density equality selects only absent bandwidths", {
  
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- di_fixture()
  bx <- npudensbw(dat = d$x, bws = c(.3,.2,.2), bandwidth.compute = FALSE)
  by <- npudensbw(dat = d$y, bws = c(.4,.3,.3), bandwidth.compute = FALSE)
  selector <- getFromNamespace("npudensbw", "np")
  selected <- list()
  local_mocked_bindings(npudensbw = function(...) {
    args <- list(...)
    selected[[length(selected) + 1L]] <<- args$dat
    selector(...)
  }, .package = "np")
  for (object in c(FALSE, TRUE)) {
    xbw <- if (object) bx else bx$bw
    ybw <- if (object) by else by$bw
    selected <- list()
    expected <- npdeneqtest(d$x, d$y, xbw, ybw, B = 9)
    expect_length(selected, 0L)
    for (side in c("x", "y")) {
      # A selected density bandwidth carries density kernel normalization;
      # a numeric explicit bandwidth retains the kernel-sum defaults.
      expected <- if (side == "x")
        npdeneqtest(d$x, d$y, xbw, by, B = 9) else
        npdeneqtest(d$x, d$y, bx, ybw, B = 9)
      selected <- list()
      actual <- if (side == "x")
        npdeneqtest(d$x, d$y, bw.x = xbw, B = 9,
          bws = by$bw, bandwidth.compute = FALSE) else
        npdeneqtest(d$x, d$y, bw.y = ybw, B = 9,
          bws = bx$bw, bandwidth.compute = FALSE)
      expect_length(selected, 1L)
      expect_identical(selected[[1L]], if (side == "x") d$y else d$x)
      expect_equal(di_payload(actual), di_payload(expected), tolerance = 2e-12)
    }
  }
  selected <- list()
  actual <- npdeneqtest(d$x, d$y, B = 9,
    bws = bx$bw, bandwidth.compute = FALSE)
  expect_identical(selected, list(d$x, d$y))
  same.y <- selector(dat = d$y, bws = bx$bw, bandwidth.compute = FALSE)
  expected <- npdeneqtest(d$x, d$y, bx, same.y, B = 9)
  expect_equal(di_payload(actual), di_payload(expected), tolerance = 2e-12)
})

test_that("density equality owns complete rows before all bootstrap routes", {
  
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
    by <- do.call(npudensbw, c(list(dat = yc, bws = c(.4,.3,.3),
      bandwidth.compute = FALSE), args))
    set.seed(71); before <- .Random.seed
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
    # A missing bandwidth is selected on its own already-clean sample.
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
    expect_error(npdeneqtest(bad$x, bad$y, c(.3,.2,.2), c(.4,.3,.3), B = 9),
      "at least two complete observations")
    expect_identical(.Random.seed, seed)
  }
  expect_s3_class(npdeneqtest(d$x, d$y, c(.3,.2,.2), c(.4,.3,.3), B = 9),
    "deneqtest")
})
