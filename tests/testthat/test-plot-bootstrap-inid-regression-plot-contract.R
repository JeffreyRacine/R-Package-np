test_that("rbandwidth plot bootstrap supports gradients across methods", {
  skip_if_not_installed("np")

  set.seed(322)
  n <- 60
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.15)
  bw <- npregbw(y ~ x, regtype = "ll", nmulti = 1)

  run_plot <- function(method) {
    suppressWarnings(
      plot(
        bw,
        output = "data",
        perspective = FALSE,
        gradients = TRUE,
        errors = "bootstrap",
        bootstrap = method,
        B = 9
      )
    )
  }

  methods <- c("inid", "fixed", "geom", "wild")
  for (m in methods) {
    set.seed(9320 + match(m, methods))
    out <- run_plot(method = m)
    expect_type(out, "list")
    expect_true(length(out) > 0)
  }
})

test_that("plot bootstrap accepts wild selector", {
  skip_if_not_installed("np")

  set.seed(323)
  n <- 80
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  bw <- npregbw(y ~ x, bws = 0.2, bandwidth.compute = FALSE)

  out <- suppressWarnings(
    plot(
      bw,
      xdat = data.frame(x = x),
      ydat = y,
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      bootstrap = "wild",
      B = 9
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("scbandwidth bootstrap non-coef path avoids npscoef refits", {
  skip_if_not_installed("np")

  set.seed(3240)
  n <- 65
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)
  xdat <- data.frame(x = x)
  zdat <- data.frame(z = z)
  bw <- npscoefbw(xdat = xdat, ydat = y, zdat = zdat, regtype = "lc", nmulti = 1)

  np.ns <- asNamespace("np")
  ctr <- new.env(parent = emptyenv())
  ctr$n <- 0L
  trace(
    what = "npscoef",
    where = np.ns,
    tracer = bquote(.(ctr)$n <- .(ctr)$n + 1L),
    print = FALSE
  )
  on.exit(untrace("npscoef", where = np.ns), add = TRUE)

  out <- suppressWarnings(
    plot(
      bw,
      xdat = xdat,
      ydat = y,
      zdat = zdat,
      coef = FALSE,
      perspective = FALSE,
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 7
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
  expect_identical(ctr$n, 0L)
})

test_that("inid lc chunked generation matches explicit chunked draws", {
  skip_if_not_installed("np")

  set.seed(324)
  n <- 30
  x <- runif(n)
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 17))
  bw <- npregbw(y ~ x, regtype = "lc", bws = c(0.2), bandwidth.compute = FALSE)
  H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  B <- 11L

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  old.chunk <- getOption("np.plot.inid.chunk.size")
  on.exit(options(np.plot.inid.chunk.size = old.chunk), add = TRUE)
  options(np.plot.inid.chunk.size = 3L)

  set.seed(7001)
  auto <- fast.fun(H = H, ydat = y, B = B, counts = NULL)

  set.seed(7001)
  counts <- matrix(NA_integer_, nrow = n, ncol = B)
  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + 3L - 1L)
    bsz <- stopi - start + 1L
    counts[, start:stopi] <- rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    start <- stopi + 1L
  }
  explicit <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  expect_equal(auto$t, explicit$t, tolerance = 1e-12)
  expect_equal(auto$t0, explicit$t0, tolerance = 1e-12)
})

test_that("lp degree is used by regression inid helper construction", {
  skip_if_not_installed("np")

  set.seed(3264)
  n <- 75
  x1 <- runif(n)
  x2 <- runif(n)
  tx <- data.frame(x1 = x1, x2 = x2)
  y <- sin(3 * x1) + 0.25 * x2^2 + rnorm(n, sd = 0.04)
  ex <- tx[seq_len(20), , drop = FALSE]
  B <- 6L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  bw2 <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(2L, 2L),
    bws = c(0.28, 0.28),
    bandwidth.compute = FALSE
  )
  bw3 <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(3L, 3L),
    bws = c(0.28, 0.28),
    bandwidth.compute = FALSE
  )

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
  out2 <- fast.fun(xdat = tx, exdat = ex, bws = bw2, ydat = y, B = B, counts = counts)
  out3 <- fast.fun(xdat = tx, exdat = ex, bws = bw3, ydat = y, B = B, counts = counts)

  fit2 <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw2, gradients = FALSE, warn.glp.gradient = FALSE)$mean
  fit3 <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw3, gradients = FALSE, warn.glp.gradient = FALSE)$mean

  expect_equal(out2$t0, fit2, tolerance = 1e-6)
  expect_equal(out3$t0, fit3, tolerance = 1e-6)
  expect_gt(max(abs(out2$t0 - out3$t0)), 1e-6)
})

test_that("npreg plot bootstrap inid supports lp basis variants", {
  skip_if_not_installed("np")

  set.seed(327)
  n <- 70
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)

  for (basis in c("additive", "tensor", "glp")) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = 2,
      basis = basis,
      bws = 0.2,
      bandwidth.compute = FALSE
    )

    out <- suppressWarnings(
      plot(
        bw,
        xdat = tx,
        ydat = y,
        output = "data",
        perspective = FALSE,
        errors = "bootstrap",
        bootstrap = "inid",
        B = 7
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0)
  }
})
