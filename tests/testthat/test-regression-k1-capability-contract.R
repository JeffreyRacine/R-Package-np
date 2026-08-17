test_that("regression k1 capability is narrow and R-owned", {
  set.seed(20260817L)
  x <- data.frame(x = sort(runif(41L, -1, 1)))
  y <- sin(2 * x$x) + x$x^2

  make_bw <- function(kernel = "gaussian", order = 2L,
                      bwtype = "generalized_nn", regtype = "lc",
                      degree = NULL, xdat = x, ydat = y,
                      ckerbound = "none", bandwidth = NULL) {
    arguments <- list(
      xdat = xdat,
      ydat = ydat,
      bws = if (is.null(bandwidth)) rep.int(1, ncol(xdat)) else bandwidth,
      bwmethod = "cv.ls",
      bwtype = bwtype,
      bwscaling = FALSE,
      ckertype = kernel,
      ckerorder = order,
      ckerbound = ckerbound,
      bandwidth.compute = FALSE,
      regtype = regtype
    )
    if (identical(regtype, "lp")) {
      arguments$degree <- rep.int(as.integer(degree), ncol(xdat))
      arguments$degree.select <- "manual"
      arguments$basis <- "glp"
    }
    do.call(npregbw, arguments)
  }

  capability <- getFromNamespace("npRegressionNnCapability", "np")
  lower <- getFromNamespace("npRegressionNnLowerBound", "np")
  prepare <- getFromNamespace(".npregbw_nomad_native_prepare_args", "np")
  setup_fun <- getFromNamespace(".npregbw_nomad_bw_setup", "np")
  bounds_fun <- getFromNamespace(".npregbw_nomad_bw_bounds", "np")
  search_lower <- getFromNamespace("npRegressionNnSearchLowerBound", "np")

  density.bw <- npudensbw(
    dat = x, bws = 1, bwtype = "generalized_nn",
    bandwidth.compute = FALSE
  )
  expect_identical(
    capability(density.bw),
    list(code = "not-regression", lower = 1L)
  )

  for (kernel in c("gaussian", "epanechnikov")) {
    for (operator in c("lc", "lp0")) {
      bw <- if (identical(operator, "lc")) {
        make_bw(kernel = kernel)
      } else {
        make_bw(kernel = kernel, regtype = "lp", degree = 0L)
      }
      expect_identical(
        capability(bw),
        list(code = "gnn-univariate-positive-lc-k1", lower = 1L)
      )
      expect_identical(lower(bw), 1L)
      expect_identical(lower(bw, owner = "lsq"), 2L)
      expect_identical(search_lower(bw, degree.candidates = list(0:2)), 2L)
      expect_identical(search_lower(bw, degree.candidates = list(0L)), 1L)

      H <- npreghat(
        bws = bw, txdat = x, leave.one.out = TRUE, output = "matrix"
      )
      applied <- npreghat(
        bws = bw, txdat = x, y = y, leave.one.out = TRUE, output = "apply"
      )
      objective <- np:::.npregbw_eval_only(x, y, bw)$objective
      expect_true(all(is.finite(H)))
      expect_true(all(attr(H, "ridge.used", exact = TRUE) == 0))
      expect_equal(as.numeric(H %*% y), as.numeric(applied), tolerance = 2e-13)
      expect_true(is.finite(objective))
      expect_equal(
        as.numeric(objective),
        mean((y - as.numeric(applied))^2),
        tolerance = 2e-13
      )

      setup <- setup_fun(xdat = x, template = bw, allow.extended.nn = FALSE)
      bounds <- bounds_fun(template = bw, setup = setup)
      expect_identical(as.numeric(bounds$lower[[1L]]), 1)
      degree.bounds <- bounds_fun(
        template = bw,
        setup = setup,
        degree.search = list(candidates = list(0:2))
      )
      expect_identical(as.numeric(degree.bounds$lower[[1L]]), 2)
      prep <- prepare(x, y, bw, invalid.penalty = "baseline")
      expect_identical(unname(tail(prep$myopti, 1L)), 1L)
    }
  }

  p2 <- data.frame(x1 = x$x, x2 = rev(x$x) + seq_len(nrow(x)) * 1e-8)
  mixed <- data.frame(x = x$x, u = factor(rep(c("a", "b"), length.out = nrow(x))))
  rejected <- list(
    make_bw(bwtype = "adaptive_nn", bandwidth = 2),
    make_bw(xdat = p2, ydat = y, bandwidth = c(2, 2)),
    make_bw(regtype = "ll", bandwidth = 2),
    make_bw(regtype = "lp", degree = 2L, bandwidth = 2),
    suppressWarnings(make_bw(kernel = "uniform", bandwidth = 2)),
    make_bw(kernel = "gaussian", order = 4L, bandwidth = 2),
    make_bw(kernel = "gaussian", ckerbound = "range", bandwidth = 2),
    make_bw(xdat = mixed, ydat = y, bandwidth = c(2, 0.5))
  )
  for (bw in rejected) {
    expect_identical(lower(bw), 2L)
    expect_identical(capability(bw)$code, "k2-only")
  }
  forged <- rejected[[1L]]
  forged$bw[forged$icon] <- 1
  forged$bandwidth$x[forged$xdati$icon] <- 1
  expect_error(npreg(bws = forged), "nearest-neighbor bandwidth must be in [2,", fixed = TRUE)
})

test_that("accepted k1 mean and derivative public owners agree", {
  set.seed(8171L)
  x <- data.frame(x = sort(runif(37L, -1.2, 1.3)))
  y <- cos(1.3 * x$x) + 0.2 * x$x
  ex <- data.frame(x = seq(-1.1, 1.2, length.out = 23L) + 1e-6)

  for (kernel in c("gaussian", "epanechnikov")) {
    bw <- npregbw(
      xdat = x, ydat = y, bws = 1,
      regtype = "lc", bwmethod = "cv.ls",
      bwtype = "generalized_nn", bwscaling = FALSE,
      ckertype = kernel, ckerorder = 2L,
      bandwidth.compute = FALSE
    )
    H <- npreghat(bws = bw, txdat = x, exdat = ex, output = "matrix")
    applied <- npreghat(
      bws = bw, txdat = x, exdat = ex, y = y, output = "apply"
    )
    derivative_H <- npreghat(
      bws = bw, txdat = x, exdat = ex, s = 1L, output = "matrix"
    )
    derivative_applied <- npreghat(
      bws = bw, txdat = x, exdat = ex, y = y, s = 1L, output = "apply"
    )
    fit <- npreg(bws = bw, exdat = ex, gradients = TRUE)

    expect_equal(as.numeric(H %*% y), as.numeric(applied), tolerance = 2e-13)
    expect_equal(drop(H %*% y), fit$mean, tolerance = 2e-13)
    expect_equal(
      as.numeric(derivative_H %*% y),
      as.numeric(derivative_applied),
      tolerance = 2e-12
    )
    expect_equal(drop(derivative_H %*% y), fit$grad[, 1L], tolerance = 2e-12)
    expect_true(all(attr(H, "ridge.used", exact = TRUE) == 0))
    expect_true(all(attr(derivative_H, "ridge.used", exact = TRUE) == 0))
  }
})

test_that("k1 zero radii retain objective-penalty and fit-error symmetry", {
  x <- data.frame(x = c(0, 0, 0.2, 0.4, 0.7, 1.1))
  y <- seq_len(nrow(x))
  bw <- npregbw(
    xdat = x, ydat = y, bws = 1,
    regtype = "lc", bwmethod = "cv.ls",
    bwtype = "generalized_nn", bwscaling = FALSE,
    ckertype = "gaussian", ckerorder = 2L,
    bandwidth.compute = FALSE
  )

  expect_identical(
    as.numeric(np:::.npregbw_eval_only(x, y, bw)$objective),
    10 * mean((y - mean(y))^2)
  )
  expect_error(npreg(bws = bw), "zero literal radius", fixed = TRUE)
  expect_error(
    npreghat(bws = bw, txdat = x),
    "zero literal radius",
    fixed = TRUE
  )
})

test_that("zero-mass compact-kernel LOO rows fail explicitly", {
  set.seed(81702L)
  x <- data.frame(
    x1 = runif(47L, -1.2, 1.4),
    x2 = runif(47L, -1.2, 1.4)
  )
  y <- sin(2.1 * x$x1) + 0.23 * x$x2^2
  bw <- npregbw(
    xdat = x,
    ydat = y,
    bws = c(4, 4),
    regtype = "lc",
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    bwscaling = FALSE,
    ckertype = "epanechnikov",
    ckerorder = 2L,
    bandwidth.compute = FALSE
  )

  expect_true(is.finite(np:::.npregbw_eval_only(x, y, bw)$objective))
  expect_error(
    npreghat(bws = bw, txdat = x, leave.one.out = TRUE),
    "leave-one-out kernel row has zero effective mass",
    fixed = TRUE
  )
  expect_error(
    npreghat(
      bws = bw,
      txdat = x,
      y = y,
      leave.one.out = TRUE,
      output = "apply"
    ),
    "leave-one-out kernel row has zero effective mass",
    fixed = TRUE
  )
})

test_that("the newly admitted discrete search endpoint is evaluated", {
  set.seed(20260828L)
  x <- sort(runif(37L))
  y <- sin(4 * pi * x) + rnorm(length(x), sd = 0.04)
  bw <- npregbw(
    xdat = data.frame(x = x), ydat = y,
    regtype = "lc", bwtype = "generalized_nn",
    ckertype = "gaussian", ckerorder = 2L,
    bwscaling = FALSE, nmulti = 1L
  )
  direct.k1 <- npregbw(
    xdat = data.frame(x = x), ydat = y, bws = 1,
    regtype = "lc", bwtype = "generalized_nn",
    ckertype = "gaussian", ckerorder = 2L,
    bwscaling = FALSE, bandwidth.compute = FALSE
  )

  expect_identical(as.integer(bw$bw[[1L]]), 1L)
  expect_identical(
    as.numeric(bw$fval),
    as.numeric(np:::.npregbw_eval_only(data.frame(x = x), y, direct.k1)$objective)
  )
})
