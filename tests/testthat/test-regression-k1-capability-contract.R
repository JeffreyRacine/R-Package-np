test_that("regression k1 capability is narrow and R-owned", {
  set.seed(20260817L)
  x <- data.frame(x = sort(runif(41L, -1, 1)))
  y <- sin(2 * x$x) + x$x^2

  make_bw <- function(kernel = "gaussian", order = 2L,
                      bwtype = "generalized_nn", regtype = "lc",
                      degree = NULL, basis = "glp", xdat = x, ydat = y,
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
      arguments$basis <- basis
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

  admitted <- list(
    list(label = "lc", regtype = "lc", degree = NULL),
    list(label = "lp0", regtype = "lp", degree = 0L),
    list(label = "lp1", regtype = "lp", degree = 1L),
    list(label = "lp3", regtype = "lp", degree = 3L),
    list(label = "lp12", regtype = "lp", degree = 12L)
  )
  for (kernel in c("gaussian", "epanechnikov")) {
    for (operator in admitted) {
      bw <- if (identical(operator$regtype, "lc")) {
        make_bw(kernel = kernel)
      } else {
        make_bw(
          kernel = kernel, regtype = operator$regtype,
          degree = operator$degree
        )
      }
      expect_identical(
        capability(bw),
        list(code = "gnn-univariate-positive-lp-k1", lower = 1L),
        info = paste(kernel, operator$label)
      )
      expect_identical(lower(bw), 1L)
      expect_identical(lower(bw, owner = "lsq"), 2L)
      expected.search.lower <- if (identical(operator$regtype, "lp")) 1L else 2L
      expect_identical(
        search_lower(bw, degree.candidates = list(0:12)),
        expected.search.lower
      )
      expect_identical(
        search_lower(bw, degree.candidates = list(0L)),
        expected.search.lower
      )

      H <- npreghat(
        bws = bw, txdat = x, leave.one.out = TRUE, output = "matrix"
      )
      applied <- npreghat(
        bws = bw, txdat = x, y = y, leave.one.out = TRUE, output = "apply"
      )
      objective <- np:::.npregbw_eval_only(x, y, bw)$objective
      expect_true(all(is.finite(H)))
      ridge <- attr(H, "ridge.used", exact = TRUE)
      expect_true(all(is.finite(ridge)))
      if (is.null(operator$degree) || operator$degree == 0L)
        expect_true(all(ridge == 0))
      expect_equal(as.numeric(H %*% y), as.numeric(applied), tolerance = 5e-9)
      expect_true(is.finite(objective))
      expect_equal(
        as.numeric(objective),
        mean((y - as.numeric(applied))^2),
        tolerance = 5e-9
      )

      setup <- setup_fun(xdat = x, template = bw, allow.extended.nn = FALSE)
      bounds <- bounds_fun(template = bw, setup = setup)
      expect_identical(as.numeric(bounds$lower[[1L]]), 1)
      degree.bounds <- bounds_fun(
        template = bw,
        setup = setup,
        degree.search = list(candidates = list(0:12))
      )
      expect_identical(
        as.numeric(degree.bounds$lower[[1L]]),
        as.numeric(expected.search.lower)
      )
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
    make_bw(regtype = "lp", degree = 1L, basis = "tensor", bandwidth = 2),
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

test_that("accepted k1 LP mean and derivative public owners agree", {
  set.seed(8171L)
  x <- data.frame(x = sort(runif(37L, -1.2, 1.3)))
  y <- cos(1.3 * x$x) + 0.2 * x$x
  ex <- data.frame(
    x = seq(min(x$x) + 1e-3, max(x$x) - 1e-3, length.out = 23L)
  )

  for (kernel in c("gaussian", "epanechnikov")) {
    for (degree in c(0L, 1L, 3L, 12L)) {
      bw <- npregbw(
        xdat = x, ydat = y, bws = 1,
        regtype = "lp", degree = degree, degree.select = "manual",
        basis = "glp", bernstein.basis = degree >= 3L,
        bwmethod = "cv.ls", bwtype = "generalized_nn",
        bwscaling = FALSE, ckertype = kernel, ckerorder = 2L,
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
        bws = bw, txdat = x, exdat = ex, y = y, s = 1L,
        output = "apply"
      )
      fit <- npreg(bws = bw, exdat = ex, gradients = TRUE)

      expect_equal(
        as.numeric(H %*% y), as.numeric(applied), tolerance = 5e-9,
        info = paste(kernel, degree, "mean apply")
      )
      expect_equal(
        drop(H %*% y), fit$mean, tolerance = 5e-9,
        info = paste(kernel, degree, "mean fit")
      )
      expect_equal(
        as.numeric(derivative_H %*% y),
        as.numeric(derivative_applied), tolerance = 5e-8,
        info = paste(kernel, degree, "derivative apply")
      )
      expect_equal(
        drop(derivative_H %*% y), fit$grad[, 1L], tolerance = 5e-8,
        info = paste(kernel, degree, "derivative fit")
      )
      expect_true(all(is.finite(attr(H, "ridge.used", exact = TRUE))))
      expect_true(all(is.finite(
        attr(derivative_H, "ridge.used", exact = TRUE)
      )))

      if (degree >= 2L) {
        for (operator in unique(c(2L, degree))) {
          operator_H <- npreghat(
            bws = bw, txdat = x, exdat = ex, s = operator,
            output = "matrix"
          )
          operator_applied <- npreghat(
            bws = bw, txdat = x, exdat = ex, y = y, s = operator,
            output = "apply"
          )
          expect_true(all(is.finite(operator_H)))
          expect_equal(
            as.numeric(operator_H %*% y), as.numeric(operator_applied),
            tolerance = 5e-8,
            info = paste(kernel, degree, "operator", operator)
          )
        }
      }
    }
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

  automatic.lp <- npregbw(
    xdat = data.frame(x = x), ydat = y,
    regtype = "lp", degree.select = "exhaustive",
    degree.min = 1L, degree.max = 1L,
    bwtype = "generalized_nn", ckertype = "gaussian",
    ckerorder = 2L, bwscaling = FALSE, nmulti = 1L,
    scale.factor.search.upper = 1
  )
  direct.lp <- npregbw(
    xdat = data.frame(x = x), ydat = y, bws = 1,
    regtype = "lp", degree = 1L, degree.select = "manual",
    basis = "glp", bernstein.basis = TRUE,
    bwmethod = "cv.ls", bwtype = "generalized_nn",
    ckertype = "gaussian", ckerorder = 2L,
    bwscaling = FALSE, bandwidth.compute = FALSE
  )
  expect_identical(as.integer(automatic.lp$bw[[1L]]), 1L)
  expect_identical(as.integer(automatic.lp$degree.engine[[1L]]), 1L)
  expect_equal(
    as.numeric(automatic.lp$fval),
    as.numeric(np:::.npregbw_eval_only(
      data.frame(x = x), y, direct.lp
    )$objective),
    tolerance = 5e-9
  )
  expect_true(all(is.finite(npreg(bws = automatic.lp, se = FALSE)$mean)))
})
