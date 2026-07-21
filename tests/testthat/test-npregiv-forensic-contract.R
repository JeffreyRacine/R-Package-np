iv_forensic_fixture <- function(n = 26L, seed = 20260721L) {
  set.seed(seed)
  w <- rnorm(n)
  v <- rnorm(n, sd = 0.22)
  z <- 0.45 * w + v
  y <- z^2 - 0.4 * v + rnorm(n, sd = 0.055)
  data.frame(y, z, w)
}

iv_forensic_skip_mpi_numerics <- function() {
  if ("npRmpi" %in% loadedNamespaces())
    skip("installed one- and multi-worker sentinels own npRmpi numerical coverage")
}

test_that("IV scalar controls fail early and deterministically", {
  dat <- iv_forensic_fixture(8L)
  expect_error(
    npregiv(y = dat$y, z = dat$z, w = dat$w,
            constant = c(0.25, 0.5)),
    "constant must be one finite number", fixed = TRUE
  )
  expect_error(
    npregiv(y = dat$y, z = dat$z, w = dat$w,
            iterate.max = 2.5),
    "iterate.max must be one integer", fixed = TRUE
  )
  expect_error(
    npregivderiv(y = dat$y, z = dat$z, w = dat$w,
                 stop.on.increase = c(TRUE, FALSE)),
    "stop.on.increase must be TRUE or FALSE", fixed = TRUE
  )
  expect_error(
    npregivderiv(y = dat$y, z = rep(1, nrow(dat)), w = dat$w),
    "at least two distinct support points", fixed = TRUE
  )
})

test_that("npregiv starting values and first-state replay are coherent", {
  iv_forensic_skip_mpi_numerics()
  dat <- iv_forensic_fixture()
  start <- seq(-0.15, 0.15, length.out = nrow(dat))

  started <- suppressWarnings(suppressMessages(npregiv(
    y = dat$y, z = dat$z, w = dat$w,
    starting.values = start, regtype = "lc",
    nmulti = 1L, iterate.max = 2L, stop.on.increase = FALSE
  )))
  expect_identical(started$starting.values.phi, start)
  expect_true(all(is.finite(started$phi.deriv.1)))

  fitted.path <- suppressWarnings(suppressMessages(npregiv(
    y = dat$y, z = dat$z, w = dat$w, regtype = "lc",
    nmulti = 1L, iterate.max = 3L, stop.on.increase = FALSE,
    iterate.diff.tol = 0
  )))
  replay.bw <- fitted.path
  replay.bw$norm.index <- 1L
  grid <- data.frame(z = seq(min(dat$z), max(dat$z), length.out = 7L))
  replay <- suppressWarnings(suppressMessages(npregiv(
    y = dat$y, z = dat$z, w = dat$w, zeval = grid,
    bw = replay.bw, regtype = "lc"
  )))

  expect_identical(replay$norm.index, 1L)
  expect_identical(replay$convergence, "BANDWIDTH_REPLAY")
  expect_identical(ncol(replay$phi.mat), 1L)
  expect_identical(length(replay$norm.stop), 1L)
  expect_identical(as.numeric(replay$phi),
                   as.numeric(fitted.path$phi.mat[, 1L]))
  expect_length(replay$phi.eval, nrow(grid))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(replay))
  expect_silent(plot(replay, TRUE, main = "IV fit",
                     xlim = range(dat$z), lty = 2))
  expect_silent(plot(replay, deriv = TRUE))
})

test_that("multivariate npregiv derivatives are coordinate-specific", {
  iv_forensic_skip_mpi_numerics()
  set.seed(812)
  n <- 38L
  w1 <- rnorm(n)
  w2 <- rnorm(n)
  v1 <- rnorm(n, sd = 0.12)
  v2 <- rnorm(n, sd = 0.12)
  z1 <- 0.7 * w1 + 0.1 * w2 + v1
  z2 <- -0.1 * w1 + 0.7 * w2 + v2
  y <- z1^2 + 0.5 * z2^2 - 0.25 * v1 + 0.1 * v2 + rnorm(n, sd = 0.04)
  z <- data.frame(z1, z2)
  w <- data.frame(w1, w2)

  orders <- getFromNamespace(".np_iv_partial_orders", "np")
  expect_identical(orders(2L, 1L), diag(1L, 2L, 2L))
  expect_identical(orders(2L, 2L), diag(2L, 2L, 2L))

  ll <- suppressWarnings(suppressMessages(npregiv(
    y = y, z = z, w = w,
    zeval = data.frame(z1 = 0, z2 = 0),
    regtype = "ll", nmulti = 1L, iterate.max = 2L,
    stop.on.increase = FALSE
  )))
  expect_identical(dim(ll$phi.deriv.1), c(n, 2L))
  expect_identical(colnames(ll$phi.deriv.1), names(z))
  expect_identical(dim(ll$phi.deriv.eval.1), c(1L, 2L))
  expect_identical(colnames(ll$phi.deriv.eval.1), names(z))
  expect_true(all(is.finite(ll$phi.deriv.1)))

  tikh <- suppressWarnings(suppressMessages(npregiv(
    y = y, z = z, w = w,
    method = "Tikhonov", regtype = "lc",
    alpha = 0.02, alpha.iter = 0.02, nmulti = 1L,
    return.weights.phi = TRUE,
    return.weights.phi.deriv.1 = TRUE
  )))
  expect_identical(dim(tikh$phi.deriv.1), c(n, 2L))
  expect_identical(colnames(tikh$phi.deriv.1), names(z))
  expect_named(tikh$phi.deriv.1.weights, names(z))
  for (j in seq_along(z)) {
    expect_equal(as.numeric(tikh$phi.deriv.1.weights[[j]] %*% y),
                 tikh$phi.deriv.1[, j], tolerance = 1e-11)
  }
})

test_that("npregiv owns categorical kernels and bandwidth scaling arguments", {
  iv_forensic_skip_mpi_numerics()
  dat <- iv_forensic_fixture(34L)
  wf <- factor(ifelse(dat$w > 0, "high", "low"))
  fit <- suppressWarnings(suppressMessages(npregiv(
    y = dat$y, z = dat$z, w = data.frame(wf),
    regtype = "lc", nmulti = 1L, iterate.max = 2L,
    ukertype = "aitchisonaitken", bandwidth.divide = FALSE
  )))
  upper <- getFromNamespace("uMaxL", "np")(nlevels(wf), "aitchisonaitken")
  expect_s3_class(fit, "npregiv")
  expect_lte(as.numeric(fit$bw.E.y.w[[1L]]), upper)
  expect_identical(fit$stage.specs$w$regtype, "lc")
})

test_that("npregivderiv evaluation is output-only and S3 accessors stay training-aligned", {
  iv_forensic_skip_mpi_numerics()
  dat <- iv_forensic_fixture()
  args <- list(nmulti = 1L, iterate.max = 3L,
               iterate.break = FALSE, stop.on.increase = FALSE,
               random.seed = 20260721L)
  set.seed(73)
  caller.seed <- .Random.seed
  training <- suppressWarnings(suppressMessages(do.call(
    npregivderiv, c(list(y = dat$y, z = dat$z, w = dat$w), args)
  )))
  expect_identical(.Random.seed, caller.seed)
  reversed <- suppressWarnings(suppressMessages(do.call(
    npregivderiv,
    c(list(y = dat$y, z = dat$z, w = dat$w, zeval = rev(dat$z)), args)
  )))
  grid <- data.frame(z = seq(min(dat$z), max(dat$z), length.out = 9L))
  evaluated <- suppressWarnings(suppressMessages(do.call(
    npregivderiv,
    c(list(y = dat$y, z = dat$z, w = dat$w, zeval = grid), args)
  )))
  formula.evaluated <- suppressWarnings(suppressMessages(do.call(
    npregivderiv,
    c(list(y ~ z | w, data = dat, newdata = grid), args)
  )))

  for (field in c("phi", "phi.prime", "phi.mat", "phi.prime.mat",
                  "norm.stop", "num.iterations")) {
    expect_identical(reversed[[field]], training[[field]])
    expect_identical(evaluated[[field]], training[[field]])
  }
  expect_identical(reversed$phi.prime.eval, rev(training$phi.prime))
  expect_identical(formula.evaluated$phi.prime.eval,
                   evaluated$phi.prime.eval)
  expect_identical(fitted(evaluated), evaluated$phi)
  expect_identical(gradients(evaluated), evaluated$phi.prime)
  expect_equal(residuals(evaluated), dat$y - evaluated$phi, tolerance = 0)
  expect_error(
    do.call(npregivderiv,
            c(list(y = dat$y, z = dat$z, w = dat$w,
                   weval = rev(dat$w)), args)),
    "weval cannot differ", fixed = TRUE
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(evaluated))
  expect_silent(plot(evaluated, phi = TRUE))
  expect_silent(plot(evaluated, plot.data = TRUE, phi = TRUE,
                     main = "IV derivative fit", lty = 2))
  expect_output(print(summary(evaluated)),
                "Bandwidth for E(y-phi(z)|w):", fixed = TRUE)
})
