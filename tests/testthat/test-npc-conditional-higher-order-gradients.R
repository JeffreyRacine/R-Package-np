library(np)

test_that("npcdens exposes lp higher-order continuous gradients via hat operators", {
  npcdenshat <- getFromNamespace("npcdenshat", "np")

  set.seed(20260514)
  n <- 70
  x <- data.frame(x = seq(-1, 1, length.out = n))
  y <- data.frame(y = 1 + x$x + x$x^2 + 0.05 * sin(3 * x$x))
  ex <- data.frame(x = seq(-0.8, 0.8, length.out = 18))
  ey <- data.frame(y = seq(min(y$y) + 0.05, max(y$y) - 0.05, length.out = 18))
  rhs <- rep.int(1.0, n)

  for (basis in c("glp", "additive", "tensor")) {
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      basis = basis,
      degree = 2L,
      bws = c(0.35, 0.45),
      bandwidth.compute = FALSE
    )

    fit <- npcdens(
      bws = bw,
      txdat = x,
      tydat = y,
      exdat = ex,
      eydat = ey,
      gradients = TRUE,
      gradient.order = 2L
    )
    H2 <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, s = 2L)

    expect_equal(as.vector(fit$congrad[, 1L]), as.vector(H2 %*% rhs),
                 tolerance = 1e-9, info = basis)
    expect_true(all(is.na(fit$congerr[, 1L])), info = basis)
    expect_equal(gradients(fit, gradient.order = 2L), fit$congrad,
                 tolerance = 0, info = basis)
    expect_error(gradients(fit, gradient.order = 1L),
                 "differs from the derivative order stored",
                 info = basis)
  }
})

test_that("npcdist exposes lp higher-order continuous gradients via hat operators", {
  npcdisthat <- getFromNamespace("npcdisthat", "np")

  set.seed(20260514)
  n <- 70
  x <- data.frame(x = seq(-1, 1, length.out = n))
  y <- data.frame(y = 1 + x$x + x$x^2 + 0.05 * sin(3 * x$x))
  ex <- data.frame(x = seq(-0.8, 0.8, length.out = 18))
  ey <- data.frame(y = seq(min(y$y) + 0.05, max(y$y) - 0.05, length.out = 18))
  rhs <- rep.int(1.0, n)

  for (basis in c("glp", "additive", "tensor")) {
    bw <- npcdistbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      basis = basis,
      degree = 2L,
      bws = c(0.35, 0.45),
      bandwidth.compute = FALSE
    )

    fit <- npcdist(
      bws = bw,
      txdat = x,
      tydat = y,
      exdat = ex,
      eydat = ey,
      gradients = TRUE,
      gradient.order = 2L
    )
    H2 <- npcdisthat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, s = 2L)

    expect_equal(as.vector(fit$congrad[, 1L]), as.vector(H2 %*% rhs),
                 tolerance = 1e-9, info = basis)
    expect_true(all(is.na(fit$congerr[, 1L])), info = basis)
    expect_equal(gradients(fit, gradient.order = 2L), fit$congrad,
                 tolerance = 0, info = basis)
    expect_error(gradients(fit, gradient.order = 1L),
                 "differs from the derivative order stored",
                 info = basis)
  }
})

test_that("conditional higher-order gradient requests validate scope", {
  set.seed(20260514)
  x <- data.frame(x = seq(0, 1, length.out = 24))
  y <- data.frame(y = x$x + rnorm(24, sd = 0.01))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE
  )
  expect_warning(
    fit.ll <- npcdens(bws = bw.ll, txdat = x, tydat = y, gradients = TRUE,
                      gradient.order = 2L),
    "exceed polynomial degree"
  )
  expect_warning(
    expect_true(all(is.na(gradients(fit.ll, gradient.order = 2L)[, 1L]))),
    "exceed polynomial degree"
  )

  bw.lp1 <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE
  )
  expect_warning(
    fit <- npcdist(bws = bw.lp1, txdat = x, tydat = y, gradients = TRUE,
                   gradient.order = 2L),
    "exceed polynomial degree"
  )
  expect_warning(
    expect_true(all(is.na(gradients(fit, gradient.order = 2L)[, 1L]))),
    "exceed polynomial degree"
  )
})

test_that("formula routes keep gradient.order out of bandwidth selection", {
  set.seed(20260514)
  dat <- data.frame(
    x = seq(0, 1, length.out = 28),
    y = seq(0, 1, length.out = 28)^2 + rnorm(28, sd = 0.02)
  )

  fit.dens <- npcdens(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 2L,
    gradients = TRUE,
    gradient.order = 2L,
    nmulti = 1L
  )
  expect_equal(fit.dens$gradient.order, 2L)

  fit.dist <- npcdist(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 2L,
    gradients = TRUE,
    gradient.order = 2L,
    nmulti = 1L
  )
  expect_equal(fit.dist$gradient.order, 2L)
})

test_that("conditional plot data honors higher-order gradient requests", {
  set.seed(20260514)
  n <- 30
  x <- data.frame(x = seq(-1, 1, length.out = n))
  y <- data.frame(y = 1 + x$x + x$x^2)

  for (constructor in list(npcdensbw, npcdistbw)) {
    bw <- constructor(
      xdat = x,
      ydat = y,
      regtype = "lp",
      degree = 2L,
      basis = "tensor",
      bernstein.basis = FALSE,
      bws = c(10, 10),
      bandwidth.compute = FALSE
    )

    plot.out <- plot(
      bw,
      xdat = x,
      ydat = y,
      gradients = TRUE,
      gradient.order = 2L,
      xq = 0.5,
      yq = 0.5,
      neval = 8,
      perspective = FALSE,
      view = "fixed",
      plot.behavior = "data",
      plot.errors.method = "none"
    )

    slice <- plot.out[[1L]]
    eval.fun <- if (inherits(bw, "condbandwidth")) npcdist else npcdens
    fit <- eval.fun(
      bws = bw,
      txdat = x,
      tydat = y,
      exdat = slice$xeval,
      eydat = slice$yeval,
      gradients = TRUE,
      gradient.order = 2L
    )

    expect_equal(slice$gradient.order, 2L)
    expect_equal(as.vector(slice$congrad[, 1L]), as.vector(fit$congrad[, 1L]),
                 tolerance = 1e-10)
    expect_true(all(is.na(slice$congerr[, 1L])))
  }
})
