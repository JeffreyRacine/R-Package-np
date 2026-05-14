test_that("npreg lp higher-order gradients match npreghat across lp bases", {
  set.seed(99)
  n <- 140
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.02)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 35))

  for (basis in c("glp", "additive", "tensor")) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = 2L,
      basis = basis,
      bws = 0.25,
      bandwidth.compute = FALSE
    )

    fit <- suppressWarnings(npreg(
      txdat = tx,
      tydat = y,
      exdat = ex,
      bws = bw,
      gradients = TRUE,
      gradient.order = 2L
    ))

    H2 <- npreghat(bws = bw, txdat = tx, exdat = ex, s = 2L)
    g2 <- as.vector(H2 %*% y)

    expect_false(all(is.na(fit$grad[, 1])), info = basis)
    expect_equal(as.vector(fit$grad[, 1]), g2, tolerance = 1e-6, info = basis)
  }
})

test_that("npreg lp supports per-variable derivative orders", {
  set.seed(123)
  n <- 180
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + x2^3 + rnorm(n, sd = 0.03)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- data.frame(x1 = seq(0.05, 0.95, length.out = 20),
                   x2 = seq(0.1, 0.9, length.out = 20))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    degree = c(2L, 3L),
    basis = "glp",
    bws = c(0.2, 0.2),
    bandwidth.compute = FALSE
  )

  gorder <- c(2L, 1L)
  fit <- suppressWarnings(npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = gorder
  ))

  Hx1 <- npreghat(bws = bw, txdat = tx, exdat = ex, s = c(2L, 0L))
  Hx2 <- npreghat(bws = bw, txdat = tx, exdat = ex, s = c(0L, 1L))
  expect_equal(as.vector(fit$grad[, 1]), as.vector(Hx1 %*% y), tolerance = 1e-6)
  expect_equal(as.vector(fit$grad[, 2]), as.vector(Hx2 %*% y), tolerance = 1e-6)
})

test_that("npreg lp gradients accessor honors stored derivative order", {
  set.seed(20260514)
  n <- 55
  x <- seq(-1, 1, length.out = n)
  y <- 1 + 2 * x + 3 * x^2
  tx <- data.frame(x = x)

  bw2 <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    degree = 2L,
    basis = "glp",
    bws = 10,
    bandwidth.compute = FALSE
  )

  fit1 <- npreg(
    txdat = tx,
    tydat = y,
    bws = bw2,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )

  expect_error(
    gradients(fit1, gradient.order = 2L),
    "differs from the derivative order stored"
  )
  expect_error(
    gradients(fit1, errors = TRUE, gradient.order = 2L),
    "differs from the derivative order stored"
  )

  fit2 <- npreg(
    txdat = tx,
    tydat = y,
    bws = bw2,
    gradients = TRUE,
    gradient.order = 2L,
    warn.glp.gradient = FALSE
  )

  expect_equal(gradients(fit2, gradient.order = 2L), fit2$grad, tolerance = 0)
  expect_equal(gradients(fit2, errors = TRUE, gradient.order = 2L),
               fit2$gerr, tolerance = 0)

  bw1 <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    basis = "glp",
    bws = 10,
    bandwidth.compute = FALSE
  )
  fit.degree1 <- npreg(
    txdat = tx,
    tydat = y,
    bws = bw1,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )
  expect_warning(
    undefined <- gradients(fit.degree1, gradient.order = 2L),
    "exceed polynomial degree"
  )
  expect_true(all(is.na(undefined[, 1L])))
})

test_that("npreg lp Bernstein derivatives are returned on original scale", {
  n <- 50
  x <- seq(-1, 1, length.out = n)
  y <- 1 + 2 * x + 3 * x^2
  tx <- data.frame(x = x)

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    degree = 2L,
    basis = "tensor",
    bernstein.basis = TRUE,
    bws = 100,
    bandwidth.compute = FALSE
  )

  fit <- npreg(
    txdat = tx,
    tydat = y,
    bws = bw,
    gradients = TRUE,
    gradient.order = 2L,
    warn.glp.gradient = FALSE
  )

  expect_equal(as.vector(fitted(fit)), y, tolerance = 1e-8)
  expect_equal(as.vector(gradients(fit, gradient.order = 2L)[, 1L]),
               rep(6, n), tolerance = 1e-8)
  expect_true(all(is.finite(gradients(fit, errors = TRUE,
                                      gradient.order = 2L)[, 1L])))
})
