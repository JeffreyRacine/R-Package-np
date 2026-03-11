test_that("npreg lp higher-order gradients match npreghat across lp bases", {
  skip_if_not(spawn_mpi_slaves(1))
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old.opts), add = TRUE)

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
  skip_if_not(spawn_mpi_slaves(1))
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old.opts), add = TRUE)

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
