library(npRmpi)

bounded_scalar_kernel_nprmpi <- function(y.eval, y.train, h, lb, ub) {
  denom <- h * (pnorm((ub - y.eval) / h) - pnorm((lb - y.eval) / h))
  dnorm((y.eval - y.train) / h) / denom
}

test_that("bounded conditional cv.ls fixed-point objective matches direct delete-one quadrature", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260320)
  n <- 24L
  x <- runif(n)
  y <- rbeta(n, 1, 1)
  h <- c(0.15, 0.12)
  x.lb <- min(x)
  x.ub <- max(x)
  y.lb <- min(y)
  y.ub <- max(y)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = h,
    bwmethod = "cv.ls",
    regtype = "lc",
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range",
    cvls.quadrature.grid = "uniform"
  )

  Kx <- outer(seq_len(n), seq_len(n), Vectorize(function(i, j) {
    denom <- h[2] * (pnorm((x.ub - x[i]) / h[2]) - pnorm((x.lb - x[i]) / h[2]))
    dnorm((x[i] - x[j]) / h[2]) / denom
  }))
  Ky <- outer(seq_len(n), seq_len(n), Vectorize(function(i, j) {
    denom <- h[1] * (pnorm((y.ub - y[i]) / h[1]) - pnorm((y.lb - y[i]) / h[1]))
    dnorm((y[i] - y[j]) / h[1]) / denom
  }))

  w.full <- Kx / rowSums(Kx)
  w.loo <- Kx
  diag(w.loo) <- 0
  w.loo <- w.loo / rowSums(w.loo)

  f.loo <- rowSums(Ky * w.loo)
  q <- bw$cvls.quadrature.points[1]
  grid <- seq(y.lb, y.ub, length.out = q)
  trap <- diff(grid)[1] * c(0.5, rep(1, q - 2L), 0.5)
  Ky.grid <- outer(
    grid,
    y,
    Vectorize(function(y.eval, y.train) {
      bounded_scalar_kernel_nprmpi(y.eval, y.train, h[1], y.lb, y.ub)
    })
  )

  intsq <- vapply(
    seq_len(n),
    function(i) {
      fit <- as.numeric(Ky.grid %*% w.loo[i, ])
      sum(trap * fit^2)
    },
    numeric(1)
  )

  exact.obj <- -(mean(intsq) - 2 * mean(f.loo))
  np.obj <- npRmpi:::.npcdensbw_eval_only(data.frame(x = x), y, bw)$objective

  expect_equal(np.obj, exact.obj, tolerance = 1e-6)
})
