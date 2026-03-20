library(npRmpi)

bounded_pair_convolution_exact_nprmpi <- function(y.eval, y.train, h, lb, ub) {
  gy.eval <- pnorm((ub - y.eval) / h) - pnorm((lb - y.eval) / h)
  gy.train <- pnorm((ub - y.train) / h) - pnorm((lb - y.train) / h)

  integrate(
    function(u) dnorm(u) * dnorm(u + (y.eval - y.train) / h),
    lower = (lb - y.eval) / h,
    upper = (ub - y.eval) / h,
    subdivisions = 4000L,
    rel.tol = 1e-10
  )$value / (h * gy.eval * gy.train)
}

test_that("bounded conditional cv.ls fixed-point objective matches direct integration", {
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
    cykerbound = "range"
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
  Ky.conv <- outer(
    y,
    y,
    Vectorize(function(y.eval, y.train) {
      bounded_pair_convolution_exact_nprmpi(y.eval, y.train, h[1], y.lb, y.ub)
    })
  )

  intsq <- vapply(
    seq_len(n),
    function(i) as.numeric(w.full[i, ] %*% Ky.conv %*% w.full[i, ]),
    numeric(1)
  )

  exact.obj <- -(mean(intsq) - 2 * mean(f.loo))
  np.obj <- npRmpi:::.npcdensbw_eval_only(data.frame(x = x), y, bw)$objective

  expect_equal(np.obj, exact.obj, tolerance = 1e-6)
})
