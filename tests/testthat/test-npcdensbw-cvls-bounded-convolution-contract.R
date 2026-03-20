library(npRmpi)

bounded_pair_convolution_exact_nprmpi <- function(y.eval, y.train, h, lb, ub) {
  integrate(
    function(t) {
      denom <- h * (pnorm((ub - t) / h) - pnorm((lb - t) / h))
      dnorm((t - y.eval) / h) / denom * dnorm((t - y.train) / h) / denom
    },
    lower = lb,
    upper = ub,
    subdivisions = 4000L,
    rel.tol = 1e-10
  )$value
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
  intsq <- vapply(seq_len(n), function(i) {
    fhat <- function(t) {
      vapply(t, function(tt) {
        denom <- h[1] * (pnorm((y.ub - tt) / h[1]) - pnorm((y.lb - tt) / h[1]))
        ky.t <- dnorm((tt - y) / h[1]) / denom
        sum(w.full[i, ] * ky.t)
      }, numeric(1))
    }
    integrate(
      function(t) fhat(t)^2,
      lower = y.lb,
      upper = y.ub,
      subdivisions = 4000L,
      rel.tol = 1e-10
    )$value
  }, numeric(1))

  exact.obj <- -(mean(intsq) - 2 * mean(f.loo))
  np.obj <- npRmpi:::.npcdensbw_eval_only(data.frame(x = x), y, bw)$objective

  expect_equal(np.obj, exact.obj, tolerance = 1e-6)
})
