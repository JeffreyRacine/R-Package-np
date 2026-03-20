library(np)

empty_shadow_matrix <- function(n) matrix(numeric(0), nrow = n, ncol = 0)

shadow_rbw_conditional <- function(bw) {
  c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )
}

shadow_safe_call <- function(name, ...) {
  on.exit(
    tryCatch(.Call("C_np_shadow_reset_state", PACKAGE = "np"),
             error = function(e) NULL),
    add = TRUE
  )
  .Call(name, ..., PACKAGE = "np")
}

bounded_pair_convolution_exact <- function(y.eval, y.train, h, lb, ub) {
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

test_that("bounded conditional y-convolution row matches direct integration at the upper boundary", {
  set.seed(20260320)
  n <- 28L
  x <- runif(n)
  y <- rbeta(n, 1, 1)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.15, 0.12),
    bwmethod = "cv.ls",
    regtype = "lc",
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )

  rbw <- as.double(shadow_rbw_conditional(bw))
  row_idx <- which.max(y)
  y.lb <- min(y)
  y.ub <- max(y)

  yconv.row <- shadow_safe_call(
    "C_np_shadow_cv_yrow_conditional",
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(data.frame(y = y)),
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(data.frame(x = x)),
    rbw, 0L, 0L, 0L, 0L, FALSE, as.integer(row_idx), 1L,
    as.double(y.lb), as.double(y.ub)
  )

  exact.row <- vapply(
    y,
    bounded_pair_convolution_exact,
    numeric(1),
    y.eval = y[row_idx],
    h = bw$ybw[1],
    lb = y.lb,
    ub = y.ub
  )

  expect_equal(as.numeric(yconv.row), exact.row, tolerance = 1e-6)
})

test_that("bounded conditional cv.ls fixed-point objective matches direct integration", {
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
  np.obj <- np:::.npcdensbw_eval_only(data.frame(x = x), y, bw)$objective

  expect_equal(np.obj, exact.obj, tolerance = 1e-6)
})
