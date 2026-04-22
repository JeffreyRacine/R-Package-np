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

bounded_scalar_kernel <- function(y.eval, y.train, h, lb, ub) {
  denom <- h * (pnorm((ub - y.eval) / h) - pnorm((lb - y.eval) / h))
  dnorm((y.eval - y.train) / h) / denom
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

test_that("bounded conditional cv.ls fixed-point objective matches direct delete-one quadrature", {
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

  w.loo <- Kx
  diag(w.loo) <- 0
  w.loo <- w.loo / rowSums(w.loo)

  f.loo <- rowSums(Ky * w.loo)
  q <- 81L
  grid <- seq(y.lb, y.ub, length.out = q)
  trap <- diff(grid)[1] * c(0.5, rep(1, q - 2L), 0.5)
  Ky.grid <- outer(
    grid,
    y,
    Vectorize(function(y.eval, y.train) {
      bounded_scalar_kernel(y.eval, y.train, h[1], y.lb, y.ub)
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
  np.obj <- np:::.npcdensbw_eval_only(data.frame(x = x), y, bw)$objective

  expect_equal(np.obj, exact.obj, tolerance = 1e-6)
})
