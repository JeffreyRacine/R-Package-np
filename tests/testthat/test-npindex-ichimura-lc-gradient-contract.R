test_that("npRmpi npindex ichimura lc gradients preserve the scalar-index fit contract", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 1000L
  x1 <- runif(n, min = -1, max = 1)
  x2 <- runif(n, min = -1, max = 1)
  y <- x1 - x2 + rnorm(n)
  dat <- data.frame(x1 = x1, x2 = x2, y = y)

  bw <- npindexbw(formula = y ~ x1 + x2, data = dat)
  fit <- npindex(bws = bw, gradients = TRUE)

  index.df <- data.frame(index = as.vector(as.matrix(dat[c("x1", "x2")]) %*% bw$beta))
  oracle <- do.call(npreg, list(
    txdat = index.df,
    tydat = dat$y,
    exdat = index.df,
    bws = bw$bw,
    bwtype = bw$type,
    ckertype = bw$ckertype,
    ckerorder = bw$ckerorder,
    regtype = bw$regtype,
    gradients = TRUE,
    warn.glp.gradient = FALSE
  ))

  expect_true(all(is.finite(fit$grad)))
  expect_true(all(is.finite(fit$betavcov)))
  expect_equal(as.numeric(fit$mean), as.numeric(oracle$mean), tolerance = 0)
  expect_equal(
    unname(as.matrix(fit$grad)),
    unname(as.matrix(oracle$grad) %*% t(as.vector(bw$beta))),
    tolerance = 0
  )
})

test_that("npRmpi Ichimura covariance scales score columns by observation", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260626)
  n <- 250L
  x <- matrix(runif(n * 3L, min = -1, max = 1), nrow = n)
  colnames(x) <- c("x1", "x2", "x3")
  x <- as.data.frame(x)
  beta.true <- c(1, 0.5, -0.25)
  index.true <- as.vector(as.matrix(x) %*% beta.true)
  y <- sin(2 * index.true) + 0.5 * index.true + rnorm(n, sd = 0.35)

  fit <- npindex(
    txdat = x,
    tydat = y,
    method = "ichimura",
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "epanechnikov",
    optim.method = "BFGS",
    nmulti = 1L,
    random.seed = 20260626L,
    gradients = TRUE,
    residuals = TRUE
  )

  bw <- fit$bws
  index <- data.frame(index = as.vector(as.matrix(x) %*% bw$beta))
  link <- npreg(
    txdat = index,
    tydat = y,
    exdat = index,
    bws = bw$bw,
    bwtype = bw$type,
    ckertype = bw$ckertype,
    ckerorder = bw$ckerorder,
    regtype = bw$regtype,
    gradients = TRUE,
    warn.glp.gradient = FALSE
  )

  x.free <- as.matrix(x)[, -1L, drop = FALSE]
  weighted.sum <- npksum(
    txdat = index,
    tydat = rep(1, n),
    weights = x.free,
    bws = bw$bw,
    bwtype = bw$type,
    ckertype = bw$ckertype,
    ckerorder = bw$ckerorder
  )$ksum
  density.sum <- npksum(
    txdat = index,
    bws = bw$bw,
    bwtype = bw$type,
    ckertype = bw$ckertype,
    ckerorder = bw$ckerorder
  )$ksum

  centered <- vapply(
    seq_len(n),
    function(i) x.free[i, ] - weighted.sum[, i] / density.sum[i],
    numeric(2L)
  )
  link.gradient <- as.numeric(link$grad[, 1L])
  residual <- as.numeric(y - link$mean)

  score <- matrix(NA_real_, nrow = 2L, ncol = n)
  information <- meat <- matrix(0, nrow = 2L, ncol = 2L)
  for (i in seq_len(n)) {
    score[, i] <- link.gradient[i] * centered[, i]
    score.outer <- score[, i, drop = FALSE] %*%
      t(score[, i, drop = FALSE])
    information <- information + score.outer
    meat <- meat + residual[i]^2 * score.outer
  }
  information.inverse <- chol2inv(chol(information))
  oracle <- information.inverse %*% meat %*% information.inverse

  expect_equal(
    unname(fit$betavcov[-1L, -1L, drop = FALSE]),
    oracle,
    tolerance = 128 * .Machine$double.eps * n
  )
  expect_equal(unname(fit$betavcov[1L, ]), rep(0, 3L), tolerance = 0)
  expect_equal(unname(fit$betavcov[, 1L]), rep(0, 3L), tolerance = 0)
})
