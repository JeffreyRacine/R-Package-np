conditional_beta_bw_weights <- function(train, evaluation, bandwidth,
                                        bwtype, kertype = "beta",
                                        order = 2L,
                                        operator = "normal") {
  args <- list(
    txdat = data.frame(value = train),
    exdat = data.frame(value = evaluation),
    bws = bandwidth,
    bwtype = bwtype,
    ckertype = kertype,
    ckerorder = order,
    operator = operator,
    return.kernel.weights = TRUE
  )
  if (identical(kertype, "beta")) {
    args$ckerbound <- "fixed"
    args$ckerlb <- 0
    args$ckerub <- 1
  }
  do.call(npksum, args)$kw
}

conditional_beta_manual_bw <- function(x, y, hx, hy,
                                       method = "cv.ml",
                                       bwtype = "fixed",
                                       xkernel = "beta",
                                       ykernel = "beta",
                                       xorder = 2L,
                                       yorder = 2L,
                                       scaling = FALSE,
                                       quadrature.grid = "uniform",
                                       quadrature.points = c(21L, 11L)) {
  args <- list(
    xdat = data.frame(x = x), ydat = data.frame(y = y),
    bws = c(hy, hx), bandwidth.compute = FALSE,
    bwmethod = method, bwtype = bwtype, bwscaling = scaling,
    cxkertype = xkernel, cxkerorder = xorder,
    cykertype = ykernel, cykerorder = yorder,
    cvls.quadrature.grid = quadrature.grid,
    cvls.quadrature.points = quadrature.points
  )
  if (identical(xkernel, "beta")) {
    args$cxkerbound <- "fixed"
    args$cxkerlb <- 0
    args$cxkerub <- 1
  }
  if (identical(ykernel, "beta")) {
    args$cykerbound <- "fixed"
    args$cykerlb <- 0
    args$cykerub <- 1
  }
  do.call(npcdensbw, args)
}

test_that("conditional beta objectives agree with independent fixed-bandwidth oracles", {
  x <- c(0.02, 0.08, 0.19, 0.37, 0.62, 0.81, 0.96)
  y <- c(0.03, 0.14, 0.24, 0.43, 0.59, 0.79, 0.97)
  hx <- 0.17
  hy <- 0.15
  wx <- conditional_beta_bw_weights(x, x, hx, "fixed")
  wy <- conditional_beta_bw_weights(y, y, hy, "fixed")
  diagonal <- cbind(seq_along(x), seq_along(x))
  wx[diagonal] <- 0
  wy[diagonal] <- 0

  density.ml <- conditional_beta_manual_bw(x, y, hx, hy, "cv.ml")
  fitted.loo <- colSums(wx * wy) / colSums(wx)
  oracle.ml <- -sum(log(fitted.loo))
  native.ml <- npRmpi:::.npcdensbw_eval_only(
    data.frame(x = x), data.frame(y = y), density.ml
  )$objective
  expect_equal(native.ml, -oracle.ml, tolerance = 2e-10)

  grid <- seq(0, 1, length.out = 9L)
  wy.cdf <- conditional_beta_bw_weights(
    y, grid, hy, "fixed", operator = "integral"
  )
  oracle.dist <- 0
  for (i in seq_along(x)) {
    keep <- setdiff(seq_along(x), i)
    for (g in seq_along(grid)) {
      estimate <- sum(wx[keep, i] * wy.cdf[keep, g]) / sum(wx[keep, i])
      oracle.dist <- oracle.dist + ((y[i] <= grid[g]) - estimate)^2
    }
  }
  oracle.dist <- oracle.dist / (length(x) * length(grid))
  distribution.bw <- npcdistbw(
    xdat = data.frame(x = x), ydat = data.frame(y = y),
    bws = c(hy, hx), bandwidth.compute = FALSE,
    bwmethod = "cv.ls", bwtype = "fixed", bwscaling = FALSE,
    cxkertype = "beta", cxkerorder = 2,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerorder = 2,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1
  )
  native.dist <- npRmpi:::.npcdistbw_eval_only(
    data.frame(x = x), data.frame(y = y), bws = distribution.bw,
    gydat = data.frame(y = grid)
  )$objective
  expect_equal(native.dist, oracle.dist, tolerance = 2e-10)

  wy.train.cdf <- conditional_beta_bw_weights(
    y, y, hy, "fixed", operator = "integral"
  )
  oracle.train.dist <- 0
  for (i in seq_along(x)) {
    keep <- setdiff(seq_along(x), i)
    for (g in seq_along(y)) {
      if (g == i)
        next
      estimate <- sum(wx[keep, i] * wy.train.cdf[keep, g]) /
        sum(wx[keep, i])
      oracle.train.dist <- oracle.train.dist + ((y[i] <= y[g]) - estimate)^2
    }
  }
  oracle.train.dist <- oracle.train.dist / length(x)^2
  native.train.dist <- npRmpi:::.npcdistbw_eval_only(
    data.frame(x = x), data.frame(y = y), bws = distribution.bw,
    do.full.integral = TRUE
  )$objective
  expect_equal(native.train.dist, oracle.train.dist, tolerance = 2e-10)

  qgrid <- seq(0, 1, length.out = 21L)
  qweights <- rep(diff(qgrid)[1L], length(qgrid))
  qweights[c(1L, length(qweights))] <-
    qweights[c(1L, length(qweights))] / 2
  wy.grid <- conditional_beta_bw_weights(y, qgrid, hy, "fixed")
  oracle.ls <- 0
  for (i in seq_along(x)) {
    keep <- setdiff(seq_along(x), i)
    denominator <- sum(wx[keep, i])
    fit.grid <- colSums(wx[keep, i] * wy.grid[keep, , drop = FALSE]) /
      denominator
    fit.at.yi <- sum(wx[keep, i] * wy[keep, i]) / denominator
    oracle.ls <- oracle.ls + sum(qweights * fit.grid^2) - 2 * fit.at.yi
  }
  oracle.ls <- oracle.ls / length(x)
  density.ls <- conditional_beta_manual_bw(x, y, hx, hy, "cv.ls")
  native.ls <- npRmpi:::.npcdensbw_eval_only(
    data.frame(x = x), data.frame(y = y), density.ls
  )$objective
  expect_equal(native.ls, -oracle.ls, tolerance = 3e-10)
})

test_that("conditional beta CVLS honors the bounded quadrature grid controls", {
  x <- c(0.02, 0.07, 0.13, 0.26, 0.44, 0.67, 0.85, 0.97)
  y <- c(0.01, 0.03, 0.08, 0.19, 0.42, 0.7, 0.91, 0.99)
  evaluate <- function(grid) {
    bw <- conditional_beta_manual_bw(
      x, y, 0.18, 0.15, method = "cv.ls",
      quadrature.grid = grid,
      quadrature.points = c(19L, 11L)
    )
    npRmpi:::.npcdensbw_eval_only(
      data.frame(x = x), data.frame(y = y), bw
    )$objective
  }
  objectives <- vapply(c("uniform", "hybrid", "sample"), evaluate, numeric(1L))

  expect_true(all(is.finite(objectives)))
  expect_gt(max(objectives) - min(objectives), 1e-10)
})

test_that("sample-grid conditional beta CVLS uses its actual generalized-NN grid", {
  x <- c(0.02, 0.07, 0.13, 0.26, 0.44, 0.67, 0.85, 0.97)
  y <- c(0.01, 0.03, 0.08, 0.19, 0.42, 0.7, 0.91, 0.99)
  hx <- 3
  hy <- 3
  grid <- sort(unique(y))
  weights <- c(
    0.5 * (grid[2L] - grid[1L]),
    0.5 * (grid[3:length(grid)] - grid[1:(length(grid) - 2L)]),
    0.5 * (grid[length(grid)] - grid[length(grid) - 1L])
  )
  wx <- conditional_beta_bw_weights(x, x, hx, "generalized_nn")
  wy.train <- conditional_beta_bw_weights(y, y, hy, "generalized_nn")
  wy.grid <- conditional_beta_bw_weights(y, grid, hy, "generalized_nn")
  oracle <- 0
  for (i in seq_along(x)) {
    keep <- setdiff(seq_along(x), i)
    denominator <- sum(wx[keep, i])
    fit.grid <- colSums(wx[keep, i] * wy.grid[keep, , drop = FALSE]) /
      denominator
    fit.at.yi <- sum(wx[keep, i] * wy.train[keep, i]) / denominator
    oracle <- oracle + sum(weights * fit.grid^2) - 2 * fit.at.yi
  }
  oracle <- oracle / length(x)
  bw <- conditional_beta_manual_bw(
    x, y, hx, hy, method = "cv.ls", bwtype = "generalized_nn",
    quadrature.grid = "sample", quadrature.points = c(19L, 11L)
  )
  native <- npRmpi:::.npcdensbw_eval_only(
    data.frame(x = x), data.frame(y = y), bw
  )$objective

  expect_equal(native, -oracle, tolerance = 5e-9)
})

test_that("conditional beta CVML agrees with nearest-neighbor weight ratios", {
  x <- c(0.02, 0.07, 0.16, 0.31, 0.5, 0.69, 0.84, 0.96)
  y <- c(0.03, 0.12, 0.22, 0.39, 0.57, 0.72, 0.87, 0.98)

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    hx <- 3
    hy <- 3
    wx <- conditional_beta_bw_weights(x, x, hx, bwtype)
    wy <- conditional_beta_bw_weights(y, y, hy, bwtype)
    diagonal <- cbind(seq_along(x), seq_along(x))
    wx[diagonal] <- 0
    wy[diagonal] <- 0
    oracle <- -sum(log(colSums(wx * wy) / colSums(wx)))
    bw <- conditional_beta_manual_bw(
      x, y, hx, hy, method = "cv.ml", bwtype = bwtype
    )
    native <- npRmpi:::.npcdensbw_eval_only(
      data.frame(x = x), data.frame(y = y), bw
    )$objective
    expect_equal(native, -oracle, tolerance = 5e-9)
  }
})

test_that("conditional beta scaling uses the matching X and Y standard deviations", {
  x <- c(0.02, 0.05, 0.11, 0.22, 0.41, 0.65, 0.83, 0.97)
  y <- c(0.2, 0.22, 0.27, 0.35, 0.48, 0.63, 0.81, 0.94)
  scale.x <- 0.8
  scale.y <- 1.1
  standard.deviation <- npRmpi:::EssDee(data.frame(x = x, y = y))
  normalization <- length(x)^(-1 / (2 * 2 + 2))
  hx <- scale.x * standard.deviation[1L] * normalization
  hy <- scale.y * standard.deviation[2L] * normalization
  wx <- conditional_beta_bw_weights(x, x, hx, "fixed")
  wy <- conditional_beta_bw_weights(y, y, hy, "fixed")
  diagonal <- cbind(seq_along(x), seq_along(x))
  wx[diagonal] <- 0
  wy[diagonal] <- 0
  oracle <- -sum(log(colSums(wx * wy) / colSums(wx)))
  bw <- conditional_beta_manual_bw(
    x, y, scale.x, scale.y, method = "cv.ml", scaling = TRUE
  )
  native <- npRmpi:::.npcdensbw_eval_only(
    data.frame(x = x), data.frame(y = y), bw
  )$objective

  expect_equal(native, -oracle, tolerance = 4e-9)
})

test_that("conditional beta objectives cover all orders and mixed sides", {
  x <- c(0.03, 0.09, 0.2, 0.38, 0.61, 0.78, 0.93)
  y <- c(0.04, 0.13, 0.27, 0.45, 0.58, 0.76, 0.95)

  for (order in c(2L, 4L, 6L, 8L)) {
    bw <- conditional_beta_manual_bw(
      x, y, 0.18, 0.16, method = "cv.ml",
      xorder = order, yorder = order
    )
    expect_true(is.finite(npRmpi:::.npcdensbw_eval_only(
      data.frame(x = x), data.frame(y = y), bw
    )$objective))
  }

  beta.x <- conditional_beta_manual_bw(
    x, qnorm(y), 0.18, 0.3, method = "cv.ml",
    xkernel = "beta", ykernel = "gaussian",
    xorder = 4L, yorder = 4L
  )
  beta.y <- conditional_beta_manual_bw(
    qnorm(x), y, 0.3, 0.16, method = "cv.ml",
    xkernel = "gaussian", ykernel = "beta",
    xorder = 4L, yorder = 4L
  )
  expect_true(is.finite(npRmpi:::.npcdensbw_eval_only(
    data.frame(x = x), data.frame(y = qnorm(y)), beta.x
  )$objective))
  expect_true(is.finite(npRmpi:::.npcdensbw_eval_only(
    data.frame(x = qnorm(x)), data.frame(y = y), beta.y
  )$objective))
})

test_that("automatic conditional beta selection returns usable bandwidths", {
  x <- data.frame(x = c(0.02, 0.06, 0.13, 0.24, 0.39, 0.57, 0.73, 0.86, 0.96))
  y <- data.frame(y = c(0.03, 0.1, 0.19, 0.3, 0.48, 0.61, 0.77, 0.89, 0.98))
  common <- list(
    xdat = x, ydat = y, bwscaling = FALSE,
    cxkertype = "beta", cxkerorder = 2,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerorder = 2,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1,
    nmulti = 1L, itmax = 8L, bwsolver = "powell"
  )
  density <- do.call(npcdensbw, c(common, list(
    bwmethod = "cv.ml", bwtype = "fixed",
    scale.factor.init = 0.18,
    scale.factor.init.lower = 0.06,
    scale.factor.init.upper = 0.4,
    scale.factor.search.lower = 0.03
  )))
  distribution <- do.call(npcdistbw, c(common, list(
    bwmethod = "cv.ls", bwtype = "generalized_nn", ngrid = 7L
  )))

  expect_true(all(is.finite(c(density$xbw, density$ybw, density$fval[1L]))))
  expect_true(all(c(density$xbw, density$ybw) > 0))
  expect_true(all(is.finite(c(
    distribution$xbw, distribution$ybw, distribution$fval[1L]
  ))))
  expect_equal(c(distribution$xbw, distribution$ybw),
               round(c(distribution$xbw, distribution$ybw)))

  skip_if_not_installed("crs", minimum_version = "0.15.46")
  mads.args <- common
  mads.args$bwmethod <- "cv.ml"
  mads.args$bwtype <- "fixed"
  mads.args$bwsolver <- "mads+powell"
  mads.args$nomad.opts <- list(MAX_BB_EVAL = 8)
  mads.args$scale.factor.init <- 0.18
  mads.args$scale.factor.init.lower <- 0.06
  mads.args$scale.factor.init.upper <- 0.4
  mads.args$scale.factor.search.lower <- 0.03
  mads <- do.call(npcdensbw, mads.args)
  expect_true(all(is.finite(c(mads$xbw, mads$ybw, mads$fval[1L]))))
})
