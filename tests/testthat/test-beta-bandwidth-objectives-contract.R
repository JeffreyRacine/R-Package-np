beta_density_objective <- function(x, bandwidth, method = c("cv.ml", "cv.ls"),
                                   order = 2L, bwtype = "fixed",
                                   bwscaling = FALSE) {
  method <- match.arg(method)
  bws <- npudensbw(
    dat = x, bws = bandwidth, bandwidth.compute = FALSE,
    bwmethod = method, bwtype = bwtype, bwscaling = bwscaling,
    ckertype = "beta", ckerorder = order,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  out <- np:::npudensbw.bandwidth(
    dat = x, bws = bws, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L, invalid.penalty = "dbmax"
  )
  -as.numeric(out$fval[1L])
}

beta_distribution_objective <- function(x, bandwidth, order = 2L,
                                        bwtype = "fixed", gdat = NULL,
                                        do.full.integral = FALSE,
                                        ngrid = 19L) {
  bws <- npudistbw(
    dat = x, bws = bandwidth, bandwidth.compute = FALSE,
    bwtype = bwtype, ckertype = "beta", ckerorder = order,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  out <- np:::npudistbw.dbandwidth(
    dat = x, bws = bws, gdat = gdat,
    bandwidth.compute = TRUE, eval.only = TRUE, nmulti = 1L,
    do.full.integral = do.full.integral, ngrid = ngrid,
    invalid.penalty = "dbmax"
  )
  as.numeric(out$fval[1L])
}

beta_regression_objective <- function(x, y, bandwidth,
                                      method = c("cv.ls", "cv.aic"),
                                      order = 2L, bwtype = "fixed") {
  method <- match.arg(method)
  bws <- npregbw(
    xdat = x, ydat = y, bws = bandwidth, bandwidth.compute = FALSE,
    regtype = "lc", bwmethod = method, bwtype = bwtype,
    ckertype = "beta", ckerorder = order,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  as.numeric(np:::.npregbw_eval_only(
    xdat = x, ydat = y, bws = bws, invalid.penalty = "dbmax"
  )$objective[1L])
}

test_that("beta CVML and local-constant regression objectives match kernel weights", {
  xval <- c(0.002, 0.012, 0.045, 0.13, 0.31, 0.58, 0.79, 0.93, 0.995)
  x <- data.frame(x = xval)
  y <- sin(3 * xval) + xval^2

  for (order in c(2L, 4L, 6L, 8L)) {
    weights_loo <- npksum(
      txdat = x, bws = 0.14, leave.one.out = TRUE,
      return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )$kw
    density_loo <- colSums(weights_loo) / (nrow(x) - 1L)
    expected_ml <- sum(vapply(density_loo, function(fit) {
      if (fit > .Machine$double.xmin)
        return(-log(fit))
      if (fit < -.Machine$double.xmin)
        return(log(-fit) - 2 * log(.Machine$double.xmin))
      -log(.Machine$double.xmin)
    }, numeric(1L)))
    fitted_loo <- colSums(weights_loo * y) / colSums(weights_loo)
    expected_ls <- mean((y - fitted_loo)^2)

    expect_equal(
      beta_density_objective(x, 0.14, "cv.ml", order),
      expected_ml,
      tolerance = 2e-11
    )
    expect_equal(
      beta_regression_objective(x, y, 0.14, "cv.ls", order),
      expected_ls,
      tolerance = 2e-11
    )

    weights_full <- npksum(
      txdat = x, exdat = x, bws = 0.14,
      return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )$kw
    denominator <- colSums(weights_full)
    fitted_full <- colSums(weights_full * y) / denominator
    trace_hat <- sum(diag(weights_full) / denominator)
    loss <- mean((y - fitted_full)^2)
    expected_aic <- log(loss) +
      (1 + trace_hat / nrow(x)) / (1 - (trace_hat + 2) / nrow(x))
    expect_equal(
      beta_regression_objective(x, y, 0.14, "cv.aic", order),
      expected_aic,
      tolerance = 2e-10
    )
  }
})

test_that("beta density CVLS uses bounded target-coordinate quadrature", {
  x <- data.frame(x = c(0.003, 0.018, 0.07, 0.19, 0.43, 0.68, 0.87, 0.982))
  grid <- data.frame(x = seq(0, 1, length.out = 81L))
  step <- 1 / 80

  for (order in c(2L, 8L)) {
    grid_weights <- npksum(
      txdat = x, exdat = grid, bws = 0.16,
      return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )$kw
    density <- colMeans(grid_weights)
    integrated_square <- step *
      (sum(density^2) - 0.5 * density[1L]^2 - 0.5 * density[81L]^2)
    loo_weights <- npksum(
      txdat = x, bws = 0.16, leave.one.out = TRUE,
      return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )$kw
    cross_term <- mean(colSums(loo_weights) / (nrow(x) - 1L))
    expected <- integrated_square - 2 * cross_term

    expect_equal(
      beta_density_objective(x, 0.16, "cv.ls", order),
      expected,
      tolerance = 3e-10
    )
  }
})

test_that("beta CDF cross-validation matches explicit leave-one-out matrices", {
  xval <- c(0.002, 0.02, 0.09, 0.24, 0.51, 0.76, 0.94, 0.997)
  x <- data.frame(x = xval)
  grid <- data.frame(x = c(0, 0.01, 0.08, 0.3, 0.65, 0.9, 1))

  for (order in c(2L, 4L, 6L, 8L)) {
    weights <- npksum(
      txdat = x, exdat = grid, bws = 0.15,
      operator = "integral", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )$kw
    fitted_loo <-
      (matrix(colSums(weights), nrow(x), nrow(grid), byrow = TRUE) - weights) /
      (nrow(x) - 1L)
    indicator <- outer(xval, grid$x, "<=")
    expected <- mean((indicator - fitted_loo)^2)

    expect_equal(
      beta_distribution_objective(x, 0.15, order, gdat = grid),
      expected,
      tolerance = 3e-11
    )
  }

  train_weights <- npksum(
    txdat = x, exdat = x, bws = 0.15,
    operator = "integral", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )$kw
  fitted_loo <-
    (matrix(colSums(train_weights), nrow(x), nrow(x), byrow = TRUE) -
       train_weights) / (nrow(x) - 1L)
  squared <- (outer(xval, xval, "<=") - fitted_loo)^2
  expected_train <- (sum(squared) - sum(diag(squared))) / nrow(x)^2
  expect_equal(
    beta_distribution_objective(
      x, 0.15, 2L, do.full.integral = TRUE
    ),
    expected_train,
    tolerance = 3e-11
  )
})

test_that("beta objectives support fixed and both nearest-neighbor modes", {
  xval <- c(0.004, 0.015, 0.04, 0.11, 0.23, 0.39, 0.58, 0.74, 0.88, 0.97)
  x <- data.frame(x = xval)
  y <- cos(2.5 * xval) + 0.1 * xval

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(bwtype, "fixed")) 0.13 else 4
    expect_true(is.finite(beta_density_objective(
      x, bandwidth, "cv.ml", 4L, bwtype
    )))
    expect_true(is.finite(beta_density_objective(
      x, bandwidth, "cv.ls", 4L, bwtype
    )))
    expect_true(is.finite(beta_distribution_objective(
      x, bandwidth, 4L, bwtype, gdat = x
    )))
    expect_true(is.finite(beta_regression_objective(
      x, y, bandwidth, "cv.ls", 4L, bwtype
    )))
  }
})

test_that("beta fixed-bandwidth objectives preserve normalized scaling semantics", {
  x <- data.frame(x = c(0.006, 0.02, 0.08, 0.2, 0.41, 0.63, 0.81, 0.96))
  scale_factor <- 0.9
  bws <- npudensbw(
    dat = x, bws = scale_factor, bandwidth.compute = FALSE,
    bwmethod = "cv.ml", bwscaling = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  raw_bandwidth <- scale_factor * bws$sdev * bws$nconfac
  weights <- npksum(
    txdat = x, bws = raw_bandwidth, leave.one.out = TRUE,
    return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )$kw
  expected <- -sum(log(colSums(weights) / (nrow(x) - 1L)))

  expect_equal(
    beta_density_objective(
      x, scale_factor, "cv.ml", 2L, "fixed", bwscaling = TRUE
    ),
    expected,
    tolerance = 2e-11
  )
})

test_that("beta automatic selectors return usable Powell and native NOMAD objects", {
  set.seed(1909)
  xval <- rbeta(24, 0.7, 2.1)
  x <- data.frame(x = xval)
  y <- sin(4 * xval) + rnorm(length(xval), sd = 0.02)
  common <- list(
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1,
    nmulti = 1L
  )

  density <- do.call(npudensbw, c(list(
    dat = x, bwmethod = "cv.ml", bwsolver = "powell", itmax = 80L
  ), common))
  distribution <- do.call(npudistbw, c(list(
    dat = x, bwsolver = "powell", itmax = 80L, ngrid = 17L
  ), common))
  regression <- do.call(npregbw, c(list(
    xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ls",
    bwsolver = "powell", itmax = 80L
  ), common))

  expect_true(all(is.finite(c(density$bw, distribution$bw, regression$bw))))
  expect_true(all(c(density$bw, distribution$bw, regression$bw) > 0))
  expect_true(all(is.finite(fitted(npudens(bws = density, tdat = x)))))
  expect_true(all(is.finite(fitted(npudist(bws = distribution, tdat = x)))))
  expect_true(all(is.finite(fitted(npreg(
    bws = regression, txdat = x, tydat = y
  )))))

  skip_if_not_installed("crs", minimum_version = "0.15.46")
  mads <- do.call(npregbw, c(list(
    xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ls",
    bwsolver = "mads", nomad.opts = list(MAX_BB_EVAL = 12)
  ), common))
  expect_true(is.finite(mads$bw) && mads$bw > 0)
  expect_true(is.finite(mads$fval))
})
