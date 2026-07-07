test_that("single-index fixed wild bootstrap training mean preserves hat apply", {
  set.seed(20260702)
  n <- 60L
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * tx$x1 - tx$x2) + rnorm(n, sd = 0.05)

  old.apply <- getFromNamespace(".np_plot_singleindex_hat_apply_index", "np")
  train.mean <- getFromNamespace(".np_plot_singleindex_training_mean_index", "np")

  cfgs <- list(
    list(regtype = "lc", label = "lc"),
    list(regtype = "ll", label = "ll"),
    list(regtype = "lp", basis = "tensor", degree = 2L, label = "lp")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 0.35),
      bandwidth.compute = FALSE,
      regtype = cfg$regtype,
      bwtype = "fixed"
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- cfg$degree
    }
    bw <- do.call(npindexbw, bw.args)
    idx <- data.frame(index = as.vector(as.matrix(tx) %*% bw$beta))

    expect_equal(
      train.mean(bws = bw, idx.train = idx, y = y),
      as.vector(old.apply(bws = bw, idx.train = idx, idx.eval = idx, y = y)),
      tolerance = 1e-10,
      info = cfg$label
    )
  }
})

test_that("single-index fixed wild bootstrap block path matches dense path", {
  set.seed(20260703)
  n <- 95L
  tx <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  y <- sin(tx$x1 + tx$x2 - tx$x3) + rnorm(n, sd = 0.06)
  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 1, 1, 0.45),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed"
  )
  fit <- npindex(bws = bw, txdat = tx, tydat = y, gradients = TRUE)

  old.threshold <- getOption("np.plot.wild.apply.operator.threshold.bytes")
  on.exit(options(np.plot.wild.apply.operator.threshold.bytes = old.threshold),
          add = TRUE)

  options(np.plot.wild.apply.operator.threshold.bytes = Inf)
  set.seed(20260704)
  dense <- suppressWarnings(plot(
    fit,
    errors = "bootstrap",
    bootstrap = "wild",
    B = 5L,
    neval = 23L,
    output = "data",
    perspective = FALSE,
    data_overlay = FALSE
  ))[[1L]]

  options(np.plot.wild.apply.operator.threshold.bytes = 0)
  set.seed(20260704)
  block <- suppressWarnings(plot(
    fit,
    errors = "bootstrap",
    bootstrap = "wild",
    B = 5L,
    neval = 23L,
    output = "data",
    perspective = FALSE,
    data_overlay = FALSE
  ))[[1L]]

  expect_equal(block$mean, dense$mean, tolerance = 1e-12)
  expect_equal(block$merr, dense$merr, tolerance = 1e-12)
})

test_that("single-index range-bounded wild bootstrap uses scalar-index bounds", {
  set.seed(32322)
  n <- 36L
  tx <- data.frame(x1 = runif(n), x2 = 2 * runif(n) - 1)
  y <- sin(tx$x1 + 0.4 * tx$x2) + rnorm(n, sd = 0.05)

  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    method = "ichimura",
    ckerbound = "range",
    bws = c(1, 0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npindex(bws = bw, txdat = tx, tydat = y)

  for (object in list(bw, fit)) {
    out <- suppressWarnings(plot(
      object,
      xdat = tx,
      ydat = y,
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      bootstrap = "wild",
      B = 5L,
      band = "pointwise"
    ))

    expect_type(out, "list")
    expect_true(length(out) >= 1L)
    expect_true(all(is.finite(out[[1L]]$mean)))
    expect_true(all(is.finite(out[[1L]]$merr)))
  }
})
