test_that("npindex basic functionality works", {
  set.seed(42)
  n <- 100
  x1 <- runif(n)
  x2 <- runif(n)
  # Single index model: y = g(x1 + x2) + e
  y <- (x1 + x2)^2 + rnorm(n, sd=0.1)
  
  mydat <- data.frame(y, x1, x2)
  bw <- npindexbw(
    y ~ x1 + x2,
    data = mydat,
    bws = c(1, 0.35, 0.45),
    bandwidth.compute = FALSE,
    method = "ichimura"
  )
  
  model <- npindex(bws=bw)
  
  expect_s3_class(model, "singleindex")
  expect_type(predict(model), "double")
  expect_output(summary(model))
})

test_that("npindexbw nearest-neighbor selection stores integer support and exact fits stay public-green", {
  set.seed(314163)
  n <- 70L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- x1 - x2 + rnorm(n, sd = 0.2)
  tx <- data.frame(x1 = x1, x2 = x2)

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw <- npindexbw(xdat = tx, ydat = y, regtype = "lc", bwtype = bt, nmulti = 1)
    fit <- npindex(bws = bw, txdat = tx, tydat = y, gradients = FALSE)

    expect_true(bw$bw >= 1, info = bt)
    expect_equal(bw$bw, as.double(as.integer(bw$bw)), tolerance = 0, info = bt)
    expect_true(all(is.finite(fit$mean)), info = bt)
  }
})

test_that("manual single-index nearest-neighbor bandwidths fail fast when not integer support", {
  tx <- data.frame(x1 = seq(0.1, 0.9, length.out = 8L), x2 = seq(0.9, 0.1, length.out = 8L))
  y <- tx$x1 - tx$x2

  expect_error(
    npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 0.5),
      bandwidth.compute = FALSE,
      bwtype = "adaptive_nn"
    ),
    "nearest-neighbor bandwidth must be an integer"
  )
})

test_that("npindex supports residual and error branches with evaluation y data", {
  set.seed(20260309)
  n <- 56L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.05)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(18L), , drop = FALSE]
  ey <- sin(ex$x1 + ex$x2) + rnorm(nrow(ex), sd = 0.05)

  for (cfg in list(
    list(regtype = "lc"),
    list(regtype = "ll"),
    list(regtype = "lp", basis = "glp", degree = 1L)
  )) {
    bw.args <- list(
      xdat = tx,
      ydat = y,
      regtype = cfg$regtype,
      bwtype = "adaptive_nn",
      nmulti = 1L
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- cfg$degree
      bw.args$bernstein.basis <- FALSE
    }

    bw <- do.call(npindexbw, bw.args)

    fit.is <- npindex(bws = bw, txdat = tx, tydat = y, gradients = FALSE, residuals = TRUE)
    expect_equal(length(fit.is$resid), nrow(tx), info = cfg$regtype)
    expect_equal(as.vector(fit.is$resid), y - as.vector(fit.is$mean), tolerance = 1e-8, info = cfg$regtype)

    fit.oos <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, eydat = ey, gradients = FALSE, errors = TRUE)
    expect_equal(length(fit.oos$mean), nrow(ex), info = cfg$regtype)
    expect_equal(length(fit.oos$merr), nrow(ex), info = cfg$regtype)
    expect_true(all(is.finite(fit.oos$mean)), info = cfg$regtype)
    expect_true(all(is.finite(fit.oos$merr)), info = cfg$regtype)

    fit.grad <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, eydat = ey, gradients = TRUE, errors = TRUE)
    expect_equal(dim(fit.grad$grad), c(nrow(ex), ncol(tx)), info = cfg$regtype)
    expect_equal(dim(fit.grad$gerr), c(nrow(ex), ncol(tx)), info = cfg$regtype)
    expect_true(all(is.finite(fit.grad$grad)), info = cfg$regtype)
    expect_true(all(is.finite(fit.grad$gerr)), info = cfg$regtype)
  }
})
