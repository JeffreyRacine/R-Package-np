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

test_that("npindexbw BFGS accounting uses bandwidth callback count", {
  target_is_num_feval <- function(expr) {
    is.call(expr) &&
      identical(expr[[1L]], quote(`$`)) &&
      identical(expr[[2L]], quote(bws)) &&
      identical(expr[[3L]], quote(num.feval))
  }

  assigns_from_counter <- function(expr) {
    ok <- FALSE
    walk <- function(node) {
      if (is.call(node)) {
        if (identical(node[[1L]], quote(`<-`)) &&
            length(node) == 3L &&
            target_is_num_feval(node[[2L]]) &&
            identical(node[[3L]], quote(bandwidth_eval_count))) {
          ok <<- TRUE
        }
        for (idx in seq_along(node)[-1L])
          walk(node[[idx]])
      } else if (is.pairlist(node) || is.expression(node)) {
        for (idx in seq_along(node))
          walk(node[[idx]])
      }
    }
    walk(expr)
    ok
  }

  increments_counter <- function(expr) {
    ok <- FALSE
    walk <- function(node) {
      if (is.call(node)) {
        if (identical(node[[1L]], quote(`<<-`)) &&
            length(node) == 3L &&
            identical(node[[2L]], quote(bandwidth_eval_count)) &&
            identical(node[[3L]], quote(bandwidth_eval_count + 1L))) {
          ok <<- TRUE
        }
        for (idx in seq_along(node)[-1L])
          walk(node[[idx]])
      } else if (is.pairlist(node) || is.expression(node)) {
        for (idx in seq_along(node))
          walk(node[[idx]])
      }
    }
    walk(expr)
    ok
  }

  fn <- getFromNamespace("npindexbw.sibandwidth", "np")
  expect_true(increments_counter(body(fn)))
  expect_true(assigns_from_counter(body(fn)))
})

test_that("npindexbw BFGS adaptive-nn accounting exceeds fast callbacks", {
  set.seed(20260627)
  n <- 45L
  x1 <- runif(n)
  x2 <- runif(n)
  x3 <- runif(n)
  y <- sin(x1 + 0.5 * x2 - 0.25 * x3) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  bw <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    method = "ichimura",
    regtype = "lc",
    bwtype = "adaptive_nn",
    nmulti = 1L,
    optim.method = "BFGS",
    optim.maxit = 80L,
    optim.maxattempts = 1L
  ))

  num.feval <- as.numeric(bw$num.feval[1L])
  num.feval.fast <- as.numeric(bw$num.feval.fast[1L])

  expect_true(is.finite(as.numeric(bw$fval[1L])))
  expect_gt(num.feval.fast, 0)
  expect_gt(num.feval, num.feval.fast)
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

test_that("npindex bootstrap error SD recovery matches covariance diagonal reduction", {
  boot.sds <- getFromNamespace(".np_plot_bootstrap_col_sds", "np")

  set.seed(20260311)
  boot.t <- matrix(rnorm(250 * 24), nrow = 250, ncol = 24)

  ref.mean <- sqrt(diag(cov(boot.t[, 1:12, drop = FALSE])))
  ref.grad <- sqrt(diag(cov(boot.t[, 13:24, drop = FALSE])))

  expect_equal(
    boot.sds(boot.t[, 1:12, drop = FALSE]),
    ref.mean,
    tolerance = 1e-15
  )
  expect_equal(
    boot.sds(boot.t[, 13:24, drop = FALSE]),
    ref.grad,
    tolerance = 1e-15
  )
})

test_that("predict.singleindex forwards explicit boot.num without changing fitted predictions", {
  set.seed(20260323)
  n <- 60L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.05)
  tx <- data.frame(x1 = x1, x2 = x2)
  nd <- tx[seq_len(8L), , drop = FALSE]

  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 0.35, 0.45),
    bandwidth.compute = FALSE,
    method = "ichimura"
  )
  fit <- npindex(bws = bw, txdat = tx, tydat = y)

  pred <- predict(fit, newdata = nd)
  pred.se <- predict(fit, newdata = nd, se.fit = TRUE, boot.num = 10)

  expect_equal(as.numeric(pred.se$fit), as.numeric(pred), tolerance = 0)
  expect_equal(length(pred.se$se.fit), nrow(nd))
})
