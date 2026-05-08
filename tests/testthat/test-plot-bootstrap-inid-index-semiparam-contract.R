test_that("npindex inid fast path matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(3231)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)
  bw <- npindexbw(xdat = tx, ydat = y, method = "ichimura", nmulti = 1)
  B <- 11L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  H <- npindexhat(bws = bw, txdat = tx, exdat = tx, output = "matrix", s = 0L)
  fast.out <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  explicit.t <- matrix(NA_real_, nrow = B, ncol = n)
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    explicit.t[b, ] <- npindex(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = tx,
      bws = bw,
      gradients = FALSE
    )$mean
  }

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-10)
  expect_equal(fast.out$t0, as.vector(H %*% y), tolerance = 1e-12)
})

test_that("npindex ll/lp inid fast path matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(32315)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
  rbw.fun <- getFromNamespace(".np_indexhat_rbw", "np")

  cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor"),
    list(regtype = "lp", basis = "glp", label = "lp-glp")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 0.25),
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- 2L
    }
    bw <- do.call(npindexbw, bw.args)

    tx.index <- data.frame(index = as.vector(as.matrix(tx) %*% bw$beta))
    rbw <- rbw.fun(bws = bw, idx.train = tx.index)

    fast.out <- fast.fun(
      xdat = tx.index,
      exdat = tx.index,
      bws = rbw,
      ydat = y,
      B = B,
      counts = counts
    )

    explicit.t <- matrix(NA_real_, nrow = B, ncol = n)
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      explicit.t[b, ] <- npindex(
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = tx,
        bws = bw,
        gradients = FALSE
      )$mean
    }

    expect_equal(fast.out$t, explicit.t, tolerance = 1e-6, info = cfg$label)
    expect_equal(
      fast.out$t0,
      npindex(txdat = tx, tydat = y, exdat = tx, bws = bw, gradients = FALSE)$mean,
      tolerance = 1e-6,
      info = cfg$label
    )
  }
})

test_that("npindex plot bootstrap inid supports ll/lp basis variants", {
  skip_if_not_installed("np")

  set.seed(3232)
  n <- 70
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)

  cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "glp", label = "lp-glp"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 0.25),
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- 2L
    }
    bw <- do.call(npindexbw, bw.args)
    out <- suppressWarnings(
      plot(
        bw,
        xdat = tx,
        ydat = y,
        output = "data",
        perspective = FALSE,
        errors = "bootstrap",
        bootstrap = "inid",
        B = 7
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0, info = cfg$label)
  }
})

test_that("npindex plot bootstrap inid fails fast for unsupported nonfixed gradients", {
  skip_if_not_installed("np")

  set.seed(32321)
  n <- 60
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      bwtype = bt
    )

    expect_error(
      suppressWarnings(
        plot(
          bw,
          xdat = tx,
          ydat = y,
          output = "data",
          perspective = FALSE,
          gradients = TRUE,
          errors = "bootstrap",
          bootstrap = "inid",
          B = 7
        )
      ),
      "gradients=FALSE",
      info = bt
    )
  }
})

test_that("npindex bounded bootstrap plot-data supports bw and fit objects", {
  skip_if_not_installed("np")

  set.seed(32322)
  n <- 36
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)

  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    method = "ichimura",
    ckerbound = "range",
    bws = c(1, 0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npindex(bws = bw, txdat = tx, tydat = y)

  run_plot <- function(obj, boot.method) {
    suppressWarnings(
      plot(
        obj,
        xdat = tx,
        ydat = y,
        output = "data",
        perspective = FALSE,
        errors = "bootstrap",
        bootstrap = boot.method,
        B = 5L,
        band = "pointwise"
      )
    )
  }

  for (boot.method in c("inid", "geom", "wild")) {
    out.bw <- run_plot(bw, boot.method)
    out.fit <- run_plot(fit, boot.method)

    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
    expect_true(length(out.bw) > 0L, info = paste("bw", boot.method))
    expect_true(length(out.fit) > 0L, info = paste("fit", boot.method))
  }
})
