test_that("npksum numeric and formula interfaces agree", {
  set.seed(1)
  n <- 30
  x <- runif(n)
  y <- rnorm(n)
  dat <- data.frame(y = y, x = x)

  k1 <- npksum(txdat = x, tydat = y, bws = 0.5)
  k2 <- npksum(y ~ x, data = dat, bws = 0.5)

  expect_s3_class(k1, "npkernelsum")
  expect_s3_class(k2, "npkernelsum")
  expect_identical(k1$ksum, k2$ksum)
})

test_that("npksum preserves 1D exdat column naming", {
  set.seed(1)
  n <- 10
  x <- runif(n)
  ex <- seq(0, 1, length.out = 3)

  k <- npksum(txdat = x, exdat = ex, bws = 0.5)
  expect_identical(colnames(k$eval), "exdat")
})

test_that("npksum positional bws dispatch matches named bws", {
  set.seed(42)
  n <- 20
  x <- runif(n)
  y <- rnorm(n)

  k_named <- npksum(txdat = x, tydat = y, bws = 0.4)
  k_pos <- npksum(0.4, txdat = x, tydat = y)

  expect_s3_class(k_pos, "npkernelsum")
  expect_identical(k_named$ksum, k_pos$ksum)
})

test_that("npksum formula subset and na.action match explicit data path", {
  set.seed(7)
  n <- 24
  dat <- data.frame(x = runif(n), y = rnorm(n))
  dat$y[c(3, 12)] <- NA_real_

  k_formula <- npksum(y ~ x, data = dat, subset = x > 0.2, na.action = na.omit, bws = 0.45)

  dat2 <- na.omit(subset(dat, x > 0.2))
  k_explicit <- npksum(txdat = dat2$x, tydat = dat2$y, bws = 0.45)

  expect_identical(k_formula$ksum, k_explicit$ksum)
})

test_that("npksum bounded nonfixed helper path works for nearest-neighbor bandwidths", {
  set.seed(20260325)
  n <- 24
  x <- runif(n)
  y <- rnorm(n)
  w <- cbind(y, 1.0)

  k_gnn <- npksum(
    txdat = x,
    tydat = w,
    weights = w,
    leave.one.out = TRUE,
    bandwidth.divide = TRUE,
    bws = 6,
    bwtype = "generalized_nn",
    ckerbound = "range"
  )
  k_adapt <- npksum(
    txdat = x,
    tydat = w,
    weights = w,
    leave.one.out = TRUE,
    bandwidth.divide = TRUE,
    bws = 7,
    bwtype = "adaptive_nn",
    ckerbound = "range"
  )

  expect_s3_class(k_gnn, "npkernelsum")
  expect_s3_class(k_adapt, "npkernelsum")
  expect_true(all(is.finite(as.numeric(k_gnn$ksum))))
  expect_true(all(is.finite(as.numeric(k_adapt$ksum))))
})

test_that("npksum validates public scalar controls before native dispatch", {
  x <- data.frame(x = seq(0.1, 0.9, length.out = 7L))
  invalid <- list(NA, logical(0), c(TRUE, FALSE))

  for (argname in c(
    "bandwidth.divide", "compute.ocg", "compute.score",
    "leave.one.out", "return.kernel.weights"
  )) {
    for (value in invalid) {
      args <- list(txdat = x, exdat = x, bws = 0.25)
      args[[argname]] <- value
      expect_error(
        do.call(npksum, args),
        sprintf("'%s' must be TRUE or FALSE", argname),
        fixed = TRUE
      )
    }
  }
})

test_that("npksum formula, numeric, and default routes share scalar validation", {
  x <- seq(0.1, 0.9, length.out = 7L)
  dat <- data.frame(y = sin(2 * pi * x), x = x)
  txdat <- dat["x"]
  bws <- np:::kbandwidth(
    bw = 0.25,
    xdati = np:::untangle(txdat),
    xnames = names(txdat)
  )

  expect_error(
    npksum(y ~ x, data = dat, bws = 0.25, compute.score = c(TRUE, FALSE)),
    "'compute.score' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    npksum(txdat = txdat, exdat = txdat, bws = 0.25,
           leave.one.out = logical(0)),
    "'leave.one.out' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    np:::npksum.default(bws = bws, txdat = txdat,
                        return.kernel.weights = NA),
    "'return.kernel.weights' must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("npksum preserves established scalar logical coercion", {
  x <- seq(0.1, 0.9, length.out = 9L)
  u <- factor(rep(c("a", "b", "c"), length.out = length(x)))
  xu <- data.frame(x = x, u = u)

  expect_identical(
    npksum(txdat = x, exdat = x, bws = 0.2, bandwidth.divide = 1)$ksum,
    npksum(txdat = x, exdat = x, bws = 0.2, bandwidth.divide = TRUE)$ksum
  )
  expect_identical(
    npksum(txdat = x, bws = 0.2, leave.one.out = 1)$ksum,
    npksum(txdat = x, bws = 0.2, leave.one.out = TRUE)$ksum
  )
  expect_identical(
    npksum(txdat = x, exdat = x, bws = 0.2, return.kernel.weights = 1)$kw,
    npksum(txdat = x, exdat = x, bws = 0.2, return.kernel.weights = TRUE)$kw
  )
  expect_identical(
    npksum(txdat = xu, exdat = xu, bws = c(0.2, 0.25), compute.score = 1)$p.ksum,
    npksum(txdat = xu, exdat = xu, bws = c(0.2, 0.25), compute.score = TRUE)$p.ksum
  )
  expect_identical(
    npksum(txdat = xu, exdat = xu, bws = c(0.2, 0.25), compute.ocg = 1)$p.ksum,
    npksum(txdat = xu, exdat = xu, bws = c(0.2, 0.25), compute.ocg = TRUE)$p.ksum
  )
})

test_that("npksum rejects unsafe kernel powers without narrowing valid powers", {
  x <- data.frame(x = seq(0.1, 0.9, length.out = 7L))

  for (kernel.pow in c(-1L, 0L, 1L, 2L)) {
    result <- npksum(txdat = x, exdat = x, bws = 0.25, kernel.pow = kernel.pow)
    expect_s3_class(result, "npkernelsum")
    expect_true(all(is.finite(result$ksum)))
  }

  invalid <- list(
    NA_integer_, NaN, Inf, -Inf, integer(0), c(1L, 2L), 1.5,
    .Machine$integer.max + 1
  )
  for (kernel.pow in invalid) {
    expect_error(
      npksum(txdat = x, exdat = x, bws = 0.25, kernel.pow = kernel.pow),
      "'kernel.pow' must be one finite integer",
      fixed = TRUE
    )
  }
})

test_that("npksum mixed-data score and OCG blocks are correct and ordered", {
  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.08, 0.92, length.out = 18L)
  u <- factor(rep(c("a", "b", "c"), length.out = length(x)))
  xu <- data.frame(u = u, x = x)
  bw <- c(0.24, 0.19)

  score <- npksum(txdat = xu, exdat = xu, bws = bw, compute.score = TRUE)
  eps <- 1e-6
  plus <- npksum(txdat = xu, exdat = xu, bws = c(bw[1L] + eps, bw[2L]))
  minus <- npksum(txdat = xu, exdat = xu, bws = c(bw[1L] - eps, bw[2L]))
  score.fd <- (plus$ksum - minus$ksum) / (2 * eps)
  expect_equal(as.numeric(score$p.ksum), as.numeric(score.fd), tolerance = 2e-7)

  ocg <- npksum(
    txdat = xu, exdat = xu, bws = bw, compute.ocg = TRUE,
    return.kernel.weights = TRUE
  )
  xu.reference <- xu
  xu.reference$u <- factor(
    rep(levels(xu$u)[1L], nrow(xu)),
    levels = levels(xu$u)
  )
  ocg.reference <- npksum(
    txdat = xu, exdat = xu.reference, bws = bw,
    return.kernel.weights = TRUE
  )
  expect_equal(as.numeric(ocg$p.ksum), as.numeric(ocg.reference$ksum), tolerance = 1e-13)

  derivative <- npksum(
    txdat = xu,
    exdat = xu,
    bws = bw,
    permutation.operator = "derivative",
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  combined <- npksum(
    txdat = xu,
    exdat = xu,
    bws = bw,
    permutation.operator = "derivative",
    compute.ocg = TRUE,
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  xu.plus <- xu
  xu.minus <- xu
  xu.plus$x <- xu.plus$x + eps
  xu.minus$x <- xu.minus$x - eps
  derivative.fd <- (
    npksum(txdat = xu, exdat = xu.plus, bws = bw)$ksum -
      npksum(txdat = xu, exdat = xu.minus, bws = bw)$ksum
  ) / (2 * eps)
  expect_equal(as.numeric(derivative$p.ksum), as.numeric(derivative.fd), tolerance = 2e-7)
  expect_identical(dim(combined$p.ksum), c(nrow(xu), 2L))
  expect_equal(combined$p.ksum[, 1L], as.numeric(ocg$p.ksum), tolerance = 1e-13)
  expect_equal(combined$p.ksum[, 2L], as.numeric(derivative$p.ksum), tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 1L], ocg.reference$kw, tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 2L], derivative$p.kw, tolerance = 1e-13)

  options(np.tree = TRUE)
  score.tree <- npksum(txdat = xu, exdat = xu, bws = bw, compute.score = TRUE)
  ocg.tree <- npksum(txdat = xu, exdat = xu, bws = bw, compute.ocg = TRUE)
  derivative.tree <- npksum(
    txdat = xu,
    exdat = xu,
    bws = bw,
    permutation.operator = "derivative"
  )
  combined.tree <- npksum(
    txdat = xu,
    exdat = xu,
    bws = bw,
    permutation.operator = "derivative",
    compute.ocg = TRUE
  )
  expect_equal(score.tree$p.ksum, score$p.ksum, tolerance = 1e-13)
  expect_equal(ocg.tree$p.ksum, ocg$p.ksum, tolerance = 1e-13)
  expect_equal(derivative.tree$p.ksum, derivative$p.ksum, tolerance = 1e-13)
  expect_equal(combined.tree$p.ksum, combined$p.ksum, tolerance = 1e-13)

  expect_error(
    npksum(
      txdat = xu,
      exdat = xu,
      bws = bw,
      permutation.operator = "derivative",
      compute.score = TRUE
    ),
    "compute.score cannot be combined with a permutation operator",
    fixed = TRUE
  )
})

test_that("npksum ordered score, OCG, and mixed derivatives match direct oracles", {
  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.08, 0.92, length.out = 18L)
  o <- ordered(rep(c("low", "middle", "high"), length.out = length(x)),
               levels = c("low", "middle", "high"))
  xo <- data.frame(x = x, o = o)
  bw <- c(0.19, 0.31)
  eps <- 1e-6

  score <- npksum(txdat = xo, exdat = xo, bws = bw, compute.score = TRUE)
  score.fd <- (
    npksum(txdat = xo, exdat = xo, bws = c(bw[1L], bw[2L] + eps))$ksum -
      npksum(txdat = xo, exdat = xo, bws = c(bw[1L], bw[2L] - eps))$ksum
  ) / (2 * eps)
  expect_equal(as.numeric(score$p.ksum), as.numeric(score.fd), tolerance = 2e-7)

  ocg <- npksum(
    txdat = xo, exdat = xo, bws = bw, compute.ocg = TRUE,
    return.kernel.weights = TRUE
  )
  reference.index <- ifelse(as.integer(xo$o) == 1L, 2L, as.integer(xo$o) - 1L)
  xo.reference <- xo
  xo.reference$o <- ordered(
    levels(xo$o)[reference.index],
    levels = levels(xo$o)
  )
  ocg.reference <- npksum(
    txdat = xo, exdat = xo.reference, bws = bw,
    return.kernel.weights = TRUE
  )
  expect_equal(as.numeric(ocg$p.ksum), as.numeric(ocg.reference$ksum), tolerance = 1e-13)

  derivative <- npksum(
    txdat = xo,
    exdat = xo,
    bws = bw,
    permutation.operator = "derivative",
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  xo.plus <- xo
  xo.minus <- xo
  xo.plus$x <- xo.plus$x + eps
  xo.minus$x <- xo.minus$x - eps
  derivative.fd <- (
    npksum(txdat = xo, exdat = xo.plus, bws = bw)$ksum -
      npksum(txdat = xo, exdat = xo.minus, bws = bw)$ksum
  ) / (2 * eps)
  expect_equal(as.numeric(derivative$p.ksum), as.numeric(derivative.fd), tolerance = 2e-7)

  combined <- npksum(
    txdat = xo,
    exdat = xo,
    bws = bw,
    permutation.operator = "derivative",
    compute.ocg = TRUE,
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  expect_equal(combined$p.ksum[, 1L], as.numeric(derivative$p.ksum), tolerance = 1e-13)
  expect_equal(combined$p.ksum[, 2L], as.numeric(ocg$p.ksum), tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 1L], derivative$p.kw, tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 2L], ocg.reference$kw, tolerance = 1e-13)

  options(np.tree = TRUE)
  score.tree <- npksum(txdat = xo, exdat = xo, bws = bw, compute.score = TRUE)
  ocg.tree <- npksum(txdat = xo, exdat = xo, bws = bw, compute.ocg = TRUE)
  derivative.tree <- npksum(
    txdat = xo,
    exdat = xo,
    bws = bw,
    permutation.operator = "derivative"
  )
  combined.tree <- npksum(
    txdat = xo,
    exdat = xo,
    bws = bw,
    permutation.operator = "derivative",
    compute.ocg = TRUE
  )
  expect_equal(score.tree$p.ksum, score$p.ksum, tolerance = 1e-13)
  expect_equal(ocg.tree$p.ksum, ocg$p.ksum, tolerance = 1e-13)
  expect_equal(derivative.tree$p.ksum, derivative$p.ksum, tolerance = 1e-13)
  expect_equal(combined.tree$p.ksum, combined$p.ksum, tolerance = 1e-13)
})

test_that("npksum mixed score and permutation buffers support NN bandwidths", {
  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260721)
  x <- sort(runif(30L, 0.05, 0.95))
  u <- factor(rep(c("a", "b", "c"), length.out = length(x)))
  xu <- data.frame(u = u, x = x)

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    common <- list(txdat = xu, exdat = xu, bws = c(0.24, 7), bwtype = bwtype)
    reference <- list(
      score = do.call(npksum, c(common, list(compute.score = TRUE))),
      ocg = do.call(npksum, c(common, list(compute.ocg = TRUE))),
      derivative = do.call(npksum, c(
        common,
        list(permutation.operator = "derivative")
      )),
      combined = do.call(npksum, c(
        common,
        list(
          permutation.operator = "derivative",
          compute.ocg = TRUE,
          return.kernel.weights = TRUE,
          return.derivative.kernel.weights = TRUE
        )
      ))
    )

    for (result in reference)
      expect_true(all(is.finite(result$p.ksum)))
    expect_identical(dim(reference$combined$p.ksum), c(nrow(xu), 2L))
    expect_identical(dim(reference$combined$p.kw), c(nrow(xu), nrow(xu), 2L))
    expect_true(all(is.finite(reference$combined$p.kw)))

    options(np.tree = TRUE)
    tree <- list(
      score = do.call(npksum, c(common, list(compute.score = TRUE))),
      ocg = do.call(npksum, c(common, list(compute.ocg = TRUE))),
      derivative = do.call(npksum, c(
        common,
        list(permutation.operator = "derivative")
      )),
      combined = do.call(npksum, c(
        common,
        list(
          permutation.operator = "derivative",
          compute.ocg = TRUE,
          return.kernel.weights = TRUE,
          return.derivative.kernel.weights = TRUE
        )
      ))
    )
    options(np.tree = FALSE)

    for (name in names(reference))
      expect_equal(tree[[name]]$p.ksum, reference[[name]]$p.ksum,
                   tolerance = 1e-13)
    expect_equal(tree$combined$p.kw, reference$combined$p.kw,
                 tolerance = 1e-13)
  }
})
