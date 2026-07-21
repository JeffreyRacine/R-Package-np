test_that("npksum validates scalar controls before MPI/native dispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  x <- data.frame(x = seq(0.1, 0.9, length.out = 7L))
  invalid.logical <- list(NA, logical(0), c(TRUE, FALSE))

  for (argname in c(
    "bandwidth.divide", "compute.ocg", "compute.score",
    "leave.one.out", "return.kernel.weights"
  )) {
    for (value in invalid.logical) {
      args <- list(txdat = x, exdat = x, bws = 0.25)
      args[[argname]] <- value
      expect_error(
        do.call(npksum, args),
        sprintf("'%s' must be TRUE or FALSE", argname),
        fixed = TRUE
      )
    }
  }

  invalid.power <- list(
    NA_integer_, NaN, Inf, -Inf, integer(0), c(1L, 2L), 1.5,
    .Machine$integer.max + 1
  )
  for (kernel.pow in invalid.power) {
    expect_error(
      npksum(txdat = x, exdat = x, bws = 0.25, kernel.pow = kernel.pow),
      "'kernel.pow' must be one finite integer",
      fixed = TRUE
    )
  }

  for (kernel.pow in c(-1L, 0L, 1L, 2L)) {
    result <- npksum(txdat = x, exdat = x, bws = 0.25,
                     kernel.pow = kernel.pow)
    expect_s3_class(result, "npkernelsum")
    expect_true(all(is.finite(result$ksum)))
  }

  dat <- data.frame(y = sin(2 * pi * x$x), x = x$x)
  bws <- npRmpi:::kbandwidth(
    bw = 0.25,
    xdati = npRmpi:::untangle(x),
    xnames = names(x)
  )
  expect_error(
    npksum(y ~ x, data = dat, bws = 0.25,
           compute.score = c(TRUE, FALSE)),
    "'compute.score' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::npksum.default(bws = bws, txdat = x,
                            return.kernel.weights = NA),
    "'return.kernel.weights' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_identical(
    npksum(txdat = x, exdat = x, bws = 0.25,
           bandwidth.divide = 1)$ksum,
    npksum(txdat = x, exdat = x, bws = 0.25,
           bandwidth.divide = TRUE)$ksum
  )

  u <- factor(rep(c("a", "b", "c"), length.out = nrow(x)))
  xu <- data.frame(u = u, x = x$x)
  expect_error(
    npksum(
      txdat = xu, exdat = xu, bws = c(0.24, 0.19),
      permutation.operator = "derivative", compute.score = TRUE
    ),
    "compute.score cannot be combined with a permutation operator",
    fixed = TRUE
  )
})

test_that("npksum unordered mixed score and permutation blocks match oracles", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.08, 0.92, length.out = 19L)
  u <- factor(rep(c("a", "b", "c"), length.out = length(x)))
  xu <- data.frame(u = u, x = x)
  bw <- c(0.24, 0.19)
  eps <- 1e-6

  score <- npksum(txdat = xu, exdat = xu, bws = bw, compute.score = TRUE)
  score.fd <- (
    npksum(txdat = xu, exdat = xu,
           bws = c(bw[1L] + eps, bw[2L]))$ksum -
      npksum(txdat = xu, exdat = xu,
             bws = c(bw[1L] - eps, bw[2L]))$ksum
  ) / (2 * eps)
  expect_equal(as.numeric(score$p.ksum), as.numeric(score.fd),
               tolerance = 2e-7)

  ocg <- npksum(
    txdat = xu, exdat = xu, bws = bw, compute.ocg = TRUE,
    return.kernel.weights = TRUE
  )
  xu.reference <- xu
  xu.reference$u <- factor(
    rep(levels(xu$u)[1L], nrow(xu)), levels = levels(xu$u)
  )
  ocg.reference <- npksum(
    txdat = xu, exdat = xu.reference, bws = bw,
    return.kernel.weights = TRUE
  )
  expect_equal(as.numeric(ocg$p.ksum), as.numeric(ocg.reference$ksum),
               tolerance = 1e-13)

  derivative <- npksum(
    txdat = xu, exdat = xu, bws = bw,
    permutation.operator = "derivative",
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
  expect_equal(as.numeric(derivative$p.ksum), as.numeric(derivative.fd),
               tolerance = 2e-7)

  combined <- npksum(
    txdat = xu, exdat = xu, bws = bw,
    permutation.operator = "derivative", compute.ocg = TRUE,
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  expect_identical(dim(combined$p.ksum), c(nrow(xu), 2L))
  expect_equal(combined$p.ksum[, 1L], as.numeric(ocg$p.ksum),
               tolerance = 1e-13)
  expect_equal(combined$p.ksum[, 2L], as.numeric(derivative$p.ksum),
               tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 1L], ocg.reference$kw,
               tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 2L], derivative$p.kw,
               tolerance = 1e-13)

  options(np.tree = TRUE)
  expect_equal(
    npksum(txdat = xu, exdat = xu, bws = bw,
           compute.score = TRUE)$p.ksum,
    score$p.ksum,
    tolerance = 1e-13
  )
  expect_equal(
    npksum(txdat = xu, exdat = xu, bws = bw,
           permutation.operator = "derivative",
           compute.ocg = TRUE)$p.ksum,
    combined$p.ksum,
    tolerance = 1e-13
  )
})

test_that("npksum ordered score and OCG paths match direct oracles", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.08, 0.92, length.out = 19L)
  o <- ordered(rep(c("low", "middle", "high"), length.out = length(x)),
               levels = c("low", "middle", "high"))
  xo <- data.frame(x = x, o = o)
  bw <- c(0.19, 0.31)
  eps <- 1e-6

  score <- npksum(txdat = xo, exdat = xo, bws = bw, compute.score = TRUE)
  score.fd <- (
    npksum(txdat = xo, exdat = xo,
           bws = c(bw[1L], bw[2L] + eps))$ksum -
      npksum(txdat = xo, exdat = xo,
             bws = c(bw[1L], bw[2L] - eps))$ksum
  ) / (2 * eps)
  expect_equal(as.numeric(score$p.ksum), as.numeric(score.fd),
               tolerance = 2e-7)

  ocg <- npksum(
    txdat = xo, exdat = xo, bws = bw, compute.ocg = TRUE,
    return.kernel.weights = TRUE
  )
  reference.index <- ifelse(as.integer(xo$o) == 1L, 2L,
                            as.integer(xo$o) - 1L)
  xo.reference <- xo
  xo.reference$o <- ordered(levels(xo$o)[reference.index],
                            levels = levels(xo$o))
  ocg.reference <- npksum(
    txdat = xo, exdat = xo.reference, bws = bw,
    return.kernel.weights = TRUE
  )
  expect_equal(as.numeric(ocg$p.ksum), as.numeric(ocg.reference$ksum),
               tolerance = 1e-13)

  derivative <- npksum(
    txdat = xo, exdat = xo, bws = bw,
    permutation.operator = "derivative",
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  combined <- npksum(
    txdat = xo, exdat = xo, bws = bw,
    permutation.operator = "derivative", compute.ocg = TRUE,
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  expect_equal(combined$p.ksum[, 1L], as.numeric(derivative$p.ksum),
               tolerance = 1e-13)
  expect_equal(combined$p.ksum[, 2L], as.numeric(ocg$p.ksum),
               tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 1L], derivative$p.kw,
               tolerance = 1e-13)
  expect_equal(combined$p.kw[, , 2L], ocg.reference$kw,
               tolerance = 1e-13)
})

test_that("npksum mixed score and permutation buffers support NN bandwidths", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
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
        list(
          permutation.operator = "derivative",
          return.kernel.weights = TRUE,
          return.derivative.kernel.weights = TRUE
        )
      )),
      combined = do.call(npksum, c(
        common,
        list(
          permutation.operator = "derivative", compute.ocg = TRUE,
          return.kernel.weights = TRUE,
          return.derivative.kernel.weights = TRUE
        )
      ))
    )
    for (result in reference)
      expect_true(all(is.finite(result$p.ksum)))
    expect_identical(dim(reference$combined$p.kw),
                     c(nrow(xu), nrow(xu), 2L))
    expect_true(all(is.finite(reference$combined$p.kw)))

    xu.ocg.reference <- xu
    xu.ocg.reference$u <- factor(
      rep(levels(xu$u)[1L], nrow(xu)), levels = levels(xu$u)
    )
    ocg.common <- common
    ocg.common$exdat <- xu.ocg.reference
    ocg.reference <- do.call(
      npksum, c(ocg.common, list(return.kernel.weights = TRUE))
    )
    expect_equal(reference$combined$p.kw[, , 1L], ocg.reference$kw,
                 tolerance = 1e-13)
    expect_equal(reference$combined$p.kw[, , 2L],
                 reference$derivative$p.kw,
                 tolerance = 1e-13)

    options(np.tree = TRUE)
    tree <- list(
      score = do.call(npksum, c(common, list(compute.score = TRUE))),
      ocg = do.call(npksum, c(common, list(compute.ocg = TRUE))),
      derivative = do.call(npksum, c(
        common, list(permutation.operator = "derivative")
      )),
      combined = do.call(npksum, c(
        common,
        list(
          permutation.operator = "derivative", compute.ocg = TRUE,
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
