test_that("conditional-moment batches equal scalar kernel contractions", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(29173)
  n <- 32L
  xdat <- data.frame(
    x = rnorm(n),
    group = factor(sample(c("a", "b", "c"), n, TRUE))
  )
  score <- matrix(rnorm(n * 5L), nrow = n)
  bw <- list(bw = c(0.58, 0.24))
  fhat <- seq(0.8, 1.2, length.out = n)
  prodh <- bw[["bw"]][[1L]]

  scalar <- lapply(seq_len(ncol(score)), function(column) {
    value <- score[, column]
    ksum <- npksum(
      txdat = xdat, tydat = value, bws = bw[["bw"]],
      leave.one.out = TRUE, bandwidth.divide = TRUE
    )[["ksum"]]
    ksum2 <- npksum(
      txdat = xdat, tydat = value^2, bws = bw[["bw"]],
      leave.one.out = TRUE, kernel.pow = 2,
      bandwidth.divide = TRUE
    )[["ksum"]]
    In <- sum(value * ksum / fhat) / n^2
    Omega.hat <- 2 * prodh * sum(value^2 * ksum2 / fhat^2) / n^2
    c(In = In, Omega.hat = Omega.hat,
      Jn = n * sqrt(prodh) * In / sqrt(Omega.hat))
  })
  scalar <- do.call(rbind, scalar)

  batch <- getFromNamespace(".np_cms_statistics_batch", "np")(
    xdat = xdat, score = score, bw = bw, fhat = fhat,
    prodh = prodh, pivot = TRUE
  )
  expect_equal(batch[["In"]], scalar[, "In"], tolerance = 2e-14)
  expect_equal(
    batch[["Omega.hat"]], scalar[, "Omega.hat"], tolerance = 2e-14
  )
  expect_equal(batch[["Jn"]], scalar[, "Jn"], tolerance = 2e-14)

  nonpivot <- getFromNamespace(".np_cms_statistics_batch", "np")(
    xdat = xdat, score = score, bw = bw, fhat = fhat,
    prodh = prodh, pivot = FALSE
  )
  expect_identical(names(nonpivot), "In")
  expect_equal(nonpivot[["In"]], scalar[, "In"], tolerance = 2e-14)

  x.nn <- data.frame(x = sort(rnorm(n)))
  bw.nn <- list(bw = 7)
  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    scalar.nn <- vapply(seq_len(ncol(score)), function(column) {
      value <- score[, column]
      ksum <- npksum(
        txdat = x.nn, tydat = value, bws = bw.nn[["bw"]],
        bwtype = bwtype, leave.one.out = TRUE,
        bandwidth.divide = TRUE
      )[["ksum"]]
      sum(value * ksum / fhat) / n^2
    }, numeric(1L))
    batch.nn <- getFromNamespace(".np_cms_statistics_batch", "np")(
      xdat = x.nn, score = score, bw = bw.nn, fhat = fhat,
      prodh = 1, pivot = FALSE,
      kernel.args = list(bwtype = bwtype)
    )
    expect_equal(batch.nn[["In"]], scalar.nn, tolerance = 2e-14)
  }
})

test_that("bootstrap chunk planners enforce bounded linear workspace", {
  cms.chunk <- getFromNamespace(".np_cms_bootstrap_chunk_size", "np")
  deneq.chunk <- getFromNamespace(".npdeneq_count_chunk_size", "np")

  expect_identical(cms.chunk(100L, 399L, FALSE), 16L)
  expect_identical(cms.chunk(100L, 399L, TRUE), 16L)
  expect_identical(cms.chunk(10^7, 399L, TRUE), 1L)
  expect_identical(deneq.chunk(100L, 399L), 16L)
  expect_identical(deneq.chunk(10^7, 399L), 1L)
})

test_that("density-equality count batches preserve literal resampling", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(
    z = c(-1.1, -0.4, 0.2, 0.2, 0.8, 1.4),
    group = factor(c("a", "b", "a", "c", "b", "a"),
                   levels = c("a", "b", "c"))
  )
  y <- data.frame(
    z = c(-0.8, -0.1, 0.5, 0.9, 1.7),
    group = factor(c("c", "a", "b", "b", "c"),
                   levels = c("a", "b", "c"))
  )
  bw.x <- c(0.53, 0.22)
  bw.y <- c(0.61, 0.28)
  pool <- data.frame(rbind(x, y))
  power12 <- getFromNamespace(".npksum_power12", "np")

  statistic <- function(sample.x, sample.y) {
    n1 <- nrow(sample.x)
    n2 <- nrow(sample.y)
    within.x <- power12(
      txdat = sample.x, bws = bw.x, leave.one.out = TRUE,
      bandwidth.divide = TRUE
    )
    within.y <- power12(
      txdat = sample.y, bws = bw.y, leave.one.out = TRUE,
      bandwidth.divide = TRUE
    )
    cross <- power12(
      txdat = sample.x, exdat = sample.y, bws = bw.x,
      bandwidth.divide = TRUE
    )
    In <- sum(within.x[["ksum"]]) / (n1 * (n1 - 1)) +
      sum(within.y[["ksum"]]) / (n2 * (n2 - 1)) -
      2 * sum(cross[["ksum"]]) / (n1 * n2)
    sigma2 <- 2 * (
      sum(within.x[["ksum.power2"]]) / (n1^2 * (n1 - 1)^2) +
      sum(within.y[["ksum.power2"]]) / (n2^2 * (n2 - 1)^2) +
      2 * sum(cross[["ksum.power2"]]) / (n1^2 * n2^2)
    )
    c(Tn = In / sqrt(sigma2), In = In)
  }

  random.seed <- 6329L
  set.seed(random.seed)
  reference <- replicate(9L, {
    sample.x <- pool[sample.int(nrow(pool), nrow(x), TRUE), , drop = FALSE]
    sample.y <- pool[sample.int(nrow(pool), nrow(y), TRUE), , drop = FALSE]
    statistic(sample.x, sample.y)
  })
  observed <- statistic(x, y)

  out <- npdeneqtest(
    x, y, bw.x = bw.x, bw.y = bw.y,
    B = 9L, random.seed = random.seed
  )
  expect_equal(out[["Tn"]], unname(observed[["Tn"]]), tolerance = 2e-13)
  expect_equal(out[["In"]], unname(observed[["In"]]), tolerance = 2e-13)
  expect_equal(out[["Tn.bootstrap"]], reference["Tn", ], tolerance = 2e-13)
  expect_equal(out[["In.bootstrap"]], reference["In", ], tolerance = 2e-13)
  expect_identical(out[["Tn.P"]], mean(reference["Tn", ] > observed[["Tn"]]))
  expect_identical(out[["In.P"]], mean(reference["In", ] > observed[["In"]]))
})
