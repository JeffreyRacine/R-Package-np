gnn_literal_radius <- function(train, evaluation, k, exclude = integer()) {
  keep <- if (length(exclude)) setdiff(seq_along(train), exclude) else seq_along(train)
  distance <- sort(abs(evaluation - train[keep]), method = "radix")
  lookup <- min(k, length(distance))
  distance[[lookup]] * k / lookup
}

gnn_literal_density <- function(train, evaluation, k, mapped = FALSE,
                                leave_one_out = FALSE) {
  vapply(seq_along(evaluation), function(i) {
    exclude <- if (mapped || leave_one_out) i else integer()
    h <- gnn_literal_radius(train, evaluation[[i]], k, exclude)
    keep <- if (leave_one_out) setdiff(seq_along(train), i) else seq_along(train)
    mean(dnorm((evaluation[[i]] - train[keep]) / h) / h)
  }, numeric(1))
}

gnn_literal_distribution <- function(train, evaluation, k, mapped = FALSE) {
  vapply(seq_along(evaluation), function(i) {
    exclude <- if (mapped) i else integer()
    h <- gnn_literal_radius(train, evaluation[[i]], k, exclude)
    mean(pnorm((evaluation[[i]] - train) / h))
  }, numeric(1))
}

test_that("ordinary generalized-NN unconditional owners delete mapped occurrences", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- c(-1.4, -0.6, -0.1, 0.2, 0.9, 1.7)
  dat <- data.frame(x = x)
  density_bw <- npudensbw(
    dat = dat, bwtype = "generalized_nn", bwmethod = "cv.ml",
    bws = 2L, bandwidth.compute = FALSE
  )
  distribution_bw <- npudistbw(
    dat = dat, bwtype = "generalized_nn", bwmethod = "cv.cdf",
    bws = 2L, bandwidth.compute = FALSE
  )

  expected_density <- gnn_literal_density(x, x, 2L, mapped = TRUE)
  expected_distribution <- gnn_literal_distribution(x, x, 2L, mapped = TRUE)

  for (tree in c(FALSE, TRUE)) {
    options(np.tree = tree)
    expect_equal(
      npudens(bws = density_bw, tdat = dat)$dens,
      expected_density,
      tolerance = 2e-10
    )
    expect_equal(
      npudist(bws = distribution_bw, tdat = dat)$dist,
      expected_distribution,
      tolerance = 2e-10
    )
  }

  options(np.tree = FALSE)
  external <- dat[c(2L, 4L), , drop = FALSE]
  expect_equal(
    npudens(bws = density_bw, tdat = dat, edat = external)$dens,
    gnn_literal_density(x, external$x, 2L),
    tolerance = 2e-10
  )

  eval_only <- np:::npudensbw.bandwidth(
    dat = dat, bws = density_bw, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L
  )
  expected_cvml <- sum(log(gnn_literal_density(
    x, x, 2L, mapped = TRUE, leave_one_out = TRUE
  )))
  expect_equal(eval_only$fval[[1L]], expected_cvml, tolerance = 2e-10)

  options(np.extendednn = TRUE)
  extended_k <- length(x) + 2L
  extended_bw <- npudensbw(
    dat = dat, bwtype = "generalized_nn", bwmethod = "cv.ml",
    bws = extended_k, bandwidth.compute = FALSE
  )
  expected_extended <- gnn_literal_density(
    x, x, extended_k, mapped = TRUE)
  for (tree in c(FALSE, TRUE)) {
    options(np.tree = tree)
    expect_equal(
      npudens(bws = extended_bw, tdat = dat)$dens,
      expected_extended,
      tolerance = 2e-10
    )
  }
})

test_that("ordinary generalized-NN unconditional fit reports zero literal radii", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  dat <- data.frame(x = c(0, 0, 0, 1, 2))
  bw <- npudensbw(
    dat = dat, bwtype = "generalized_nn", bwmethod = "cv.ml",
    bws = 2L, bandwidth.compute = FALSE
  )
  expect_error(
    npudens(bws = bw, tdat = dat),
    "zero literal radius",
    fixed = TRUE
  )
  dist_bw <- npudistbw(
    dat = dat, bwtype = "generalized_nn", bwmethod = "cv.cdf",
    bws = 2L, bandwidth.compute = FALSE
  )
  expect_error(
    npudist(bws = dist_bw, tdat = dat),
    "zero literal radius",
    fixed = TRUE
  )

  invalid <- np:::npudensbw.bandwidth(
    dat = dat, bws = bw, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L
  )
  expect_gt(sum(invalid$invalid.history), 0)

  density_ls_bw <- npudensbw(
    dat = dat, bwtype = "generalized_nn", bwmethod = "cv.ls",
    bws = 2L, bandwidth.compute = FALSE
  )
  invalid_density_ls <- np:::npudensbw.bandwidth(
    dat = dat, bws = density_ls_bw, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L
  )
  invalid_distribution <- np:::npudistbw.dbandwidth(
    dat = dat, bws = dist_bw, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L
  )
  expect_gt(sum(invalid_density_ls$invalid.history), 0)
  expect_gt(sum(invalid_distribution$invalid.history), 0)
})

test_that("mixed generalized-NN distribution search recovers from invalid candidates", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  data("Italy")
  dat <- transform(Italy[seq_len(20L), , drop = FALSE],
                   year = ordered(year))
  bw <- npudistbw(
    ~ year + gdp, data = dat, bwtype = "generalized_nn", nmulti = 1L
  )

  expect_true(all(is.finite(bw$bw)))
  expect_true(all(bw$bw > 0))
  expect_gt(sum(bw$invalid.history), 0)
})
