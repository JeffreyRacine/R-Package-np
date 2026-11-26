test_that("batched categorical inference retains near-pooled and shifted contrasts", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(930L)
  n <- 180L
  x <- data.frame(x = rnorm(n), u = factor(sample.int(3L, n, TRUE)),
                  o = ordered(sample.int(4L, n, TRUE)))
  y <- sin(x$x) + rnorm(n)
  for (near in c(FALSE, TRUE)) {
    bw <- npregbw(xdat = x, ydat = y, bandwidth.compute = FALSE,
                  bws = c(.45, .18, if (near) 1 - 1e-12 else .24))
    response <- cbind(y, rev(y), y + 1e10)
    mask <- np:::.np_npsig_structure(bw, x, 3L)
    scores <- np:::.np_npsig_scalar_batch(bw, x, 3L, response, mask, "test")
    expected <- vapply(seq_len(ncol(response)), function(j) {
      fit <- npreg(bws = bw, txdat = x, tydat = response[, j], gradients = TRUE, se = TRUE)
      np:::.np_npsig_statistic(fit, 3L, TRUE, mask)
    }, numeric(1L))
    expect_length(scores, ncol(response))
    expect_equal(scores, expected, tolerance = 2e-12)
  }
})

test_that("categorical bootstrap batching follows the fitted polynomial degree", {
  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  eligible <- np:::.np_npsig_scalar_batch_eligible
  # Use canonical engine metadata as carried by mixed-data bandwidth objects.
  bw <- list(ncon = 2L, type = "fixed", ckertype = "gaussian", ckerorder = 2L,
             ckerbound = "none", regtype.engine = "lp", degree.engine = c(0L, 0L))
  expect_true(eligible(bw, 120L))
  bw$degree.engine <- c(0L, 1L)
  expect_false(eligible(bw, 120L))
  bw$degree.engine <- c(1L, 1L)
  expect_false(eligible(bw, 120L))
})

test_that("OLS projections and shared kernel rows preserve individual CMS draws", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(42)
  n <- 240L
  x <- data.frame(x = rep(rnorm(n / 2L), each = 2L))
  y <- x$x + rnorm(n)
  fit <- lm(y ~ x, data = x, x = TRUE, y = TRUE)
  design <- np:::.np_cms_bootstrap_qr(fit)
  expect_s3_class(design, "qr")
  Ys <- matrix(rnorm(n * 9L), n)
  Ys[, 1L] <- Ys[, 1L] + 1e10
  scores <- vapply(seq_len(ncol(Ys)), function(j) {
    oracle <- np:::.np_cms_refit_residuals(fit, Ys[, j])
    actual <- np:::.np_cms_refit_residuals(fit, Ys[, j], design.qr = design)
    expect_equal(as.numeric(actual), as.numeric(oracle), tolerance = 2e-12)
    actual
  }, numeric(n))
  b <- npregbw(xdat = x, ydat = y, bws = .4, bandwidth.compute = FALSE)
  batch <- np:::.np_cms_statistics_batch
  inputs <- list(xdat = x, score = scores, bw = np:::kbandwidth(b),
                 fhat = seq(.8, 1.2, length.out = n), prodh = .4, pivot = TRUE)
  expect_equal(do.call(batch, c(inputs, list(dense.batch = TRUE))), do.call(batch, inputs),
               tolerance = 2e-12)
  expect_null(np:::.np_cms_bootstrap_qr(update(fit, weights = 1 + abs(x))))
})
