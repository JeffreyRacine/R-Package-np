test_that("categorical GNN fitted tiles retain training identity across LP owners", {
  n <- 36L
  xdat <- data.frame(
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high"), each = 12L),
                levels = c("low", "mid", "high")),
    x = seq(-1, 1, length.out = n)
  )
  ydat <- sin(2 * xdat$x) + 0.3 * (xdat$u == "b") -
    0.2 * (xdat$o == "low") + cos(seq_len(n)) * 0.05
  donor <- cbind(c(2:n, 1L), n:1L)
  null.mean <- rep(mean(ydat), n)
  residual.pool <- ydat - mean(ydat)
  response <- null.mean + matrix(residual.pool[donor], nrow = n)
  specifications <- list(
    list(regtype = "lc"),
    list(regtype = "ll"),
    list(regtype = "lp", degree = 2L, bernstein.basis = FALSE),
    list(regtype = "lp", degree = 2L, bernstein.basis = TRUE)
  )
  for (specification in specifications) {
    bw <- do.call(npregbw, c(list(
      xdat = xdat, ydat = ydat, bws = c(0.18, 0.2, 10),
      bandwidth.compute = FALSE, bwtype = "generalized_nn"
    ), specification))
    for (index in 1:2) {
      set.seed(917)
      before <- .Random.seed
      tile <- np:::.np_npsig_streamed_iid_tile(
        bw, xdat, index, donor.index = donor,
        null.mean = null.mean, residual.pool = residual.pool,
        pivotal = FALSE
      )
      expect_identical(.Random.seed, before)
      ready <- np:::.np_npsig_streamed_iid_tile(
        bw, xdat, index, response.matrix = response,
        null.mean = null.mean, residual.pool = residual.pool,
        pivotal = FALSE
      )
      expect_identical(tile, ready)
      oracle <- vapply(seq_len(ncol(response)), function(j) {
        fit <- npreg(bws = bw, txdat = xdat, tydat = response[, j],
                     gradients = TRUE, se = FALSE)
        np:::.np_npsig_statistic(fit, index, FALSE)
      }, numeric(1L))
      expect_equal(tile, oracle,
                   tolerance = if (specification$regtype == "lp") 1e-9 else 2e-10)
    }
  }
})

test_that("GNN training and external radii retain distinct tie boundaries", {
  n <- 12L
  for (multiplicity in 3:5) {
    xdat <- data.frame(
      group = factor(rep(c("a", "b"), length.out = n)),
      tied = c(rep(0, multiplicity), seq_len(n - multiplicity))
    )
    ydat <- 1 + sin(seq_len(n)) / 4
    bw <- npregbw(xdat = xdat, ydat = ydat, regtype = "lc",
                  bwtype = "generalized_nn", bws = c(0.2, 4),
                  bandwidth.compute = FALSE)
    set.seed(918)
    before <- .Random.seed
    tile <- tryCatch(np:::.np_npsig_streamed_iid_tile(
      bw, xdat, 1L, response.matrix = matrix(ydat, ncol = 1L),
      null.mean = ydat, residual.pool = ydat, pivotal = FALSE
    ), error = identity)
    expect_identical(.Random.seed, before)
    if (multiplicity <= 4L) {
      fit <- npreg(bws = bw, txdat = xdat, tydat = ydat,
                   gradients = TRUE, se = FALSE)
      expect_equal(tile, np:::.np_npsig_statistic(fit, 1L, FALSE),
                   tolerance = 2e-10)
    } else {
      # The serial private bridge retains its incumbent generic failure.
      expect_s3_class(tile, "error")
      training <- tryCatch(npreg(
        bws = bw, txdat = xdat, tydat = ydat, se = FALSE
      ), error = identity)
      expect_s3_class(training, "np_nn_zero_radius")
      expect_identical(training$variable, "tied")
      expect_identical(training$lookup.k, 4L)
      expect_identical(training$excluded, 1L)
    }
    external <- tryCatch(npreg(
      bws = bw, txdat = xdat, tydat = ydat, exdat = xdat, se = FALSE
    ), error = identity)
    if (multiplicity < 4L) {
      expect_true(all(is.finite(external$mean)))
    } else {
      expect_s3_class(external, "np_nn_zero_radius")
      expect_identical(external$variable, "tied")
      expect_identical(external$lookup.k, 4L)
      expect_identical(external$excluded, 0L)
    }
  }
})
