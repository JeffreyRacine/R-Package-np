test_that("LC derivative hats retain unsupported external rows without masking failures", {
  withr::local_options(np.messages = FALSE)
  withr::local_preserve_seed()
  set.seed(91831)
  x <- data.frame(x = seq(-1, 1, length.out = 48L),
                  z = seq(-1, 1, length.out = 48L))
  y <- sin(x$x) + rnorm(48L, sd = .1)
  ex <- rbind(x[c(4, 22, 41), ], data.frame(x = -10, z = 10))
  strict <- getFromNamespace(".npreghat_complete", "npRmpi")
  for (tree in c(FALSE, TRUE))
    for (type in c("fixed", "generalized_nn", "adaptive_nn"))
      for (regtype in c("lc", "lp")) {
        withr::local_options(np.tree = tree)
        bw <- npregbw(xdat = x, ydat = y, bws = rep(if (type == "fixed") .2 else 7, 2),
          bandwidth.compute = FALSE, regtype = regtype,
          degree = if (regtype == "lp") c(0L, 0L) else NULL,
          bwtype = type, ckertype = "epanechnikov")
        H <- suppressWarnings(npreghat(bw, txdat = x, exdat = ex,
                                      s = c(1L, 0L), output = "matrix",
                                      .np.defer.empty.rows = TRUE))
        oracle <- suppressWarnings(gradients(npreg(bw, txdat = x, tydat = y,
                                      exdat = ex, gradients = TRUE))[, 1L])
        expect_equal(as.vector(H %*% y), oracle, tolerance = 1e-10)
        one <- lapply(seq_len(nrow(ex)), function(i)
          suppressWarnings(npreghat(bw, txdat = x, exdat = ex[i, , drop = FALSE],
                                    s = c(1L, 0L), output = "matrix")))
        expect_equal(as.numeric(H), as.numeric(do.call(rbind, one)), tolerance = 1e-12)
        Y <- cbind(y, 2*y + 1)
        actual <- suppressWarnings(npreghat(bw, txdat = x, exdat = ex,
                         y = Y, s = c(1L, 0L), output = "apply"))
        expect_equal(as.numeric(actual), as.numeric(H %*% Y), tolerance = 1e-12)
        empty <- is.na(oracle)
        if (type != "generalized_nn") expect_true(empty[4L])
        if (any(empty)) {
          expect_true(all(is.na(H[empty, , drop = FALSE])))
          expect_identical(attr(H, ".np.empty.rows"), as.integer(empty))
          expect_error(strict(bw, txdat = x, exdat = ex, s = c(1L, 0L),
                              output = "matrix"), "hat helper failed")
        }
        expect_true(all(is.finite(H[!empty, , drop = FALSE])))
      }
  # A repeated-coordinate zero radius is invalid geometry, not an empty row.
  repeated <- data.frame(x = c(rep(0, 24), seq(.1, 1, length.out = 24)), z = x$z)
  bw <- npregbw(xdat = repeated, ydat = y, bws = c(2, 2),
                bandwidth.compute = FALSE, bwtype = "adaptive_nn")
  expect_error(npreghat(bw, txdat = repeated, exdat = ex,
                       s = c(1L, 0L), output = "matrix"), "[Zz]ero|radius")
})

test_that("wild derivative reuse propagates only certified empty rows", {
  withr::local_options(np.messages = FALSE)
  withr::local_preserve_seed()
  x <- data.frame(x = seq(-1, 1, length.out = 48L))
  y <- sin(2*x$x)
  ex <- data.frame(x = c(0, 10, -.3))
  bw <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE,
                ckertype = "epanechnikov")
  helper <- getFromNamespace(".np_wild_boot_from_regression_exact", "npRmpi")
  pilot <- fitted(npreg(bw, txdat = x, tydat = y))
  set.seed(852); w <- matrix(ifelse(runif(48*7) <= .5, -1, 1), 48, 7)
  seed <- .Random.seed
  expected <- t(vapply(1:7, function(j) suppressWarnings(gradients(npreg(bw,
    txdat = x, tydat = pilot + (y-pilot)*w[,j], exdat = ex, gradients = TRUE))[,1]),
    numeric(3)))
  for (threshold in c(Inf, 1)) {
    withr::local_options(np.plot.wild.apply.operator.threshold.bytes = threshold,
                        np.plot.wild.hat.block.bytes = 8*48)
    set.seed(852)
    a <- helper(xdat = x, exdat = ex, bws = bw, ydat = y, fit.mean.train = pilot,
                B = 7, gradients = TRUE)
    expect_equal(a$t, expected, tolerance = 1e-10)
    expect_true(is.na(a$t0[2]))
    expect_identical(.Random.seed, seed)
  }
})

test_that("MPI wild partial-row handling does not accept uncertified nonfinite output", {
  withr::local_options(np.messages = FALSE,
                      np.plot.wild.apply.operator.threshold.bytes = Inf)
  x <- data.frame(x = seq(-1, 1, length.out = 48L)); y <- sin(x$x)
  ex <- data.frame(x = c(0, 10, -.3))
  bw <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE,
                ckertype = "epanechnikov")
  helper <- getFromNamespace(".np_wild_boot_from_regression_exact", "npRmpi")
  for (bad in c("uncertified", "infinite")) local({
    value <- bad
    testthat::local_mocked_bindings(.np_plot_boot_from_hat_wild = function(H, ...) {
      force(H) # Complete canonical construction and its certificate first.
      t0 <- c(0, NA_real_, 0); t <- matrix(rep(t0, each = 7), 7, 3)
      if (value == "uncertified") t[1,1] <- NA_real_ else t[1,2] <- Inf
      list(t0 = t0, t = t)
    }, .package = "npRmpi")
    expect_error(helper(xdat = x, exdat = ex, bws = bw, ydat = y,
      fit.mean.train = y*.8, B = 7, gradients = TRUE), "non-finite")
  })
})
