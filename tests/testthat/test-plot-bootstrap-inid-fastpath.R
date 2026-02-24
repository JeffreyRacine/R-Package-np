test_that("inid lc fast path matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(321)
  n <- 40
  x <- runif(n)
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 21))

  bw <- npregbw(
    y ~ x,
    regtype = "lc",
    bws = c(0.2),
    bandwidth.compute = FALSE
  )

  H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  B <- 13L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  fast.out <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    fit.b <- npreg(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = ex,
      bws = bw,
      gradients = FALSE,
      warn.glp.gradient = FALSE
    )
    explicit.t[b, ] <- fit.b$mean
  }

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-10)
  expect_equal(fast.out$t0, as.vector(H %*% y), tolerance = 1e-12)
})

test_that("inid lc fast path toggle preserves plot bootstrap contract", {
  skip_if_not_installed("np")

  set.seed(322)
  n <- 60
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.15)
  bw <- npregbw(y ~ x, regtype = "lc", nmulti = 1)

  run_plot <- function(disable) {
    old <- getOption("np.plot.inid.fastpath.disable")
    on.exit(options(np.plot.inid.fastpath.disable = old), add = TRUE)
    options(np.plot.inid.fastpath.disable = disable)

    suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 9
      )
    )
  }

  out.fast <- run_plot(disable = FALSE)
  out.legacy <- run_plot(disable = TRUE)

  expect_type(out.fast, "list")
  expect_true(length(out.fast) > 0)
  expect_type(out.legacy, "list")
  expect_true(length(out.legacy) > 0)
})
