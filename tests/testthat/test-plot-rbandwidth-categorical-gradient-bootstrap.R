test_that("wild categorical regression gradient helper matches explicit refits", {
  skip_if_not_installed("np")

  library(np)

  wild_helper <- getFromNamespace(".np_plot_boot_from_hat_wild_factor_effects", "np")
  rad_draws <- getFromNamespace(".np_rademacher_draws", "np")

  set.seed(20260312)
  n <- 36L
  g <- factor(sample(c("a", "b"), n, replace = TRUE))
  x <- runif(n)
  y <- 1 + 0.5 * (g == "b") + sin(2 * pi * x) + rnorm(n, sd = 0.05)
  xdat <- data.frame(g = g, x = x)
  bw <- npregbw(
    xdat = xdat,
    ydat = y,
    regtype = "ll",
    bwtype = "fixed",
    bws = c(0.25, 0.3),
    bandwidth.compute = FALSE
  )

  exdat <- plot(
    bw,
    xdat = xdat,
    ydat = y,
    gradients = TRUE,
    output = "data",
    perspective = FALSE
  )[[1L]]$eval
  H <- npreghat(
    bws = bw,
    txdat = xdat,
    exdat = exdat,
    output = "matrix"
  )
  fit.mean <- as.vector(npreghat(
    bws = bw,
    txdat = xdat,
    exdat = xdat,
    y = y,
    output = "apply"
  ))

  B <- 7L
  set.seed(11)
  helper.out <- wild_helper(
    H = H,
    ydat = y,
    fit.mean = fit.mean,
    B = B,
    wild = "rademacher"
  )

  set.seed(11)
  draws <- rad_draws(n = n, B = B)
  explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(exdat))
  for (b in seq_len(B)) {
    ystar <- fit.mean + (y - fit.mean) * draws[, b]
    mean.b <- npreg(
      txdat = xdat,
      tydat = ystar,
      exdat = exdat,
      bws = bw,
      gradients = FALSE,
      warn.glp.gradient = FALSE
    )$mean
    explicit.t[b, ] <- mean.b - mean.b[1L]
  }
  fit0 <- npreg(
    txdat = xdat,
    tydat = y,
    exdat = exdat,
    bws = bw,
    gradients = FALSE,
    warn.glp.gradient = FALSE
  )$mean

  expect_equal(helper.out$t, explicit.t, tolerance = 1e-6)
  expect_equal(as.vector(helper.out$t0), as.vector(fit0 - fit0[1L]), tolerance = 1e-6)
})

test_that("categorical regression gradient bootstrap works for default, inid, and wild routes", {
  skip_if_not_installed("np")

  library(np)

  set.seed(20260312)
  n <- 36L
  g <- factor(sample(c("a", "b"), n, replace = TRUE))
  x <- runif(n)
  y <- 1 + 0.4 * (g == "b") + cos(2 * pi * x) + rnorm(n, sd = 0.05)
  xdat <- data.frame(g = g, x = x)
  bw <- npregbw(
    xdat = xdat,
    ydat = y,
    regtype = "ll",
    bwtype = "fixed",
    bws = c(0.25, 0.3),
    bandwidth.compute = FALSE
  )
  fit <- npreg(
    bws = bw,
    txdat = xdat,
    tydat = y,
    gradients = TRUE,
    errors = TRUE
  )

  for (boot.method in c("default", "inid", "wild")) {
    args <- list(
      fit,
      output = "data",
      perspective = FALSE,
      gradients = TRUE,
      errors = "bootstrap",
      B = 11L
    )
    if (!identical(boot.method, "default"))
      args$bootstrap <- boot.method

    out <- suppressWarnings(do.call(plot, args))
    expect_type(out, "list")
    expect_true(length(out) >= 1L, info = boot.method)
    expect_true(length(out[[1L]]$bxp) > 0L, info = boot.method)
    expect_equal(length(out[[1L]]$bxp$names), 2L, info = boot.method)
    expect_true(all(is.finite(out[[1L]]$grad)), info = boot.method)
    expect_true(all(is.finite(out[[1L]]$glerr)), info = boot.method)
    expect_true(all(is.finite(out[[1L]]$gherr)), info = boot.method)
  }
})
