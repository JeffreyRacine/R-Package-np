test_that("duplicate endpoint comparison uses only local operand scales", {
  pkg <- "np"
  f <- getFromNamespace(".npcopula_display_values", pkg)
  g <- list(row.order = 1:4, display.rows = c(1L, 3L), display.cell = c(1L, 1L, 2L, 2L))
  eps <- .Machine$double.eps
  for (a in 2^c(-400, 0, 400)) {
    # Independent binary arithmetic: exact endpoint zero, one-ulp operand perturbation.
    c <- a * c(1, 1 + eps, 2, 2)
    d <- a * c(1, 1, 1, 1)
    e <- c - d
    expect_error(f(e, g), "conflicting")
    expect_identical(f(e, g, roundoff.scale = abs(c) + abs(d)), e[c(1L, 3L)])
    expect_error(f(e, g, exact = TRUE, roundoff.scale = abs(c) + abs(d)), "conflicting")
    for (k in c(31, 32)) {
      v <- a * c(0, k*eps*2, 1, 1)
      expect_identical(f(v, g, roundoff.scale = rep(2*a, 4)), v[c(1L, 3L)])
    }
    v <- a * c(0, 33*eps*2, 1, 1)
    expect_error(f(v, g, roundoff.scale = rep(2*a, 4)), "conflicting")
    # Another cell's enormous scale must not mask this conflict.
    expect_error(f(v, g, roundoff.scale = a*c(2, 2, 1e30, 1e30)), "conflicting")
    order <- c(4L, 2L, 3L, 1L)
    perm <- list(row.order = 1:4, display.rows = c(4L, 1L), display.cell = g$display.cell[order])
    expect_identical(f(e[order], perm, roundoff.scale = (abs(c)+abs(d))[order]), e[c(1L, 4L)])
  }
  expect_error(f(c(0, 0, 1, Inf), g, roundoff.scale = rep(2, 4)), "conflicting")
  expect_error(f(rep(0, 4), g, roundoff.scale = c(1, 2)), "operand scale")
  expect_error(f(rep(0, 4), g, roundoff.scale = c(1, -1, 1, 1)), "operand scale")
  expect_error(f(rep(0, 4), g, roundoff.scale = c(1, Inf, 1, 1)), "operand scale")
})

test_that("copula duplicate endpoints tolerate cancellation without changing data", {
  if (FALSE) {
    skip_on_cran()
    if (!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES", "1"))))
      stop("required MPI pool unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  set.seed(301)
  d <- data.frame(x = rnorm(180))
  d$y <- .4*d$x + rnorm(180)
  ax <- seq(.1, .9, length.out = 10L)
  u <- data.frame(u1 = ax[c(1:10, 3L)], u2 = rev(ax)[c(1:10, 3L)])
  bw <- npudensbw(dat = d, bws = c(10, 10), bandwidth.compute = FALSE,
                   bwtype = "generalized_nn")
  fit <- npcopula(bw, data = d, u = u, n.quasi.inv = 200L, se = TRUE)
  frozen <- serialize(fit, NULL)
  payloads <- list()
  for (band in if (identical(Sys.getenv("N5_FULL"), "true")) c("pointwise", "bonferroni", "all") else "pointwise") {
    rng <- .Random.seed
    expected <- plot(fit, errors = "asymptotic", band = band, view = "fixed", output = "data")
    payloads[[band]] <- expected
    actual <- plot(fit, errors = "asymptotic", band = band, view = "fixed", output = "plot-data")
    expect_identical(actual, expected)
    expect_identical(.Random.seed, rng)
    expect_equal(nrow(actual$panels[[1L]]$values), length(fitted(fit)))
  }
  if (nzchar(Sys.getenv("N5_SNAPSHOT"))) saveRDS(payloads, Sys.getenv("N5_SNAPSHOT"))
  expect_identical(serialize(fit, NULL), frozen)
  if (identical(Sys.getenv("N5_FULL"), "true")) {
    fixed <- npcopula(npudensbw(dat = d, bws = c(.4, .5), bandwidth.compute = FALSE),
                     data = d, u = u, n.quasi.inv = 200L, se = TRUE)
    a <- plot(fixed, errors = "asymptotic", band = "pointwise", view = "fixed", output = "data")
    b <- plot(fixed, errors = "asymptotic", band = "pointwise", view = "fixed", output = "plot-data")
    expect_identical(a, b)
  }
})
