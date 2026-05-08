library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

test_that("npindex frozen mean helper stays on the exact single-index scale", {
  set.seed(42)
  n <- 200L
  B <- 25L
  x <- runif(n, -1, 1)
  z <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))
  xdat <- data.frame(x = x, z = z)

  fit <- quiet_eval(
    npindex(
      y ~ x + z,
      bwtype = "adaptive_nn",
      nmulti = 1L
    )
  )

  eval_grid_fun <- getFromNamespace(".np_plot_singleindex_eval_grid", "np")
  to_frame <- getFromNamespace("toFrame", "np")
  exact_helper <- getFromNamespace(".np_inid_boot_from_index", "np")

  eval.info <- eval_grid_fun(
    bws = fit$bws,
    xdat = to_frame(xdat),
    neval = 40L,
    trim = 0.0
  )
  counts <- stats::rmultinom(B, size = n, prob = rep.int(1 / n, n))

  frozen <- exact_helper(
    xdat = xdat,
    ydat = y,
    bws = fit$bws,
    B = B,
    counts = counts,
    frozen = TRUE,
    idx.eval = eval.info$idx.eval,
    progress.label = "test frozen"
  )
  exact <- exact_helper(
    xdat = xdat,
    ydat = y,
    bws = fit$bws,
    B = B,
    counts = counts,
    frozen = FALSE,
    idx.eval = eval.info$idx.eval,
    progress.label = "test exact"
  )

  sd.ratio <- stats::na.omit(apply(frozen$t, 2L, stats::sd) / apply(exact$t, 2L, stats::sd))

  expect_equal(frozen$t0, exact$t0, tolerance = 1e-10)
  expect_true(length(sd.ratio) > 0L)
  expect_gt(stats::median(sd.ratio), 0.75)
  expect_lt(stats::median(sd.ratio), 1.25)
})

test_that("npindex public frozen plot-data mean stays on the exact scale", {
  set.seed(42)
  n <- 120L
  x <- runif(n, -1, 1)
  z <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))

  fit <- quiet_eval(
    npindex(
      y ~ x + z,
      bwtype = "adaptive_nn",
      nmulti = 1L
    )
  )

  get_obj <- function(mode) {
    suppressWarnings(plot(
      fit,
      output = "data",
      neval = 40L,
      errors = "bootstrap",
      bootstrap = "inid",
      boot_control = np_boot_control(nonfixed = mode),
      B = 39L,
      band = "pointwise"
    ))[[1L]]
  }

  frozen <- get_obj("frozen")
  exact <- get_obj("exact")
  ratio <- stats::median(abs(exact$merr[, 1L]) / pmax(abs(frozen$merr[, 1L]), 1e-12), na.rm = TRUE)

  expect_equal(frozen$mean, exact$mean, tolerance = 1e-12)
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("npindex nonfixed frozen bootstrap still supports gradient slices", {
  set.seed(42)
  n <- 80L
  x <- runif(n)
  z <- rnorm(n)
  y <- x + rnorm(n)

  fit <- quiet_eval(
    npindex(
      y ~ x + z,
      regtype = "ll",
      bwtype = "adaptive_nn",
      nmulti = 1
    )
  )

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    unlink(tf)
  }, add = TRUE)

  expect_no_error(
    capture.output(plot(
      fit,
      neval = 20L,
      errors = "bootstrap",
      bootstrap = "inid",
      boot_control = np_boot_control(nonfixed = "frozen"),
      B = 41L,
      band = "pointwise",
      gradients = TRUE
    ))
  )
})
