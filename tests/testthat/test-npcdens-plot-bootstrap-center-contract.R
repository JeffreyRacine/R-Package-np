test_that("npcdens plot-data centers stay on predict and absolute bootstrap endpoints remain ordered", {
  skip_if_not_installed("np")

  check_mode <- function(bwtype, nonfixed_mode) {
    set.seed(42)
    n <- 100L
    x <- runif(n, -1, 1)
    y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))

    bw <- npcdensbw(y ~ x, nmulti = 1L, bwtype = bwtype)
    fit <- npcdens(bws = bw)

    out <- suppressWarnings(plot(
      fit,
      output = "data",
      view = "fixed",
      perspective = FALSE,
      neval = 40L,
      errors = "bootstrap",
      bootstrap = "inid",
      B = 39L,
      band = "pointwise",
      boot_control = np_boot_control(nonfixed = nonfixed_mode)
    ))

    expect_s3_class(out, "np_plot_data")
    for (j in seq_along(out$panels)) {
      pred <- predict(fit, exdat = plot_coordinates(out, j, "X"), eydat = plot_coordinates(out, j, "Y"))
      d <- plot_target(out, j, band = "pointwise")
      expect_equal(d$estimate, as.numeric(pred), tolerance = 1e-8)
      expect_true(all(d$lower <= d$upper + 1e-10, na.rm = TRUE))
    }

    out
  }

  fixed.exact <- check_mode("fixed", "exact")
  generalized.exact <- check_mode("generalized_nn", "exact")
  generalized.frozen <- check_mode("generalized_nn", "frozen")
  adaptive.exact <- check_mode("adaptive_nn", "exact")
  adaptive.frozen <- check_mode("adaptive_nn", "frozen")

  expect_equal(
    lapply(generalized.exact$panels, function(p) p$values$estimate),
    lapply(generalized.frozen$panels, function(p) p$values$estimate),
    tolerance = 1e-12
  )
  expect_equal(
    lapply(adaptive.exact$panels, function(p) p$values$estimate),
    lapply(adaptive.frozen$panels, function(p) p$values$estimate),
    tolerance = 1e-12
  )
  expect_true(all(vapply(fixed.exact$panels, function(p) nrow(p$intervals) == nrow(p$grid), logical(1))))
})
