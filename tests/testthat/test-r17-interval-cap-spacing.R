test_that("interval caps depend only on distinct horizontal display spacing", {
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw <- get("draw.error.bars", envir = environment(npreg))
  records <- list()
  capture <- new.env(parent = environment(draw))
  capture$lines <- function(x, y, ...) {
    records[[length(records) + 1L]] <<- list(x = x, y = y)
  }
  environment(draw) <- capture
  width <- function(ylim, lo, hi, x = c(1, 2), log = "", reverse = FALSE) {
    records <<- list()
    plot.new()
    plot.window(xlim = if (reverse) c(2.5, .5) else c(.5, 2.5),
                ylim = if (reverse) rev(ylim) else ylim, log = log)
    draw(x, rep(lo, length(x)), rep(hi, length(x)))
    expect_identical(records[[1]]$y[seq_along(x)*3-2], rep(lo, length(x)))
    expect_identical(records[[1]]$y[seq_along(x)*3-1], rep(hi, length(x)))
    abs(diff(grconvertX(records[[2]]$x[1:2], "user", "inches")))
  }
  for (log in c("", "x", "y", "xy")) for (reverse in c(FALSE, TRUE)) {
    ref <- width(c(.01, 5), 1, 2, log = log, reverse = reverse)
    expect_gt(ref, .01)
    expect_equal(width(c(.001, 5000), 1, 1.00001, log = log, reverse = reverse),
                 ref, tolerance = 1e-12)
    expect_equal(width(c(.01, 5), 1, 2, x = c(1, 1, 2, 2),
                       log = log, reverse = reverse), ref, tolerance = 1e-12)
    expect_equal(width(c(.01, 5), 1, 1, log = log, reverse = reverse), 0)
  }
  expect_gt(width(c(0, 5), 1, 2, x = 1), .01)
})

test_that("coefficient cap unification preserves the incumbent graphics style", {
  withr::local_options(np.messages = FALSE)
  withr::local_preserve_seed()
  set.seed(2)
  d <- data.frame(x = runif(40), z = runif(40))
  d$y <- d$x + d$z + rnorm(40)
  b <- npplregbw(y ~ x | z, data = d, bws = matrix(c(.4, .4), 2, 1),
                 bandwidth.compute = FALSE)
  f <- npplreg(b, se = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(col = "red", fg = "blue", lty = 2)
  seen <- NULL
  pkg <- getNamespaceName(environment(npreg))
  draw <- getFromNamespace("draw.error.bars", pkg)
  testthat::local_mocked_bindings(draw.error.bars = function(...) {
    seen <<- list(...)
    draw(...)
  }, .package = pkg)
  plot(f, coef = TRUE, errors = "asymptotic")
  expect_identical(seen$col, "blue")
  expect_identical(seen$lty, "dashed")
})
