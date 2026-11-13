test_that("interval caps use display geometry and retain vertical endpoints", {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  draw <- get("draw.error.bars", envir = environment(npreg))
  records <- list()
  capture <- new.env(parent = environment(draw))
  capture$lines <- function(x, y, ...) {
    records[[length(records) + 1L]] <<- list(x = x, y = y)
  }
  environment(draw) <- capture
  widths <- function(multiplier, log = "", reverse = FALSE, ex = c(1, 2)) {
    records <<- list()
    yl <- c(.01, 1) * multiplier
    xl <- c(.5, 2.5)
    if (reverse) { yl <- rev(yl); xl <- rev(xl) }
    graphics::plot.new(); graphics::plot.window(xlim = xl, ylim = yl, log = log)
    lo <- rep(.1, length(ex)) * multiplier
    hi <- rep(.8, length(ex)) * multiplier
    draw(ex, lo, hi)
    expect_identical(records[[1]]$x[seq_along(ex) * 3 - 2], ex)
    expect_identical(records[[1]]$y[seq_along(ex) * 3 - 2], lo)
    expect_identical(records[[1]]$y[seq_along(ex) * 3 - 1], hi)
    caps <- records[[2]]$x
    abs(diff(graphics::grconvertX(caps[1:2], from = "user", to = "inches")))
  }
  for (log in c("", "x", "y", "xy")) for (reverse in c(FALSE, TRUE)) {
    a <- widths(1, log, reverse)
    expect_gt(a, 0)
    expect_equal(widths(.001, log, reverse), a, tolerance = 1e-12)
  }
  expect_gt(widths(1, ex = 1), 0)
  expect_gt(widths(1, ex = c(1, 1)), 0)
  records <- list()
  draw(c(1, NA_real_, 2), c(0, NA_real_, 0), c(0, NA_real_, 1), hbar = FALSE)
  expect_length(records, 1)
})

test_that("independent-scale partially linear overlays use panel-local data", {
  withr::local_options(np.messages = FALSE)
  set.seed(19313)
  d <- data.frame(x = runif(30, -1, 1), z = rnorm(30), o = ordered(rep(1:3, 10)))
  d$y <- d$z + d$x + rnorm(30, sd = .3)
  b <- npplregbw(y ~ z | x + o, data = d,
    bws = matrix(c(.6, .6, .3, .3), 2, 2), bandwidth.compute = FALSE)
  f <- npplreg(b)
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  for (object in list(b, f)) for (common in c(FALSE, TRUE))
    expect_no_error(plot(object, errors = "none", common.scale = common,
                         data.overlay = TRUE, neval = 4, perspective = FALSE))
})
