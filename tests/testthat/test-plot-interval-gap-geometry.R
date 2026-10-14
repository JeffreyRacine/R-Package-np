# Renderer-only tests: no fitting, MPI, randomness or numerical-policy changes.
.gap_interval_package <- "np"
.gap_interval_helper <- function(name)
  getFromNamespace(name, .gap_interval_package)

.gap_interval_capture <- function(code) {
  captured <- list()
  testthat::with_mocked_bindings(
    force(code),
    draw.error.bands = function(ex, ely, ehy, ...) {
      captured[[length(captured) + 1L]] <<-
        list(style = "band", x = ex, lower = ely, upper = ehy)
    },
    draw.error.bars = function(ex, ely, ehy, ...) {
      captured[[length(captured) + 1L]] <<-
        list(style = "bar", x = ex, lower = ely, upper = ehy)
    },
    .np_plot_draw_all_band_legend = function(...) invisible(NULL),
    .np_plot_draw_variability_annotation = function(...) invisible(FALSE),
    .package = .gap_interval_package
  )
  captured
}

test_that("GAP01 continuous intervals retain exact tuples and only adjacent finite edges", {
  draw <- .gap_interval_helper(".np_plot_draw_errors_with_gaps")
  x <- as.double(1:8)
  lower <- c(-Inf, -.2, NA, -.4, -.5, -.6, NA, -.8)
  upper <- c(Inf, .2, .3, .4, .5, .6, NA, NaN)
  got <- .gap_interval_capture(draw(x, lower, upper, "band", "|", 8L, 2L, col = "black"))
  expect_length(got, 1L)
  keep <- is.finite(x) & is.finite(lower) & is.finite(upper)
  expect_identical(got[[1L]]$x[keep], x[keep])
  expect_identical(got[[1L]]$lower[keep], lower[keep])
  expect_identical(got[[1L]]$upper[keep], upper[keep])
  expect_true(all(is.na(got[[1L]]$x[!keep])))
  expect_true(all(is.na(got[[1L]]$lower[!keep])))
  expect_true(all(is.na(got[[1L]]$upper[!keep])))
  finite <- is.finite(got[[1L]]$x)
  edges <- which(head(finite, -1L) & tail(finite, -1L))
  expect_identical(edges, c(4L, 5L))
  all.finite <- .gap_interval_capture(draw(x, -x, x, "band", "|", 8L, 2L, col = "black"))
  expect_identical(all.finite[[1L]], list(style = "band", x = x, lower = -x, upper = x))
  empty <- .gap_interval_capture(draw(x, rep(NA_real_, 8), x,
                                      "band", "|", 8L, 2L, col = "black"))
  expect_length(empty, 0L)
})

test_that("GAP02 categorical bars preserve the original finite coordinate positions", {
  draw <- .gap_interval_helper(".np_plot_draw_errors_with_gaps")
  x <- c(1, 2, 3, 4, 5)
  lower <- c(-1, NA, -3, -4, -5)
  upper <- c(1, 2, 3, Inf, 5)
  got <- .gap_interval_capture(draw(x, lower, upper, "bar", "I", 5L, 1L, col = "black"))
  expect_identical(got, list(list(style = "bar", x = x[c(1L, 3L, 5L)],
    lower = lower[c(1L, 3L, 5L)], upper = upper[c(1L, 3L, 5L)])))
})

test_that("GAP03 estimate-error conversion preserves missing and mismatched row positions", {
  vectors <- .gap_interval_helper(".np_plot_error_vectors_with_gaps")
  x <- c(10, 20, 30, 40, 50)
  estimate <- c(1, NA, 3, 4, 5)
  err <- cbind(c(.1, .2, .3, NA, .5), c(.2, .3, .4, .5, .6), NA_real_)
  got <- vectors(x, estimate, err)
  expect_identical(got$x, c(10, NA_real_, 30, NA_real_, 50))
  expect_identical(got$lower[c(1, 3, 5)], estimate[c(1, 3, 5)] - err[c(1, 3, 5), 1L])
  expect_identical(got$upper[c(1, 3, 5)], estimate[c(1, 3, 5)] + err[c(1, 3, 5), 2L])
  expect_identical(got$finite, c(TRUE, FALSE, TRUE, FALSE, TRUE))
})

test_that("GAP04 every all-band layer retains its own gaps and validated ragged rows", {
  draw <- .gap_interval_helper(".np_plot_draw_all_errors_with_gaps")
  center <- c(.2, .4, .6, .8)
  err <- cbind(c(.1, NA, .2, .3), c(.15, .25, .35, .45))
  whole <- cbind(rep(.2, 4), rep(.3, 4))
  bands <- list(pointwise = err, simultaneous = matrix(NA_real_, 4, 2), bonferroni = whole)
  got <- .gap_interval_capture(draw(1:4, center, bands, add.legend = FALSE,
                                    plot.errors.method = "asymptotic"))
  expect_length(got, 2L)
  expect_identical(got[[1L]]$x, c(1, NA_real_, 3, 4))
  expect_identical(got[[2L]]$x, as.double(1:4))
  expect_identical(got[[1L]]$lower[c(1, 3, 4)], (center - err[, 1L])[c(1, 3, 4)])
  padded <- .gap_interval_capture(draw(1:4, c(center, NA, NA), bands,
    add.legend = FALSE, plot.errors.method = "asymptotic"))
  expect_identical(padded, got)
  expect_error(.gap_interval_capture(draw(1:4, c(center, 9), bands,
    add.legend = FALSE, plot.errors.method = "asymptotic")), "rows do not match")
  bars <- .gap_interval_capture(draw(1:4, center, bands, xi.factor = TRUE,
    add.legend = FALSE, plot.errors.method = "asymptotic"))
  expect_identical(bars[[1L]]$style, "bar")
  expect_identical(bars[[1L]]$x, c(1, 3, 4))
})

test_that("GAP05 conmode continuous intervals break while categorical intervals stay segments", {
  panel <- .gap_interval_helper(".np_plot_conmode_panel")
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  captured.lines <- captured.segments <- list()
  capture.panel <- function(dat) {
    captured.lines <<- captured.segments <<- list()
    testthat::with_mocked_bindings(
      panel(dat, level = "chosen"),
      lines = function(x, y, ...) {
        captured.lines[[length(captured.lines) + 1L]] <<- list(x = x, y = y)
      },
      segments = function(x0, y0, x1, y1, ...) {
        captured.segments[[length(captured.segments) + 1L]] <<-
          list(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
      }, .package = "graphics"
    )
  }
  dat <- data.frame(x = c(3, 1, 2, 4, 5), variable = "x", view = "fixed",
                    probability = c(.3, .1, NA, .4, .5),
                    lower = c(.25, .05, NA, .35, .45),
                    upper = c(.35, .15, NA, .45, .55))
  capture.panel(dat)
  expect_length(captured.lines, 2L)
  expect_identical(captured.lines[[1L]]$x, c(1, NA_real_, 3, 4, 5))
  expect_identical(captured.lines[[1L]]$y, c(.05, NA_real_, .25, .35, .45))
  expect_length(captured.segments, 0L)
  dat <- data.frame(x = factor(c("a", "b", "c")), variable = "x", view = "fixed",
                    probability = c(.2, NA, .6), lower = c(.1, NA, .5),
                    upper = c(.3, NA, .7))
  capture.panel(dat)
  expect_length(captured.lines, 0L)
  expect_length(captured.segments, 1L)
  expect_identical(captured.segments[[1L]]$x0, c(1, 3))
  expect_identical(captured.segments[[1L]]$y0, c(.1, .5))
  expect_identical(captured.segments[[1L]]$y1, c(.3, .7))
  dat$lower <- dat$upper <- NA_real_
  capture.panel(dat)
  expect_length(captured.lines, 0L)
  expect_length(captured.segments, 0L)
})
