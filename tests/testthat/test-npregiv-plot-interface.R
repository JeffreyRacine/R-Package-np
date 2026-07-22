iv_plot_package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"

iv_plot_object <- function() {
  structure(
    list(
      z = data.frame(z = c(-2, 0, 2)),
      y = c(-10, 0, 10),
      phi = c(-1.5, 0, 1.5),
      phi.deriv.1 = matrix(c(0.25, 0.5, 0.75), ncol = 1L),
      zeval = data.frame(z = c(-0.5, 0.5)),
      phi.eval = c(-1, 1),
      phi.deriv.eval.1 = matrix(c(0.4, 0.6), ncol = 1L)
    ),
    class = "npregiv"
  )
}

iv_deriv_plot_object <- function() {
  structure(
    list(
      z = data.frame(z = c(-2, 0, 2)),
      y = c(-10, 0, 10),
      phi = c(-1.5, 0, 1.5),
      phi.prime = c(0.25, 0.5, 0.75),
      zeval = data.frame(z = c(-0.5, 0.5)),
      phi.prime.eval = c(0.4, 0.6)
    ),
    class = "npregivderiv"
  )
}

iv_plot_spec <- function(object, family, gradients, data_overlay, data_rug,
                         dots = list()) {
  getFromNamespace(".np_iv_plot_spec", iv_plot_package)(
    object = object,
    family = family,
    gradients = gradients,
    data_overlay = data_overlay,
    data_rug = data_rug,
    dots = dots
  )
}

test_that("IV plot methods expose the canonical interface and class defaults", {
  expect_identical(
    names(formals(getS3method("plot", "npregiv"))),
    c("x", "gradients", "data_overlay", "data_rug", "...")
  )
  expect_identical(formals(getS3method("plot", "npregiv"))$gradients, FALSE)
  expect_identical(formals(getS3method("plot", "npregiv"))$data_overlay, TRUE)
  expect_identical(formals(getS3method("plot", "npregiv"))$data_rug, FALSE)

  expect_identical(
    names(formals(getS3method("plot", "npregivderiv"))),
    c("x", "gradients", "data_overlay", "data_rug", "...")
  )
  expect_identical(formals(getS3method("plot", "npregivderiv"))$gradients, TRUE)
  expect_identical(formals(getS3method("plot", "npregivderiv"))$data_overlay, TRUE)
  expect_identical(formals(getS3method("plot", "npregivderiv"))$data_rug, FALSE)
})

test_that("IV plot specifications select exact training and evaluation fields", {
  iv <- iv_plot_object()
  level <- iv_plot_spec(iv, "npregiv", FALSE, TRUE, FALSE)
  expect_true(level$evaluated)
  expect_identical(level$curve.args$x, iv$zeval$z)
  expect_identical(level$curve.args$y, iv$phi.eval)
  expect_identical(level$curve.args$xlim, range(c(iv$zeval$z, iv$z$z)))
  expect_identical(level$curve.args$ylim, range(c(iv$phi.eval, iv$y)))
  expect_true(level$overlay)

  gradient <- iv_plot_spec(iv, "npregiv", TRUE, TRUE, FALSE)
  expect_true(gradient$evaluated)
  expect_identical(gradient$curve.args$x, iv$zeval$z)
  expect_identical(gradient$curve.args$y,
                   as.vector(iv$phi.deriv.eval.1[, 1L]))
  expect_false(gradient$overlay)
  expect_identical(gradient$curve.args$ylim,
                   range(iv$phi.deriv.eval.1[, 1L]))

  deriv <- iv_deriv_plot_object()
  derivative <- iv_plot_spec(deriv, "npregivderiv", TRUE, TRUE, FALSE)
  expect_true(derivative$evaluated)
  expect_identical(derivative$curve.args$x, deriv$zeval$z)
  expect_identical(derivative$curve.args$y, deriv$phi.prime.eval)
  expect_false(derivative$overlay)

  structural <- iv_plot_spec(deriv, "npregivderiv", FALSE, TRUE, FALSE)
  expect_false(structural$evaluated)
  expect_identical(structural$curve.args$x, deriv$z$z)
  expect_identical(structural$curve.args$y, deriv$phi)
  expect_true(structural$overlay)
  expect_identical(structural$curve.args$ylim, range(c(deriv$phi, deriv$y)))
})

test_that("IV automatic limits include exactly the active drawable layers", {
  iv <- iv_plot_object()
  curve.only <- iv_plot_spec(iv, "npregiv", FALSE, FALSE, FALSE)
  expect_identical(curve.only$curve.args$xlim, range(iv$zeval$z))
  expect_identical(curve.only$curve.args$ylim, range(iv$phi.eval))

  rugged.gradient <- iv_plot_spec(iv, "npregiv", TRUE, TRUE, TRUE)
  expect_true(rugged.gradient$rug)
  expect_identical(rugged.gradient$curve.args$xlim,
                   range(c(iv$zeval$z, iv$z$z)))
  expect_identical(rugged.gradient$curve.args$ylim,
                   range(iv$phi.deriv.eval.1[, 1L]))

  user.xlim <- c(-100, 100)
  user.ylim <- c(-200, 200)
  explicit <- iv_plot_spec(
    iv, "npregiv", FALSE, TRUE, TRUE,
    dots = list(xlim = user.xlim, ylim = user.ylim)
  )
  expect_identical(explicit$curve.args$xlim, user.xlim)
  expect_identical(explicit$curve.args$ylim, user.ylim)

  positive <- iv
  positive$z$z <- c(-5, 2, 7)
  positive$zeval$z <- c(-1, 3)
  positive$y <- c(-20, 4, 12)
  positive$phi.eval <- c(-2, 6)
  logged <- iv_plot_spec(
    positive, "npregiv", FALSE, TRUE, TRUE,
    dots = list(log = "xy")
  )
  expect_identical(logged$curve.args$xlim, c(2, 7))
  expect_identical(logged$curve.args$ylim, c(4, 12))
  nonpositive <- iv
  nonpositive$phi.eval <- c(-2, -1)
  expect_error(
    iv_plot_spec(nonpositive, "npregiv", FALSE, FALSE, FALSE,
                 dots = list(log = "y")),
    "finite positive"
  )
})

test_that("IV plot specification is pure and rejects malformed stored geometry", {
  iv <- iv_plot_object()
  set.seed(1729)
  seed <- .Random.seed
  options.before <- options()
  device.before <- grDevices::dev.cur()
  invisible(iv_plot_spec(iv, "npregiv", FALSE, TRUE, TRUE))
  expect_identical(.Random.seed, seed)
  expect_identical(options(), options.before)
  expect_identical(grDevices::dev.cur(), device.before)

  missing.coordinate <- iv
  missing.coordinate$zeval <- NULL
  expect_error(
    iv_plot_spec(missing.coordinate, "npregiv", FALSE, TRUE, FALSE),
    "evaluation coordinate 'zeval' is missing",
    fixed = TRUE
  )

  mismatched <- iv
  mismatched$phi.eval <- c(1, 2, 3)
  expect_error(
    iv_plot_spec(mismatched, "npregiv", FALSE, TRUE, FALSE),
    "different lengths",
    fixed = TRUE
  )

  multivariate <- iv
  multivariate$z <- data.frame(z1 = 1:3, z2 = 4:6)
  expect_error(
    iv_plot_spec(multivariate, "npregiv", FALSE, TRUE, FALSE),
    "only univariate z is supported",
    fixed = TRUE
  )
})

test_that("IV renderer routes frame, overlay, curve, and rug arguments once", {
  iv <- iv_plot_object()
  spec <- iv_plot_spec(
    iv, "npregiv", FALSE, TRUE, TRUE,
    dots = list(type = "b", col = "red", lty = 2, lwd = 3,
                pch = 4, cex = 1.25, bg = "white", main = "IV")
  )
  calls <- new.env(parent = emptyenv())
  calls$order <- character()
  calls$plot <- calls$overlay <- calls$lines <- calls$rug <- NULL

  local_mocked_bindings(
    plot = function(...) {
      calls$order <- c(calls$order, "plot")
      calls$plot <- list(...)
      invisible(NULL)
    },
    lines = function(...) {
      calls$order <- c(calls$order, "lines")
      calls$lines <- list(...)
      invisible(NULL)
    },
    .np_plot_overlay_points_1d = function(...) {
      calls$order <- c(calls$order, "overlay")
      calls$overlay <- list(...)
      invisible(TRUE)
    },
    .np_plot_draw_rug_1d = function(...) {
      calls$order <- c(calls$order, "rug")
      calls$rug <- list(...)
      invisible(TRUE)
    },
    .package = iv_plot_package
  )

  value <- withVisible(
    getFromNamespace(".np_iv_plot_render", iv_plot_package)(spec)
  )
  expect_identical(value, list(value = NULL, visible = FALSE))
  expect_identical(calls$order, c("plot", "overlay", "lines", "rug"))
  expect_identical(calls$plot$type, "n")
  expect_identical(calls$plot$main, "IV")
  expect_null(calls$plot$col)
  expect_null(calls$plot$lty)
  expect_null(calls$plot$lwd)
  expect_identical(calls$overlay$col, "red")
  expect_identical(calls$overlay$pch, 4)
  expect_identical(calls$overlay$cex, 1.25)
  expect_identical(calls$overlay$bg, "white")
  expect_identical(calls$lines$type, "b")
  expect_identical(calls$lines$col, "red")
  expect_identical(calls$lines$lty, 2)
  expect_identical(calls$lines$lwd, 3)
  expect_null(calls$lines$pch)
  expect_identical(calls$rug$x, iv$z$z)
})

test_that("retired IV plot controls fail with canonical migration guidance", {
  iv <- iv_plot_object()
  deriv <- iv_deriv_plot_object()
  expect_error(plot(iv, plot.data = TRUE), "use data_overlay", fixed = TRUE)
  expect_error(plot(iv, deriv = TRUE), "use gradients", fixed = TRUE)
  expect_error(plot(iv, phi = 30), "not a supported IV plot argument",
               fixed = TRUE)
  expect_error(plot(deriv, phi = TRUE), "use gradients = FALSE", fixed = TRUE)
  expect_error(plot(deriv, plot.data = TRUE), "use data_overlay", fixed = TRUE)
  expect_error(plot(iv, 1, 2, 3, 4), "unnamed plot arguments", fixed = TRUE)
  expect_error(plot(iv, not_a_plot_argument = TRUE), "unused plot argument",
               fixed = TRUE)
  expect_error(plot(iv, col = 1, col = 2), "supplied more than once",
               fixed = TRUE)
})

test_that("all IV plot flags are validated before object fields or devices", {
  empty.iv <- structure(list(), class = "npregiv")
  expect_error(plot(empty.iv, gradients = NA), "'gradients' must be TRUE or FALSE",
               fixed = TRUE)
  expect_error(plot(empty.iv, data_overlay = NA),
               "'data_overlay' must be TRUE or FALSE", fixed = TRUE)
  expect_error(plot(empty.iv, data_rug = NA),
               "'data_rug' must be TRUE or FALSE", fixed = TRUE)
})

test_that("real devices show expanded IV geometry without state leakage", {
  iv <- iv_plot_object()
  deriv <- iv_deriv_plot_object()
  iv.before <- iv
  deriv.before <- deriv
  device.before <- grDevices::dev.cur()
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    if (grDevices::dev.cur() != device.before)
      grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  stable.par <- graphics::par(c("mfrow", "mar", "oma"))

  value <- withVisible(plot(iv, xaxs = "i", yaxs = "i", data_rug = TRUE))
  expect_identical(value, list(value = NULL, visible = FALSE))
  expect_equal(graphics::par("usr"), c(-2, 2, -10, 10), tolerance = 0)
  expect_identical(graphics::par(c("mfrow", "mar", "oma")), stable.par)

  derivative.value <- withVisible(plot(deriv, data_rug = TRUE))
  expect_identical(derivative.value, list(value = NULL, visible = FALSE))
  level.value <- withVisible(plot(deriv, gradients = FALSE,
                                  data_overlay = FALSE))
  expect_identical(level.value, list(value = NULL, visible = FALSE))
  expect_identical(iv, iv.before)
  expect_identical(deriv, deriv.before)

  grDevices::dev.off()
  expect_identical(grDevices::dev.cur(), device.before)
})

test_that("categorical IV displays retain base behavior and omit numerical rugs", {
  categorical <- structure(
    list(
      z = data.frame(z = factor(c("a", "b", "a"))),
      y = c(1, 2, 3),
      phi = c(1.2, 2.2, 2.8),
      phi.deriv.1 = matrix(0, nrow = 3L, ncol = 1L)
    ),
    class = "npregiv"
  )
  explicit.xlim <- c(0.5, 2.5)
  spec <- iv_plot_spec(
    categorical, "npregiv", FALSE, TRUE, TRUE,
    dots = list(xlim = explicit.xlim)
  )
  expect_false(spec$continuous)
  expect_false(spec$rug)
  expect_identical(spec$curve.args$xlim, explicit.xlim)

  calls <- new.env(parent = emptyenv())
  calls$plot <- NULL
  local_mocked_bindings(
    plot = function(...) {
      calls$plot <- list(...)
      invisible(NULL)
    },
    lines = function(...) invisible(NULL),
    .package = iv_plot_package
  )
  getFromNamespace(".np_iv_plot_render", iv_plot_package)(spec)
  expect_identical(calls$plot$xlim, explicit.xlim)

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  expect_silent(plot(categorical, data_rug = TRUE, xlim = explicit.xlim))
})
