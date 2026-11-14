test_that("automatic RGL clearance preserves user limits and unusual ranges", {
  automatic <- .np_plot_rgl_automatic_limits
  pad <- .np_plot_rgl_pad_automatic_limits
  expect_identical(automatic(list()), c("xlim", "ylim", "zlim"))
  expect_identical(automatic(list(xlim = c(0, 1), zlim = c(0, NA))), "ylim")
  expect_identical(automatic(list(xlim = NULL)), c("xlim", "ylim", "zlim"))
  original <- list(xlim = c(-2, 2), ylim = c(0, 1), zlim = c(-1, 3))
  expect_identical(pad(original, character(), list()), original)
  expect_identical(pad(original, names(original), original), original)
  expect_identical(pad(original, names(original), list(ylim = NULL)),
    list(xlim = c(-2.004, 2.004), ylim = c(0, 1), zlim = c(-1.004, 3.004)))
  for (value in list(NULL, c(0, NA_real_), c(1, 1), c(2, 0), c(0, Inf),
                     c(-1e308, 1e308))) {
    a <- list(zlim = value)
    expect_identical(pad(a, "zlim", list()), a)
  }
  for (scale in c(1e-100, 1, 1e100)) {
    got <- pad(list(zlim = c(-scale, scale)), "zlim", list())$zlim
    expect_equal(got / scale, c(-1.002, 1.002), tolerance = 1e-14)
  }
})

test_that("RGL clearance changes clipping planes without changing geometry or RNG", {
  skip_if_not_installed("rgl")
  withr::local_options(rgl.useNULL = TRUE)
  withr::local_envvar(RGL_USE_NULL = "TRUE")
  render <- .np_plot_render_surface_rgl
  pkg <- getNamespaceName(environment(render))
  seen <- new.env(parent = emptyenv())
  original.persp <- rgl::persp3d
  testthat::local_mocked_bindings(persp3d = function(x, y, z, ...) {
    seen$args <- list(...)
    seen$geometry <- list(x = x, y = y, z = z)
    ans <- original.persp(x, y, z, ...)
    seen$bbox <- as.numeric(rgl::par3d("bbox"))
    ans
  }, rglwidget = function(x, ...) list(scene = x), .package = "rgl")
  testthat::local_mocked_bindings(.np_plot_rgl_hidpi_widget = function(widget, ...) widget,
                                .package = pkg)
  axis <- c(0, .5, 1)
  surfaces <- list(floor = outer(axis, axis, "*"),
                   ceiling = 1 - outer(axis, axis, "*"),
                   negative = outer(axis, axis, "-") - 1,
                   constant = matrix(.5, 3L, 3L))
  before.devices <- as.integer(rgl::rgl.dev.list())
  set.seed(184L); before.seed <- .Random.seed
  for (z in surfaces) {
    a <- list(x = axis, y = axis, z = z, xlab = "x", ylab = "y", zlab = "z",
              main = "", theta = 0, phi = 20,
              xlim = c(0, 1), ylim = c(0, 1), zlim = c(-2, 2))
    capture.output(do.call(render, c(a, list(automatic.limits = c("xlim", "ylim", "zlim")))))
    expect_equal(seen$bbox, c(-.001, 1.001, -.001, 1.001, -2.004, 2.004), tolerance = 1e-6)
    expect_identical(seen$geometry, a[c("x", "y", "z")])
    capture.output(do.call(render, c(a, list(automatic.limits = "zlim"))))
    expect_identical(seen$args$xlim, c(0, 1))
    expect_identical(seen$args$ylim, c(0, 1))
    capture.output(do.call(render, c(a, list(automatic.limits = "zlim",
      persp3d.args = list(zlim = c(0, 1))))))
    expect_identical(seen$args$zlim, c(0, 1))
    expect_identical(.Random.seed, before.seed)
    expect_identical(as.integer(rgl::rgl.dev.list()), before.devices)
  }
})

test_that("fresh and replay surfaces retain automatic-limit provenance", {
  skip_if_not_installed("rgl")
  withr::local_options(np.messages = FALSE, rgl.useNULL = TRUE)
  withr::local_envvar(RGL_USE_NULL = "TRUE")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  seen <- list()
  pkg <- getNamespaceName(environment(.np_plot_render_surface_rgl))
  # Capture the actual shared-renderer entry; geometry is exercised above.
  capture <- function(...) {
    seen[[length(seen) + 1L]] <<- as.list(environment())
    invisible(NULL)
  }
  # Keep the real matching contract, including inherited abbreviated arguments.
  formals(capture) <- formals(.np_plot_render_surface_rgl)
  testthat::local_mocked_bindings(.np_plot_render_surface_rgl = capture, .package = pkg)
  fits <- evalq({
    d <- data.frame(x = seq(.05, .95, length.out = 30L),
                    w = cos(seq(.1, 5, length.out = 30L)))
    d$y <- sin(3 * d$x) + .2 * d$w
    d$g <- factor(rep(c("a", "b"), 15L))
    list(density = npudens(~x + w, data = d, bws = c(.3, .4), bandwidth.compute = FALSE),
      distribution = npudist(~x + w, data = d, bws = c(.3, .4), bandwidth.compute = FALSE),
      cdensity = npcdens(y ~ x, data = d, bws = c(.3, .4), bandwidth.compute = FALSE),
      cdistribution = npcdist(y ~ x, data = d, bws = c(.3, .4), bandwidth.compute = FALSE),
      regression = npreg(y ~ x + w, data = d, bws = c(.3, .4), bandwidth.compute = FALSE),
      plreg = npplreg(y ~ x | w, data = d, bws = matrix(c(.4, .4), 2L, 1L), bandwidth.compute = FALSE),
      scoef = npscoef(y ~ x | w, data = d, bws = .4, bandwidth.compute = FALSE),
      copula = npcopula(~x + w, data = d, bws = c(.3, .4), neval = 4L),
      conmode = npconmode(g ~ x + w, data = d, bws = c(.1, .3, .4), probabilities = TRUE))
  }, new.env(parent = globalenv()))
  for (family in names(fits)) {
    fit <- fits[[family]]
    args <- list(x = fit, neval = 4L, view = "fixed")
    if (family == "copula") args$neval <- NULL
    if (family == "conmode") {
      args$level <- "a"
      args$perspective <- TRUE
    }
    fitted.before <- serialize(fit, NULL)
    set.seed(184L); before.seed <- .Random.seed
    baseline <- do.call(plot, c(args, list(output = "data")))
    seen <- list()
    drawn <- do.call(plot, c(args, list(renderer = "rgl", output = "plot-data")))
    expect_length(seen, 1L)
    expect_setequal(seen[[1L]]$automatic.limits, c("xlim", "ylim", "zlim"))
    expect_identical(baseline$panels[[1L]][c("grid", "values", "intervals")],
                     drawn$panels[[1L]][c("grid", "values", "intervals")])
    saved <- serialize(drawn, NULL)
    seen <- list(); plot(drawn, renderer = "rgl", zlim = c(0, .2))
    expect_length(seen, 1L)
    expect_setequal(seen[[1L]]$automatic.limits, c("xlim", "ylim"))
    expect_identical(seen[[1L]]$zlim, c(0, .2))
    expect_identical(serialize(drawn, NULL), saved)
    for (axis in c("x", "y")) {
      seen <- list()
      shortened <- setNames(list(c(0, 1)), paste0(axis, "li"))
      do.call(plot, c(list(x = drawn, renderer = "rgl"), shortened))
      expect_identical(seen[[1L]][[paste0(axis, "lim")]], c(0, 1))
      expect_false(paste0(axis, "lim") %in% seen[[1L]]$automatic.limits)
      expect_identical(serialize(drawn, NULL), saved)
    }
    seen <- list()
    explicit <- do.call(plot, c(args, list(renderer = "rgl",
      output = "plot-data", xlim = c(-1, 2), zlim = c(0, .2))))
    expect_setequal(seen[[1L]]$automatic.limits, "ylim")
    seen <- list(); plot(explicit, renderer = "rgl")
    expect_setequal(seen[[1L]]$automatic.limits, "ylim")
    expect_identical(seen[[1L]]$xlim, c(-1, 2))
    expect_identical(seen[[1L]]$zlim, c(0, .2))
    if (family == "copula") {
      # Preserve limits already accepted by copula's inherited list-prefix lookup.
      for (axis in c("x", "y")) {
        seen <- list()
        alternate <- setNames(list(c(0, 1)), paste0(axis, "limit"))
        do.call(plot, c(args, list(renderer = "rgl", output = "plot-data"), alternate))
        expect_identical(seen[[1L]][[paste0(axis, "lim")]], c(0, 1))
        expect_false(paste0(axis, "lim") %in% seen[[1L]]$automatic.limits)
      }
    }
    expect_identical(serialize(fit, NULL), fitted.before)
    expect_identical(.Random.seed, before.seed)
  }
})
