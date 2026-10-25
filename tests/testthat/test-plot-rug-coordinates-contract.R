# Unique contract: floor rugs share surface-grid metric coordinates. Existing
# RGL ownership/HiDPI tests cover resources; these tiny renderer-boundary probes
# add no fitting, bandwidth search, bootstrap or MPI work to the default lane.
.rug_coordinate_package <- if (isNamespaceLoaded('npRmpi')) 'npRmpi' else 'np'

test_that('floor rugs map ordered levels to surface coordinates in both renderers', {
  skip_if_not(suppressWarnings(requireNamespace('rgl', quietly = TRUE)))
  base <- getFromNamespace('.np_plot_draw_floor_rug_persp', .rug_coordinate_package)
  rgl <- getFromNamespace('.np_plot_draw_floor_rug_rgl', .rug_coordinate_package)
  cases <- list(
    numeric = list(x = c(1951, 1970, NA, 1998, Inf), expected = c(1951, 1970, NA, 1998, Inf)),
    metric = list(x = ordered(c('1951', '1970', NA, '1998', '1970'),
                             levels = c('1951', '1960', '1970', '1998')),
                  expected = c(1951, 1970, NA, 1998, 1970)),
    labels = list(x = ordered(c('low', 'medium', NA, 'high', 'medium'),
                             levels = c('low', 'unused', 'medium', 'high')),
                  expected = c(1, 3, NA, 4, 3)))
  for (case in cases) for (axis in 1:2) {
    input <- case$x; saved <- input; other <- c(.1, .3, .5, .7, NA)
    x1 <- if (axis == 1L) input else other
    x2 <- if (axis == 2L) input else other
    a <- if (axis == 1L) case$expected else other
    b <- if (axis == 2L) case$expected else other
    ok <- is.finite(a) & is.finite(b); capture <- new.env(parent = emptyenv())
    testthat::with_mocked_bindings(
      expect_true(base(x1, x2, c(0, 2), diag(4))),
      segments = function(...) capture$base <- list(...), .package = 'graphics')
    testthat::with_mocked_bindings(
      expect_true(rgl(x1, x2, c(0, 2))),
      segments3d = function(...) capture$rgl <- list(...), .package = 'rgl')
    expect_identical(input, saved)
    expect_equal(unname(capture$base$x0), a[ok], tolerance = 0)
    expect_equal(unname(capture$base$y0), b[ok], tolerance = 0)
    expect_equal(unname(capture$base$x1), a[ok], tolerance = 0)
    expect_equal(unname(capture$base$y1), b[ok], tolerance = 0)
    expect_equal(capture$rgl$x, rbind(a[ok], a[ok]), tolerance = 0)
    expect_equal(capture$rgl$y, rbind(b[ok], b[ok]), tolerance = 0)
    expect_equal(capture$rgl$z, rbind(rep(0, sum(ok)), rep(.04, sum(ok))), tolerance = 0)
  }
})

test_that('floor rugs retain empty-input and caller style behavior', {
  skip_if_not(suppressWarnings(requireNamespace('rgl', quietly = TRUE)))
  base <- getFromNamespace('.np_plot_draw_floor_rug_persp', .rug_coordinate_package)
  rgl <- getFromNamespace('.np_plot_draw_floor_rug_rgl', .rug_coordinate_package)
  missing <- ordered(c(NA, NA), levels = c('1951','1998'))
  expect_false(base(missing, 1:2, c(0,1), diag(4)))
  expect_false(rgl(missing, 1:2, c(0,1)))
  capture <- new.env(parent = emptyenv())
  testthat::with_mocked_bindings(
    rgl(1:2, 3:4, c(2,2), list(color='red',lwd=3)),
    segments3d = function(...) capture$args <- list(...), .package='rgl')
  expect_identical(capture$args$color, 'red')
  expect_identical(capture$args$lwd, 3)
  expect_equal(capture$args$z, rbind(c(2,2),c(2.04,2.04)), tolerance=0)
})
