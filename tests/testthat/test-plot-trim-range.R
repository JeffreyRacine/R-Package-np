test_that("public trim ranges apply to every numeric X margin", {
  d <- data.frame(x = seq(-1, 1, length.out = 30), z = cos(seq_len(30)))
  b <- npudensbw(~ x + z, data = d, bws = c(.3, .3), bandwidth.compute = FALSE)
  for (surface in c(FALSE, TRUE)) {
    out <- plot(b, output = "data", neval = 5L, perspective = surface,
                grid_control = np_grid_control(xtrim = c(.1, .8)))
    for (p in out$panels) for (j in seq_along(p$grid)) {
      value <- p$grid[[j]]
      role <- if (surface) b$xnames[j] else p$variables$name[j]
      if (length(unique(value)) <= 1L) next
      expected <- seq(quantile(d[[role]], .1), quantile(d[[role]], .8), length.out = 5L)
      expect_equal(sort(unique(value)), expected, tolerance = 1e-12)
    }
  }
})

test_that("trim adapters preserve legacy and full-range percentile families", {
  wrap <- function(x) structure(list(probs = x), class = "np_plot_trim_range")
  d <- data.frame(x = seq(-1, 1, length.out = 30))
  expect_identical(.np_plot_recycle_xtrim(.1, 3L), rep(.1, 3L))
  expect_identical(trim.quantiles(d$x, .1), quantile(d$x, c(.1, .9)))
  expect_identical(trim.quantiles(d$x, .np_plot_recycle_xtrim(wrap(c(.1, .8)), 1L)),
                   quantile(d$x, c(.1, .8)))
  grid <- data.frame(x = seq(-1, 1, length.out = 5))
  state <- list(x = d)
  design <- list(vary = 1L, xtrim = .np_plot_recycle_xtrim(wrap(c(0, 1)), 1L), xq = .5)
  expect_null(.np_udens_calibrated_family(state, grid, design)$builder)
  design$xtrim <- .np_plot_recycle_xtrim(wrap(c(.1, .8)), 1L)
  builder <- .np_udens_calibrated_family(state, grid, design)$builder
  outer <- list(x = data.frame(x = d$x * 2 + .2), category = list())
  expected <- seq(quantile(outer$x$x, .1), quantile(outer$x$x, .8), length.out = 5)
  expect_equal(builder(outer)$x, expected, tolerance = 1e-12)
})
