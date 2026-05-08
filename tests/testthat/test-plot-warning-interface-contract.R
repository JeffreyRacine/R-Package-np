test_that("bootstrap coverage warning uses current public plot option names", {
  quantile_bounds <- getFromNamespace("compute.bootstrap.quantile.bounds", "np")

  set.seed(20260507L)
  boot.t <- matrix(rnorm(25L * 3L), nrow = 25L, ncol = 3L)
  msg <- NULL

  withCallingHandlers(
    quantile_bounds(
      boot.t = boot.t,
      alpha = 0.05,
      band.type = "all",
      warn.coverage = TRUE
    ),
    warning = function(w) {
      msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  expect_type(msg, "character")
  expect_match(msg, "B=25", fixed = TRUE)
  expect_match(msg, "band=\"bonferroni/simultaneous/all\"", fixed = TRUE)
  expect_match(msg, "B >= ceiling(2*m/alpha - 1)", fixed = TRUE)
  expect_false(grepl("plot.errors.boot.num", msg, fixed = TRUE))
  expect_false(grepl("plot.errors.type", msg, fixed = TRUE))
})

test_that("common plot option validation reports current public names", {
  normalize <- getFromNamespace(".np_plot_normalize_common_options", "np")

  expect_error(
    normalize(
      plot.behavior = "data",
      plot.errors.method = "none",
      plot.errors.boot.method = "wild",
      plot.errors.boot.wild = "rademacher",
      plot.errors.boot.blocklen = NULL,
      plot.errors.center = "estimate",
      plot.errors.type = "pointwise",
      plot.errors.alpha = 0.5,
      plot.errors.style = "band",
      plot.errors.bar = "|",
      xdat = data.frame(x = 1:3),
      common.scale = FALSE,
      ylim = NULL
    ),
    "alpha must lie in \\(0, 0\\.5\\)"
  )

  msg <- NULL
  out <- withCallingHandlers(
    normalize(
      plot.behavior = "data",
      plot.errors.method = "none",
      plot.errors.boot.method = "wild",
      plot.errors.boot.wild = "rademacher",
      plot.errors.boot.blocklen = NULL,
      plot.errors.center = "estimate",
      plot.errors.type = "all",
      plot.errors.alpha = 0.05,
      plot.errors.style = "band",
      plot.errors.bar = "|",
      xdat = data.frame(x = 1:3),
      common.scale = FALSE,
      ylim = NULL
    ),
    warning = function(w) {
      msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(out$plot.errors.method, "bootstrap")
  expect_match(msg, "band=\"all\" requires bootstrap errors", fixed = TRUE)
  expect_match(msg, "errors=\"bootstrap\"", fixed = TRUE)
  expect_false(grepl("plot.errors.type", msg, fixed = TRUE))
  expect_false(grepl("plot.errors.method", msg, fixed = TRUE))
})

test_that("snake_case plot controls normalize to engine controls", {
  normalize <- getFromNamespace(".np_plot_normalize_public_dots", "np")

  dots <- normalize(
    list(
      output = "both",
      data_overlay = FALSE,
      data_rug = TRUE,
      layout = "current",
      factor_boxplot = TRUE,
      boxplot_outliers = FALSE,
      coef_index = 2L,
      gradient_order = c(1L, 2L),
      common_scale = FALSE,
      proper_method = "isotonic",
      proper_control = list(mode = "slice"),
      errors = "asymptotic",
      render_control = np_render_control(style = "bar", bar = "I", bar_num = 7L)
    )
  )

  expect_identical(dots$plot.errors.method, "asymptotic")
  expect_identical(dots$plot.behavior, "plot-data")
  expect_false(dots$plot.data.overlay)
  expect_true(dots$plot.rug)
  expect_false(dots$plot.par.mfrow)
  expect_true(dots$plot.bxp)
  expect_false(dots$plot.bxp.out)
  expect_identical(dots$coef.index, 2L)
  expect_identical(dots$gradient.order, c(1L, 2L))
  expect_false(dots$common.scale)
  expect_identical(dots$proper.method, "isotonic")
  expect_identical(dots$proper.control, list(mode = "slice"))
  expect_identical(dots$plot.errors.style, "bar")
  expect_identical(dots$plot.errors.bar, "I")
  expect_identical(dots$plot.errors.bar.num, 7L)
})

test_that("snake_case plot controls fail cleanly on conflicts", {
  normalize <- getFromNamespace(".np_plot_normalize_public_dots", "np")

  expect_error(
    normalize(list(output = "data", plot.behavior = "plot")),
    "conflicting plot arguments: output and plot.behavior"
  )
  expect_error(
    normalize(list(data_rug = TRUE, plot.rug = FALSE)),
    "conflicting plot arguments: data_rug and plot.rug"
  )
  expect_error(
    normalize(list(layout = "maybe")),
    "layout must be one of"
  )
  expect_error(
    np_render_control(bar.num = 3L),
    "unused argument"
  )
})
