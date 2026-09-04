.plot_annotation_package <- if (isNamespaceLoaded("npRmpi")) "npRmpi" else "np"

test_that("variability descriptors distinguish bootstrap and asymptotic constructions", {
  descriptor <- getFromNamespace(".np_plot_variability_descriptor",
                                 .plot_annotation_package)

  bootstrap.pmzsd <- descriptor("bootstrap", "pmzsd", 0.05, "estimate")
  expect_identical(bootstrap.pmzsd$level, "95%")
  expect_match(bootstrap.pmzsd$candidates[[1L]], "bootstrap SD", fixed = TRUE)
  expect_match(bootstrap.pmzsd$candidates[[1L]], "symmetric", fixed = TRUE)

  bootstrap.pointwise <- descriptor("bootstrap", "pointwise", 0.10, "estimate")
  expect_identical(bootstrap.pointwise$level, "90%")
  expect_match(bootstrap.pointwise$candidates[[1L]], "quantiles", fixed = TRUE)
  expect_false(any(grepl("symmetric", bootstrap.pointwise$candidates, fixed = TRUE)))

  bootstrap.simultaneous <- descriptor("bootstrap", "simultaneous", 0.025,
                                       "estimate")
  expect_identical(bootstrap.simultaneous$level, "97.5%")
  expect_match(bootstrap.simultaneous$candidates[[1L]], "rank-based", fixed = TRUE)

  asym.pmzsd <- descriptor("asymptotic", "pmzsd", 0.05, "estimate")
  asym.pointwise <- descriptor("asymptotic", "pointwise", 0.05, "estimate")
  expect_identical(asym.pmzsd$candidates, asym.pointwise$candidates)
  expect_match(asym.pmzsd$candidates[[1L]], "asymptotic SE", fixed = TRUE)

  asym.simultaneous <- descriptor("asymptotic", "simultaneous", 0.05,
                                  "estimate")
  expect_length(asym.simultaneous$candidates, 0L)
})

test_that("variability level formatting is compact and numerically stable", {
  format_level <- getFromNamespace(".np_plot_variability_level",
                                   .plot_annotation_package)
  expect_identical(format_level(0.05), "95%")
  expect_identical(format_level(0.10), "90%")
  expect_identical(format_level(0.025), "97.5%")
  expect_identical(format_level(0.033333333333), "96.66667%")
  expect_false(grepl("e", format_level(0.123456789), ignore.case = TRUE))
})

test_that("bias-corrected bootstrap descriptors identify centering without changing method", {
  descriptor <- getFromNamespace(".np_plot_variability_descriptor",
                                 .plot_annotation_package)
  pointwise <- descriptor("bootstrap", "pointwise", 0.05, "bias-corrected")
  expect_true(all(grepl("bias-corrected", pointwise$candidates, fixed = TRUE)))
  expect_match(pointwise$candidates[[1L]], "quantiles", fixed = TRUE)

  all.bands <- descriptor("bootstrap", "all", 0.05, "bias-corrected")
  expect_length(all.bands$candidates, 2L)
  expect_true(all(grepl("bias-corrected center", all.bands$candidates,
                        fixed = TRUE)))
})

test_that("label selection is ordered and suppresses exact no-fit cells", {
  select_label <- getFromNamespace(".np_plot_variability_label_select",
                                   .plot_annotation_package)
  candidates <- c("preferred", "compact", "minimal")
  widths <- c(6, 4, 2)

  expect_identical(select_label(candidates, widths, 6), "preferred")
  expect_identical(select_label(candidates, widths, 4), "compact")
  expect_identical(select_label(candidates, widths, 2), "minimal")
  expect_null(select_label(candidates, widths, 1.999999))
  expect_null(select_label(candidates, c(NA, Inf, 3), 2))
})

test_that("annotation specs preserve missing-versus-supplied subtitle state", {
  make_spec <- getFromNamespace(".np_plot_variability_annotation_spec",
                                .plot_annotation_package)
  args <- list(xlab = "x", cex.sub = 0.9)

  expect_null(make_spec("none", "pmzsd", 0.05, "estimate", FALSE, args))
  expect_null(make_spec("bootstrap", "pointwise", 0.05, "estimate", TRUE, args))
  expect_null(make_spec("bootstrap", "pointwise", 0.05, "estimate", FALSE,
                        args, eligible = FALSE))
  spec <- make_spec("bootstrap", "pointwise", 0.05, "estimate", FALSE, args)
  expect_identical(spec$plot.args, args)
  expect_identical(spec$plot.errors.method, "bootstrap")
})

test_that("finite-band eligibility matches paired drawable coordinates", {
  finite_pair <- getFromNamespace(".np_plot_finite_band_pair",
                                  .plot_annotation_package)
  expect_true(finite_pair(c(NA, 1), c(NA, 2)))
  expect_false(finite_pair(c(NA, Inf), c(NA, 2)))
  expect_false(finite_pair(c(1, 2), c(NA)))
  expect_false(finite_pair(NULL, 1))
})

test_that("all-band surface ranges ignore unavailable families without finite drift", {
  surface_range <- getFromNamespace(".np_plot_all_surface_range",
                                    .plot_annotation_package)
  lower <- c(1, 2)
  upper <- c(4, 5)
  lower.all <- list(pointwise = c(0, 1),
                    simultaneous = c(NA_real_, NA_real_),
                    bonferroni = c(-1, 0))
  upper.all <- list(pointwise = c(5, 6),
                    simultaneous = c(NA_real_, NA_real_),
                    bonferroni = c(6, 7))
  expect_identical(
    surface_range(c(2, 3), lower, upper, lower.all, upper.all),
    c(-1, 7)
  )

  finite.lower <- lapply(lower.all, function(value) replace(value,
                                                            !is.finite(value), 0))
  finite.upper <- lapply(upper.all, function(value) replace(value,
                                                            !is.finite(value), 8))
  expect_identical(
    surface_range(c(2, 3), lower, upper, finite.lower, finite.upper),
    c(min(c(unlist(finite.lower), lower)),
      max(c(unlist(finite.upper), upper)))
  )
})

test_that("default all-band legends identify construction and omit absent bands", {
  labels <- getFromNamespace(".np_plot_all_band_legend_labels",
                             .plot_annotation_package)
  legend_args <- getFromNamespace(".np_plot_all_band_legend_args",
                                  .plot_annotation_package)

  expect_identical(
    unname(labels("bootstrap")),
    c("Pointwise (quantiles)", "Simultaneous (rank-based)",
      "Bonferroni (quantiles)")
  )
  asym <- legend_args(
    legend.value = TRUE,
    x = "topright",
    lty = 2,
    lwd = 1,
    bty = "n",
    plot.errors.method = "asymptotic",
    drawn = c(pointwise = TRUE, simultaneous = FALSE, bonferroni = TRUE)
  )
  expect_identical(
    asym$legend,
    c("Pointwise (z × SE)",
      "Bonferroni (adjusted z × SE)")
  )
  expect_length(asym$col, 2L)
})

test_that("explicit all-band legend labels remain unmodified and unfiltered", {
  legend_args <- getFromNamespace(".np_plot_all_band_legend_args",
                                  .plot_annotation_package)
  custom <- c("First", "Second", "Third")
  args <- legend_args(
    legend.value = list(legend = custom),
    x = "topright",
    lty = 2,
    lwd = 1,
    bty = "n",
    plot.errors.method = "asymptotic",
    drawn = c(pointwise = TRUE, simultaneous = FALSE, bonferroni = TRUE)
  )
  expect_identical(args$legend, custom)
  expect_length(args$col, 3L)
})

test_that("base annotation uses existing subtitle space without changing geometry", {
  draw <- getFromNamespace(".np_plot_draw_variability_annotation",
                           .plot_annotation_package)
  make_spec <- getFromNamespace(".np_plot_variability_annotation_spec",
                                .plot_annotation_package)

  grDevices::pdf(file = tempfile(fileext = ".pdf"), width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(1:5, 1:5, type = "l", xlab = "x", ylab = "estimate")
  before <- graphics::par(c("mar", "oma", "fig", "fin", "pin", "plt", "usr"))
  spec <- make_spec(
    "bootstrap", "pointwise", 0.05, "estimate", FALSE,
    list(xlab = "x", cex.sub = graphics::par("cex.sub"),
         font.sub = graphics::par("font.sub"), family = graphics::par("family"))
  )
  expect_true(draw(spec))
  after <- graphics::par(c("mar", "oma", "fig", "fin", "pin", "plt", "usr"))
  expect_identical(after, before)
})

test_that("annotation placement follows actual geometry and figure ownership", {
  placement <- getFromNamespace(".np_plot_variability_placement", .plot_annotation_package)
  context <- getFromNamespace(".np_plot_variability_context", .plot_annotation_package)
  begin <- getFromNamespace(".np_plot_variability_panel_begin", .plot_annotation_package)
  grDevices::pdf(NULL, width = 9, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(1:3)
  expect_identical(placement(context(FALSE)), "subtitle")
  expect_identical(placement(context(TRUE)), "subtitle")
  graphics::par(mfrow = c(2, 2))
  caller <- context(FALSE)
  begin(caller)
  expect_identical(placement(caller), "none")
  owned <- context(TRUE)
  expect_identical(placement(owned), "none")
  begin(owned)
  expect_identical(placement(owned), "shared")
  graphics::par(mfrow = c(1, 1), fig = c(0, 0.5, 0, 1))
  expect_identical(placement(caller), "none")
})

test_that("shared annotations draw once per page including missing-first panels", {
  context <- getFromNamespace(".np_plot_variability_context", .plot_annotation_package)(TRUE)
  begin <- getFromNamespace(".np_plot_variability_panel_begin", .plot_annotation_package)
  spec <- getFromNamespace(".np_plot_variability_annotation_spec", .plot_annotation_package)
  draw <- getFromNamespace(".np_plot_draw_variability_annotation", .plot_annotation_package)
  grDevices::pdf(NULL, width = 9, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2))
  # No interval at each page's first panel; draw once at the second instead.
  for (i in seq_len(10L)) {
    begin(context)
    graphics::plot(1:3, xlab = "x", ylab = "value")
    if (i %% 4L == 1L) next
    before <- graphics::par(c("mar", "oma", "fig", "fin", "pin", "plt", "usr", "mfg", "cex"))
    annotation <- spec("bootstrap", "simultaneous", 0.05, "estimate",
                        FALSE, list(xlab = "x"), context = context)
    expect_identical(draw(annotation), i %% 4L == 2L)
    expect_identical(graphics::par(names(before)), before)
    expect_false(draw(annotation))
  }
})

test_that("individual pages retain subtitles including categorical intervals", {
  context <- getFromNamespace(".np_plot_variability_context", .plot_annotation_package)(FALSE)
  spec <- getFromNamespace(".np_plot_variability_annotation_spec", .plot_annotation_package)
  draw <- getFromNamespace(".np_plot_draw_variability_annotation", .plot_annotation_package)
  descriptor <- getFromNamespace(".np_plot_variability_descriptor", .plot_annotation_package)
  grDevices::pdf(NULL, width = 9, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (continuous in c(TRUE, FALSE, TRUE)) {
    graphics::plot(1:3, xlab = "x")
    expect_true(draw(spec("bootstrap", "pointwise", 0.05, "estimate",
                          FALSE, list(xlab = "x"), context = context,
                          continuous = continuous)))
  }
  categorical <- descriptor("bootstrap", "pointwise", 0.05,
                             continuous = FALSE)
  expect_match(categorical$candidates[[1L]], "variability interval", fixed = TRUE)
  shared <- descriptor("bootstrap", "simultaneous", 0.05,
                        "bias-corrected", shared = TRUE)
  expect_true(all(grepl("within each panel", shared$candidates, fixed = TRUE)))
  expect_true(all(grepl("rank-based", shared$candidates, fixed = TRUE)))
  expect_true(all(grepl("bias-corrected", shared$candidates, fixed = TRUE)))
})

test_that("shared annotations respect title clearance and no-fit typography", {
  context <- getFromNamespace(".np_plot_variability_context", .plot_annotation_package)(TRUE)
  begin <- getFromNamespace(".np_plot_variability_panel_begin", .plot_annotation_package)
  spec <- getFromNamespace(".np_plot_variability_annotation_spec", .plot_annotation_package)
  draw <- getFromNamespace(".np_plot_draw_variability_annotation", .plot_annotation_package)
  grDevices::pdf(NULL, width = 9, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2))
  begin(context)
  graphics::plot(1:3, main = "User title")
  make <- function(args) spec("bootstrap", "pointwise", 0.05, "estimate",
                               FALSE, args, context = context)
  expect_false(draw(make(list(main = expression(beta)))))
  expect_false(draw(make(list(main = "User title", cex.main = 12))))
  expect_false(draw(make(list(cex.sub = 50))))
  expect_false(draw(make(list(main = "Positioned", line = 1))))
  expect_true(draw(make(list(main = "User title"))))
  graphics::par(mfrow = c(2, 2), mar = c(5.1, 4.1, 0.2, 2.1))
  begin(context)
  graphics::plot(1:3)
  expect_false(draw(make(list())))
})

test_that("multi-quantile annotations consume one finite draw and state their scope", {
  descriptor <- getFromNamespace(".np_plot_variability_descriptor", .plot_annotation_package)
  labels <- descriptor("bootstrap", "simultaneous", 0.05,
                       shared = TRUE, per.quantile = TRUE)$candidates
  expect_true(all(grepl("within each quantile curve", labels, fixed = TRUE)))
  expect_false(any(grepl("within each panel", labels, fixed = TRUE)))
  helper <- getFromNamespace(".np_plot_draw_multi_tau_errors", .plot_annotation_package)
  # Only intercept the final drawing seam; both actual finite-bound gates
  # and the multi-tau loop remain the production implementations.
  env <- new.env(parent = environment(helper))
  draw.errors <- getFromNamespace("draw.errors", .plot_annotation_package)
  environment(draw.errors) <- env
  env$draw.errors <- draw.errors
  env$.np_plot_draw_variability_annotation <- function(annotation) {
    if (is.null(annotation)) return(invisible(FALSE))
    env$count <- env$count + 1L
    invisible(TRUE)
  }
  environment(helper) <- env
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(1:3, type = "n")
  value <- matrix(1, 3, 2)
  err <- array(0.1, c(3, 3, 2))
  for (band in c("pointwise", "all")) {
    for (finite in c(FALSE, TRUE)) {
      env$count <- 0L
      err[,,1] <- NA_real_
      err[,,2] <- if (finite) 0.1 else NA_real_
      all.err <- lapply(1:2, function(j) list(pointwise = err[,1:2,j],
                                             bonferroni = err[,1:2,j]))
      helper(1:3, value, err, all.err, FALSE, TRUE, c(1,2),
             band, "band", "|", 3, annotation = list())
      expect_identical(env$count, if (finite) 1L else 0L)
    }
  }
})
