# This map is confined to the public plot boundary. Estimator and control-
# constructor spellings (including gradient.order and bar.num) do not change.
.np_plot_argument_map <- function() {
  c(bias.control = "bias_control", boot.control = "boot_control",
    grid.control = "grid_control", render.control = "render_control",
    gaussian.draws = "gaussian_draws", random.seed = "random_seed",
    gradient.order = "gradient_order", coef.index = "coef_index",
    data.overlay = "data_overlay", data.rug = "data_rug",
    common.scale = "common_scale", factor.boxplot = "factor_boxplot",
    boxplot.outliers = "boxplot_outliers", proper.method = "proper_method",
    proper.control = "proper_control", plot.vars = "plot_vars",
    fit.mean.train = "fit_mean_train")
}

.np_plot_public_names <- function(names) {
  map <- .np_plot_argument_map()
  hit <- match(names, names(map))
  names[!is.na(hit)] <- unname(map[hit[!is.na(hit)]])
  names
}

.np_plot_internal_names <- function(dots) {
  map <- .np_plot_argument_map()
  hit <- match(names(dots), unname(map))
  names(dots)[!is.na(hit)] <- names(map)[hit[!is.na(hit)]]
  if (anyDuplicated(names(dots)))
    stop("duplicate plot argument", call. = FALSE)
  dots
}

.np_plot_check_public_names <- function(names) {
  map <- .np_plot_argument_map()
  legacy <- c(gradient = "gradients", persp = "perspective", behavior = "output",
    plot.behavior = "output", plot.par.mfrow = "layout",
    plot.data.overlay = "data_overlay", plot.rug = "data_rug",
    plot.bxp = "factor_boxplot", plot.bxp.out = "boxplot_outliers",
    plot.errors.method = "errors", plot.errors.type = "band",
    plot.errors.alpha = "alpha", plot.errors.center = "center",
    plot.errors.boot.method = "bootstrap", plot.errors.boot.num = "B",
    plot.errors.boot.nonfixed = "boot_control = np_boot_control(nonfixed = ...)",
    plot.errors.boot.wild = "boot_control = np_boot_control(wild = ...)",
    plot.errors.boot.blocklen = "boot_control = np_boot_control(blocklen = ...)",
    plot.errors.boot.calibration = "boot_control = np_boot_control(calibration = ...)",
    plot.errors.style = "render_control = np_render_control(style = ...)",
    plot.errors.bar = "render_control = np_render_control(bar = ...)",
    plot.errors.bar.num = "render_control = np_render_control(bar.num = ...)")
  replacements <- c(map, legacy)
  bad <- intersect(names, names(replacements))
  if (length(bad))
    stop(paste0("retired plot argument: ", bad[1L], "; use ", replacements[[bad[1L]]]), call. = FALSE)
  if (anyDuplicated(names)) stop("duplicate plot argument", call. = FALSE)
  invisible(NULL)
}

# Request provenance exposes only public vocabulary, even after engine
# normalization. This separate map must not alter admission or dispatch.
.np_plot_data_supplied_names <- function(x) {
  map <- c(plot.errors.method="errors", plot.errors.type="band",
    plot.errors.alpha="alpha", plot.errors.center="center",
    plot.errors.boot.method="bootstrap", plot.errors.boot.num="B",
    plot.data.overlay="data_overlay", plot.rug="data_rug",
    plot.par.mfrow="layout", plot.bxp="factor_boxplot", plot.bxp.out="boxplot_outliers")
  i <- match(x, names(map)); x[!is.na(i)] <- unname(map[i[!is.na(i)]])
  x[startsWith(x, "plot.errors.boot.")] <- "boot_control"
  x[startsWith(x, "plot.errors.")] <- "render_control"
  unique(.np_plot_public_names(x[!x %in% c("plot.behavior", ".np.empty.report")]))
}
