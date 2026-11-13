# Public compound plot arguments use dotted names, as do their engines.
# Retired names are rejection metadata only; they are never translated.
.np_plot_check_public_names <- function(names, replay = FALSE) {
  retired <- c(bias_control = "bias.control", boot_control = "boot.control",
    grid_control = "grid.control", render_control = "render.control",
    gaussian_draws = "gaussian.draws", random_seed = "random.seed",
    gradient_order = "gradient.order", coef_index = "coef.index",
    data_overlay = "data.overlay", data_rug = "data.rug",
    common_scale = "common.scale", factor_boxplot = "factor.boxplot",
    boxplot_outliers = "boxplot.outliers", proper_method = "proper.method",
    proper_control = "proper.control", plot_vars = "plot.vars",
    fit_mean_train = "fit.mean.train")
  # Replay only consumes graphics overrides; do not change its incumbent
  # duplicate/legacy-graphics handling or advertise fresh inference controls.
  if (replay) {
    bad <- intersect(names, names(retired))
    if (length(bad))
      stop(paste0("retired plot argument: ", bad[1L],
        "; plot-data replay accepts graphical overrides, not estimation controls"),
        call. = FALSE)
    return(invisible(NULL))
  }
  legacy <- c(gradient = "gradients", persp = "perspective", behavior = "output",
    plot.behavior = "output", plot.par.mfrow = "layout",
    plot.data.overlay = "data.overlay", plot.rug = "data.rug",
    plot.bxp = "factor.boxplot", plot.bxp.out = "boxplot.outliers",
    plot.errors.method = "errors", plot.errors.type = "band",
    plot.errors.alpha = "alpha", plot.errors.center = "center",
    plot.errors.boot.method = "bootstrap", plot.errors.boot.num = "B",
    plot.errors.boot.nonfixed = "boot.control = np_boot_control(nonfixed = ...)",
    plot.errors.boot.wild = "boot.control = np_boot_control(wild = ...)",
    plot.errors.boot.blocklen = "boot.control = np_boot_control(blocklen = ...)",
    plot.errors.boot.calibration = "boot.control = np_boot_control(calibration = ...)",
    plot.errors.style = "render.control = np_render_control(style = ...)",
    plot.errors.bar = "render.control = np_render_control(bar = ...)",
    plot.errors.bar.num = "render.control = np_render_control(bar.num = ...)")
  replacements <- c(retired, legacy)
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
    plot.data.overlay="data.overlay", plot.rug="data.rug",
    plot.par.mfrow="layout", plot.bxp="factor.boxplot", plot.bxp.out="boxplot.outliers")
  i <- match(x, names(map)); x[!is.na(i)] <- unname(map[i[!is.na(i)]])
  x[startsWith(x, "plot.errors.boot.")] <- "boot.control"
  x[startsWith(x, "plot.errors.")] <- "render.control"
  unique(x[!x %in% c("plot.behavior", ".np.empty.report")])
}

.np_plot_unique_args <- function(dots) {
  if (anyDuplicated(names(dots)))
    stop("duplicate plot argument", call. = FALSE)
  dots
}
