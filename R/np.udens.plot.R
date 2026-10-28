# Family-specific controls; regression iterations and pilot controls do not apply.
np_udens_bias_control <- function(construction = c("joint-scale", "positive-covariance", "legacy")) {
  construction <- match.arg(construction)
  version <- if (identical(construction, "positive-covariance")) 1L else 3L
  structure(list(construction = construction, calibration = "absolute", version = version),
            class = "np_udens_bias_control")
}

.np_udens_plot_prepare <- function(control, bws, x, center, method, bootstrap, band) {
  if (is.null(control)) {
    if (!identical(center, "bias-corrected")) return(NULL)
    control <- np_udens_bias_control()
  }
  version <- if (is.list(control)) control[["version", exact = TRUE]] else NULL
  if (!inherits(control, "np_udens_bias_control") || !is.list(control) ||
      !is.character(control[["construction", exact = TRUE]]) ||
      length(control[["construction", exact = TRUE]]) != 1L ||
      !control[["construction", exact = TRUE]] %in% c("joint-scale", "positive-covariance", "legacy") ||
      length(version) != 1L || !is.integer(version) ||
      is.na(version) || !version %in% c(1L, 2L, 3L))
    stop("bias.control must be created by np_udens_bias_control()", call. = FALSE)
  if (identical(control[["construction", exact = TRUE]], "joint-scale") && version != 3L)
    stop("joint-scale density controls require version 3", call. = FALSE)
  expected <- np_udens_bias_control(control[["construction", exact = TRUE]])
  expected$version <- version
  if (!identical(control, expected) || !is.integer(version))
    stop("bias.control must be created by np_udens_bias_control()", call. = FALSE)
  if (!identical(center, "bias-corrected") || !method %in% c("asymptotic", "bootstrap"))
    stop("density bias.control requires center='bias-corrected' and asymptotic or bootstrap errors", call. = FALSE)
  if (identical(control[["construction", exact = TRUE]], "legacy")) {
    if (!identical(method, "bootstrap"))
      stop("legacy density bias correction requires errors='bootstrap'", call. = FALSE)
    return(NULL)
  }
  if (version == 1L && !identical(method, "bootstrap"))
    stop("saved version-1 density controls require errors='bootstrap'", call. = FALSE)
  if (identical(method, "bootstrap") && !identical(bootstrap, "inid"))
    stop("positive-covariance density bounds require bootstrap='inid'; other laws remain unqualified", call. = FALSE)
  if (!band %in% c("pmzsd", "pointwise", "simultaneous", "bonferroni", "all"))
    stop("unsupported bias-corrected density band", call. = FALSE)
  state <- .np_udens_bc_prepare(bws, x)
  state$plot.control.version <- version
  state$plot.law <- if (identical(method, "asymptotic")) "gaussian" else "inid"
  state
}

.np_udens_plot_errors <- function(state, grid, alpha, B, band) {
  names(grid) <- names(state$x)
  provider.band <- if (band == "pmzsd") "pointwise" else band
  inference <- .np_plot_activity_run("Constructing bias-corrected density intervals",
    if (identical(state$plot.control.version, 3L))
      .np_udens_scale_inference(state, grid, law = state$plot.law,
        alpha = alpha, B = B,
        bands = if (provider.band == "all") c("pointwise", "simultaneous", "bonferroni") else provider.band)
    else if (identical(state$plot.control.version, 2L))
      .np_udens_reference_inference(state, grid, law = state$plot.law,
        alpha = alpha, B = B,
        bands = if (provider.band == "all") c("pointwise", "simultaneous", "bonferroni") else provider.band)
    else .np_udens_b4_inference(state, grid, law = "inid", alpha = alpha, B = B))
  inference$settings$primary <- "absolute"
  inference$bands <- inference$absolute.intervals
  if (band == "pmzsd") inference <- .np_plot_bias_normal(inference, alpha)
  all <- lapply(inference$bands, function(z)
    cbind(inference$center - z[, 1L, 1L], z[, 2L, 1L] - inference$center))
  width <- all[[if (band == "all") "pointwise" else band]]
  list(boot.err = cbind(width, inference$center),
       boot.all.err = if (band == "all") all else NULL,
       inference = inference, bxp = list())
}

# Density-specific payload: the Gaussian provider has already calibrated its
# complete corrected process. Do not pass it through legacy SE calibration.
.np_udens_plot_payload <- function(raw) {
  list(err = raw$boot.err, all.err = raw$boot.all.err,
       center = raw$inference$center, bxp = raw$bxp)
}
