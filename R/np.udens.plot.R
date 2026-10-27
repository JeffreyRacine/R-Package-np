# Family-specific controls; regression iterations and pilot controls do not apply.
np_udens_bias_control <- function(construction = c("positive-covariance", "legacy")) {
  construction <- match.arg(construction)
  structure(list(construction = construction, calibration = "absolute", version = 1L),
            class = "np_udens_bias_control")
}

.np_udens_plot_prepare <- function(control, bws, x, center, method, bootstrap, band) {
  if (is.null(control)) return(NULL)
  if (!inherits(control, "np_udens_bias_control") || !is.list(control) ||
      !is.character(control[["construction", exact = TRUE]]) ||
      length(control[["construction", exact = TRUE]]) != 1L ||
      !control[["construction", exact = TRUE]] %in% c("positive-covariance", "legacy") ||
      !identical(control, np_udens_bias_control(control[["construction", exact = TRUE]])))
    stop("bias.control must be created by np_udens_bias_control()", call. = FALSE)
  if (!identical(center, "bias-corrected") || !identical(method, "bootstrap"))
    stop("density bias.control requires center='bias-corrected' and errors='bootstrap'", call. = FALSE)
  if (identical(control[["construction", exact = TRUE]], "legacy")) return(NULL)
  if (!identical(bootstrap, "inid"))
    stop("positive-covariance density bounds require bootstrap='inid'; other laws remain unqualified", call. = FALSE)
  if (!band %in% c("pointwise", "simultaneous", "bonferroni", "all"))
    stop("positive-covariance density bounds require an explicit quantile band", call. = FALSE)
  .np_udens_bc_prepare(bws, x)
}

.np_udens_plot_errors <- function(state, grid, alpha, B, band) {
  names(grid) <- names(state$x)
  inference <- .np_plot_activity_run("Constructing bias-corrected density intervals",
    .np_udens_b4_inference(state, grid, law = "inid", alpha = alpha, B = B))
  inference$settings$primary <- "absolute"
  inference$bands <- inference$absolute.intervals
  all <- lapply(inference$bands, function(z)
    cbind(inference$center - z[, 1L, 1L], z[, 2L, 1L] - inference$center))
  width <- all[[if (band == "all") "pointwise" else band]]
  list(boot.err = cbind(width, inference$center),
       boot.all.err = if (band == "all") all else NULL,
       inference = inference, bxp = list())
}
