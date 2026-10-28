# Private distribution-specific inference controls; no regression pilot state.
np_udist_bias_control <- function(construction = c("joint-process", "cdf", "legacy")) {
  structure(list(construction = match.arg(construction), version = 1L),
            class = "np_udist_bias_control")
}

.np_udist_plot_prepare <- function(control, bws, x, center, method, bootstrap, band) {
  if (is.null(control)) {
    if (!identical(center, "bias-corrected")) return(NULL)
    control <- np_udist_bias_control()
  }
  construction <- if (is.list(control)) control[["construction", exact = TRUE]] else NULL
  if (!inherits(control, "np_udist_bias_control") || length(construction) != 1L ||
      !is.character(construction) || is.na(construction) ||
      !construction %in% c("joint-process", "cdf", "legacy") ||
      !identical(control, np_udist_bias_control(construction)))
    stop("bias.control must be created by np_udist_bias_control()", call. = FALSE)
  if (!identical(center, "bias-corrected") || !method %in% c("asymptotic", "bootstrap"))
    stop("distribution bias.control requires center='bias-corrected' and asymptotic or bootstrap errors", call. = FALSE)
  if (identical(construction, "legacy")) {
    if (!identical(method, "bootstrap"))
      stop("legacy distribution bias correction requires errors='bootstrap'", call. = FALSE)
    return(NULL)
  }
  if (!identical(bws[["type", exact = TRUE]], "fixed"))
    stop("new distribution bounds for nearest-neighbor bandwidths remain unqualified", call. = FALSE)
  if (identical(method, "bootstrap") && !identical(bootstrap, "inid"))
    stop("new distribution bounds require bootstrap='inid'; dependent laws remain unqualified", call. = FALSE)
  if (!band %in% c("pointwise", "simultaneous", "bonferroni", "all"))
    stop("new distribution bounds require an explicit quantile band", call. = FALSE)
  state <- .np_udist_bc_prepare(bws, x)
  if (identical(construction, "joint-process") &&
      (!state$q %in% 0:2 || !state$s %in% 0:2))
    stop("joint CDF process currently admits at most two numeric and two ordered variables; use an explicit qualified construction", call. = FALSE)
  state$plot.construction <- construction
  state$plot.law <- if (identical(method, "asymptotic")) "gaussian" else "inid"
  state
}

.np_udist_plot_errors <- function(state, grid, alpha, B, band) {
  names(grid) <- names(state$x)
  infer <- if (identical(state$plot.construction, "joint-process"))
    .np_udist_joint_inference else .np_udist_bc_inference
  inference <- .np_plot_activity_run("Constructing bias-corrected distribution intervals",
    infer(state, grid, law = state$plot.law, alpha = alpha, B = B,
      bands = if (band == "all") c("pointwise", "simultaneous", "bonferroni") else band))
  inference$bands <- inference$intervals
  all <- lapply(inference$bands, function(z)
    cbind(inference$center - z[, 1L, 1L], z[, 2L, 1L] - inference$center))
  width <- all[[if (band == "all") "pointwise" else band]]
  list(err = cbind(width, inference$center),
       all.err = if (band == "all") all else NULL,
       center = inference$center, inference = inference, bxp = list())
}
