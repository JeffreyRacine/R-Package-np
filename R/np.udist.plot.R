# Private distribution-specific inference controls; no regression pilot state.
np_udist_bias_control <- function(construction = c("joint-process", "cdf", "legacy")) {
  structure(list(construction = match.arg(construction), version = 1L),
            class = "np_udist_bias_control")
}

.np_udist_plot_prepare <- function(control, bws, x, center, method, bootstrap, band, blocklen = NULL, block.automatic = FALSE) {
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
  nearest <- !identical(bws[["type", exact = TRUE]], "fixed")
  if (nearest && !identical(construction, "joint-process"))
    stop("explicit cdf construction for nearest-neighbor bandwidths remains unqualified", call. = FALSE)
  if (nearest && identical(method, "bootstrap") && !bootstrap %in% c("inid", "fixed", "geom"))
    stop("selected nearest-neighbor CDF bootstrap law remains unqualified", call. = FALSE)
  if (identical(method, "bootstrap") &&
      !bootstrap %in% if (construction == "joint-process") c("inid", "fixed", "geom") else "inid")
    stop("selected bootstrap law is not qualified for this distribution construction", call. = FALSE)
  if (!band %in% c("pmzsd", "pointwise", "simultaneous", "bonferroni", "all"))
    stop("unsupported bias-corrected distribution band", call. = FALSE)
  state <- .np_udist_bc_prepare(bws, x)
  if (identical(construction, "joint-process") &&
      (!state$q %in% 0:2 || !state$s %in% 0:2))
    stop("joint CDF process currently admits at most two numeric and two ordered variables; use an explicit qualified construction", call. = FALSE)
  state$plot.construction <- construction
  state$plot.law <- if (identical(method, "asymptotic")) "gaussian" else bootstrap
  if (state$plot.law %in% c("fixed", "geom")) {
    state$plot.block.selection <- .np_udist_bc_block_request(x, state$plot.law, blocklen)
    if (isTRUE(block.automatic) && !is.null(blocklen)) {
      state$plot.block.selection$automatic <- TRUE
      state$plot.block.selection$source <- "b.star(xdat, round=TRUE)[1,1]"
      state$plot.block.selection$owner <- "shared plot normalization"
    }
  }
  state
}

.np_udist_plot_errors <- function(state, grid, alpha, B, band) {
  names(grid) <- names(state$x)
  provider.band <- if (band == "pmzsd") "pointwise" else band
  infer <- if (identical(state$plot.construction, "joint-process"))
    .np_udist_joint_inference else .np_udist_bc_inference
  inference <- .np_plot_activity_run("Constructing bias-corrected distribution intervals",
    do.call(infer, c(list(state = state, grid = grid, law = state$plot.law,
      alpha = alpha, B = B,
      bands = if (provider.band == "all") c("pointwise", "simultaneous", "bonferroni") else provider.band),
      if (identical(state$plot.construction, "joint-process"))
        list(selection = state$plot.block.selection) else list())))
  inference$bands <- inference$intervals
  if (band == "pmzsd") inference <- .np_plot_bias_normal(inference, alpha)
  all <- lapply(inference$bands, function(z)
    cbind(inference$center - z[, 1L, 1L], z[, 2L, 1L] - inference$center))
  width <- all[[if (band == "all") "pointwise" else band]]
  list(err = cbind(width, inference$center),
       all.err = if (band == "all") all else NULL,
       center = inference$center, inference = inference, bxp = list())
}
