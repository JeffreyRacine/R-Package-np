# Plot-only normal-quantile reduction of a family-owned corrected scale.
# Keep the provider's calibrated bands and provenance; pmzsd is a distinct rule.
.np_plot_bias_normal <- function(inference, alpha) {
  center <- inference[["center", exact = TRUE]]
  scale <- inference[["se", exact = TRUE]]
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= .5)
    stop("normal bias intervals require alpha in (0, 0.5)", call. = FALSE)
  if (!is.numeric(center) || !length(center) || any(!is.finite(center)) ||
      !is.numeric(scale) || length(scale) != length(center) ||
      any(!is.finite(scale)) || any(scale < 0))
    stop("normal bias intervals require finite centers and nonnegative corrected scales", call. = FALSE)
  critical <- qnorm(alpha/2, lower.tail = FALSE)
  width <- critical * scale
  endpoints <- array(c(center - width, center + width),
    c(length(center), 2L, 1L),
    dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
  if (any(!is.finite(endpoints)))
    stop("nonfinite normal bias interval endpoints", call. = FALSE)
  inference$bands$pmzsd <- endpoints
  inference$normal.band <- list(construction = "normal-quantile corrected-scale",
    alpha = alpha, critical = critical, scale = "inference$se",
    calibrated.envelope = FALSE)
  inference$settings$display.band <- "pmzsd"
  inference
}
