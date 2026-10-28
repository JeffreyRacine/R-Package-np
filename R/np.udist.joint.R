# CDF absolute-process calibration. Critical values are in CDF units.
# The complete centered observation factor is shared across every target tile.
.np_udist_joint_inference <- function(state, grid, law = "gaussian", alpha = .05,
                                     B = 9999L, tile = 25L,
                                     bands = c("pointwise", "simultaneous", "bonferroni"),
                                     corrected = TRUE) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  entry <- if (existed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); success <- FALSE
  on.exit(if (!success) {
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", entry, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  law <- match.arg(law, c("gaussian", "inid"))
  if (!identical(state$type, "fixed"))
    stop("absolute CDF process requires fixed bandwidths", call. = FALSE)
  bands <- match.arg(bands, c("pointwise", "simultaneous", "bonferroni"), several.ok = TRUE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= .5) || anyDuplicated(alpha))
    stop("invalid CDF alpha", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "CDF draws", .Machine$integer.max)
  if (B < 2L) stop("CDF draws must be at least two", call. = FALSE)
  .np_udist_bc_grid(state, grid)
  m <- nrow(grid); budget <- .np_udist_bc_budget(state, m, B, tile)
  stochastic <- law != "gaussian" || any(bands != "pointwise")
  if (stochastic && !existed) runif(1L)
  start <- if (stochastic) get(".Random.seed", .GlobalEnv, inherits = FALSE) else entry
  end <- NULL
  center <- original <- se <- numeric(m)
  limits <- .np_bias_limits_prepare(rep("cdf", m), alpha, law == "gaussian", B)
  for (block in .np_bias_target_blocks(m, 1L, anchors = tile)) {
    ii <- block$grid
    rows <- .np_udist_bc_rows(state, grid[ii, , drop = FALSE], corrected)
    center[ii] <- rows$center; original[ii] <- rows$original.center; se[ii] <- rows$se
    if (stochastic) {
      stream <- .np_bias_stream(start,
        .np_bias_process_draws(rows$factor, B, law, NULL))
      if (!is.null(end) && !identical(end, stream$end))
        stop("absolute CDF tile RNG mismatch", call. = FALSE)
      end <- stream$end
      limits <- .np_bias_limits_update(limits, abs(stream$value), ii)
    }
  }
  if (any(!is.finite(se)) || !any(se > 0))
    stop("absolute CDF process unavailable: no positive observed variance", call. = FALSE)
  raw <- if (stochastic) unname(quantile(limits$maximum, 1-alpha, type = 1L)) else NULL
  if (law == "gaussian") {
    point <- outer(se, qnorm(alpha/2, lower.tail = FALSE))
    bonf <- outer(se, qnorm(alpha/(2*m), lower.tail = FALSE))
    joint <- if (stochastic) pmin(apply(bonf, 2L, max),
      pmax(apply(point, 2L, max), raw)) else apply(point, 2L, max)
  }
  if (law == "inid") {
    point <- limits$point; bonf <- limits$bonf; joint <- raw
  }
  joint <- matrix(rep(joint, each = m), m, length(alpha))
  raw.bonf <- bonf; bonf <- pmax(bonf, joint)
  widths <- list(pointwise = point, simultaneous = joint, bonferroni = bonf)
  intervals <- lapply(widths[bands], function(w) {
    out <- array(0, c(m, 2L, length(alpha)),
      dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    out[, 1L, ] <- center-w; out[, 2L, ] <- center+w
    out
  })
  if (stochastic) assign(".Random.seed", end, .GlobalEnv) else end <- entry
  success <- TRUE
  list(center = center, original.center = original, se = se,
    intervals = intervals, critical = widths[bands], raw.max = raw,
    raw.bonf = raw.bonf, envelope = bonf > raw.bonf,
    process.maximum = if (stochastic) drop(limits$maximum) else NULL,
    zero.marginals = which(se == 0), degeneracy = integer(m), cancellation = integer(m),
    settings = list(construction = "cdf-absolute-process-v1", critical.units = "CDF",
      bwtype = state$type, h = state$h, law = law, n = state$n, alpha = alpha,
      B = B, corrected = corrected, budget = budget,
      representation = "bounded observation factors", metadata = state$metadata,
      target.ids = rownames(grid), scope = "declared finite family"),
    rng = list(kind = kind, entry = entry, start = start, end = end))
}
