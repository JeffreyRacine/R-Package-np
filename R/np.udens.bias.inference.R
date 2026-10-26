# Private entry, intentionally absent from all public call/dispatch paths.
.np_udens_bc_inference <- function(state, grid, law = "gaussian", alpha = c(.05, .1),
                                    B = 9999L, representation = c("covariance", "factor"),
                                    tile = 25L, bands = c("pointwise", "simultaneous", "bonferroni"),
                                    corrected = TRUE) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (existed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); success <- FALSE
  on.exit({
    if (!success) {
      do.call(RNGkind, as.list(kind))
      if (existed) assign(".Random.seed", previous, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  law <- match.arg(law, "gaussian")
  representation <- match.arg(representation)
  bands <- match.arg(bands, c("pointwise", "simultaneous", "bonferroni"), several.ok = TRUE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= .5) || anyDuplicated(alpha))
    stop("density alpha values must be unique and in (0, .5)", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "density draws", .Machine$integer.max)
  if (B < 2L) stop("density calibration requires at least two draws", call. = FALSE)
  .np_udens_bc_grid(state, grid)
  budget <- .np_udens_bc_budget(state, nrow(grid), B, tile, representation == "covariance")
  if (!existed) runif(1L)
  start <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  value <- if (representation == "covariance")
    .np_udens_bc_gaussian(state, grid, alpha, B, tile, corrected, start) else
    .np_udens_bc_factor(state, grid, alpha, B, tile, corrected, start)
  m <- nrow(grid)
  intervals <- lapply(value$critical, function(critical) {
    out <- array(0, c(m, 2L, length(alpha)),
                 dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    out[, 1L, ] <- value$center - value$se * critical
    out[, 2L, ] <- value$center + value$se * critical
    out
  })
  unbounded <- lapply(value$critical, function(z) is.infinite(z))
  result <- list(original.center = value$original.center, center = value$center,
    se = value$se, intervals = intervals[bands], critical = value$critical[bands],
    raw.max = value$raw.max, raw.bonf = value$raw.bonf,
    envelope = value$critical$bonferroni > value$raw.bonf,
    unbounded = unbounded[bands], degeneracy = value$degeneracy,
    cancellation = value$cancellation, covariance = value$covariance,
    map = data.frame(row = seq_len(m), family = rep("density", m)),
    settings = list(family = "npudens", target = "joint density/mass level",
      construction = state$construction, corrected = corrected, law = law,
      h = state$h, auxiliary.h = state$h, ratio = 1, n = state$n,
      q = state$q, r = state$r, s = state$s, alpha = alpha, B = B,
      representation = representation, covariance = "empirical, n-squared",
      metadata = state$metadata, budget = budget),
    rng = list(kind = kind, entry = previous, start = start, end = value$end))
  assign(".Random.seed", value$end, .GlobalEnv)
  success <- TRUE
  result
}
