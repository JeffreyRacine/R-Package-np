# Experimental CDS1 univariate CDF percentile owner. Existing defaults do not enter.
# The caller owns fresh outer-sample CV; original fitted geometry is immutable.
.np_udist_smoothed_support <- function() {
  a <- stats::pnorm(-sqrt(3)) - sqrt(3) * stats::dnorm(sqrt(3)) / 2
  c(lower = a, upper = 1 - a)
}

.np_udist_smoothed_prepare <- function(fit, x, max.bytes) {
  if (!is.list(fit) || is.null(fit[["bws", exact = TRUE]]) ||
      !is.data.frame(x) || ncol(x) != 1L || !is.numeric(x[[1L]]) ||
      nrow(x) < 2L || any(!is.finite(x[[1L]])) || diff(range(x[[1L]])) <= 0 ||
      !identical(fit[["bws", exact = TRUE]][["type", exact = TRUE]], "fixed"))
    stop("CDF smoothed-percentile requires a fixed continuous univariate fit", call. = FALSE)
  .np_udist_bc_prepare(fit[["bws", exact = TRUE]], x, max.bytes)
}

.np_udist_smoothed_chunk <- function(task, x, h, grid, B, law, blocklen) {
  draw <- .np_udens_smoothed_samples(x, h, task$start, task$bsz, B,
                                     law, blocklen, task$streams)
  if (!identical(draw$end, task$end))
    stop("smoothed CDF replay consumed inconsistent streams", call. = FALSE)
  y <- data.frame(as.vector(draw$y)); names(y) <- names(grid)
  rows <- .np_udist_bc_continuous(y, grid, h, corrected = TRUE)$corrected
  out <- matrix(0, nrow(grid), task$bsz)
  for (j in seq_len(nrow(grid)))
    out[j, ] <- colMeans(matrix(rows[j, ], nrow(x), task$bsz))
  if (any(!is.finite(out))) stop("nonfinite smoothed CDF curves", call. = FALSE)
  out
}

.np_udist_smoothed_draws <- function(tasks, x, h, grid, B, law, blocklen) {
  out <- matrix(0, nrow(grid), B)
  for (task in tasks) {
    jj <- seq.int(task$start, task$start + task$bsz - 1L)
    out[, jj] <- .np_udist_smoothed_chunk(task, x, h, grid, B, law, blocklen)
  }
  out
}

.np_udist_smoothed_stream_check <- function(streams) {
  if (!is.list(streams) || !length(streams) ||
      !all(vapply(streams, function(s) is.integer(s) && length(s) == 7L &&
        !anyNA(s) && s[1L] %% 100L == 7L, logical(1L))))
    stop("CDF smoothed-percentile requires explicit CMRG streams", call. = FALSE)
  invisible(TRUE)
}

.np_udist_smoothed_provider <- function(state, grid, B, streams, law = "inid",
                                       blocklen = NULL, tile = 25L, chunk = 32L) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (existed) base::get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", previous, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  .np_udist_bc_grid(state, grid)
  B <- .np_plot_gaussian_integer(B, "smoothed CDF draws", 9999L)
  tile <- .np_plot_gaussian_integer(tile, "smoothed CDF tile", 25L)
  chunk <- .np_plot_gaussian_integer(chunk, "smoothed CDF chunk", 32L)
  law <- match.arg(law, c("inid", "fixed", "geom"))
  if ((law == "inid" && !is.null(blocklen)) ||
      (law != "inid" && !identical(as.numeric(blocklen), 1)))
    stop("CDF smoothed-percentile admits inid or fixed/geom block length one", call. = FALSE)
  if (!identical(names(streams), c("count", "jitter")))
    stop("invalid smoothed CDF stream names", call. = FALSE)
  .np_udist_smoothed_stream_check(streams)
  n <- state$n; m <- nrow(grid)
  budget <- .np_udist_bc_budget(state, m, B, tile)
  extra <- 8 * (8 * min(m, tile) * n * chunk + 6 * n * chunk + 12 * B + 12 * m) +
    32 * ceiling(B / chunk) * (2 * length(streams$count) + 32)
  budget$per.process <- budget$per.process + extra
  budget$retained.outer.bytes <- 8 * 39 * (16 * m + 600)
  budget$bytes <- budget$per.process * budget$processes + budget$retained.outer.bytes
  if (!is.finite(budget$bytes) || budget$bytes > state$max.bytes)
    stop("smoothed CDF resident memory forecast exceeds budget", call. = FALSE)
  if (!existed) assign(".Random.seed", streams$count, .GlobalEnv)
  x <- as.matrix(state$x)
  tasks <- .np_udens_smoothed_plan(x, state$h, B, law, blocklen, streams, chunk, FALSE)
  get <- function(ii) {
    if (!is.numeric(ii) || !length(ii) || length(ii) > tile || anyNA(ii) ||
        any(ii != as.integer(ii)) || any(ii < 1L | ii > m))
      stop("invalid smoothed CDF target tile", call. = FALSE)
    # Fanout bookkeeping must not advance the caller's random stream.
    existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    previous <- if (existed) base::get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
    kind <- RNGkind()
    on.exit({
      do.call(RNGkind, as.list(kind))
      if (existed) assign(".Random.seed", previous, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    if (!existed) assign(".Random.seed", streams$count, .GlobalEnv)
    .np_udist_smoothed_draws(tasks, x, state$h, grid[ii, , drop = FALSE], B, law, blocklen)
  }
  list(get = get, m = m, B = B, tile = tile, bound = .np_udist_smoothed_support(),
       start = streams, end = attr(tasks, "streams_end", exact = TRUE), budget = budget)
}

.np_udist_smoothed_percentile <- function(fit, x, grid, streams, outer, outer.fit,
                                         B = 9999L, alpha = .05, law = "inid",
                                         blocklen = NULL, max.bytes = 256 * 1024^2,
                                         tile = 25L, chunk = 32L, warn = TRUE,
                                         progress = NULL, prepared = NULL) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (existed) base::get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", previous, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= .5 ||
      !is.logical(warn) || length(warn) != 1L || is.na(warn) ||
      !is.function(outer.fit) || (!is.null(progress) && !is.function(progress)))
    stop("invalid smoothed CDF calibration controls", call. = FALSE)
  if (!is.list(outer) || length(outer) != 39L || !all(vapply(outer, function(z)
      identical(names(z), c("data", "cv", "count", "jitter")), logical(1L))))
    stop("smoothed CDF calibration requires 39 explicit outer stream records", call. = FALSE)
  all.streams <- c(streams, unlist(outer, recursive = FALSE))
  .np_udist_smoothed_stream_check(all.streams)
  if (anyDuplicated(vapply(all.streams, paste, "", collapse = ":")))
    stop("smoothed CDF stream starts must be distinct", call. = FALSE)
  if (!existed) assign(".Random.seed", streams$count, .GlobalEnv)
  state <- if (is.null(prepared)) .np_udist_smoothed_prepare(fit, x, max.bytes) else {
    if (!is.list(prepared) || !identical(prepared$x, x) ||
        !identical(prepared$bws, fit[["bws", exact=TRUE]]) ||
        !identical(prepared$max.bytes, max.bytes))
      stop("prepared CDF percentile state does not match its fitted inputs", call.=FALSE)
    prepared
  }
  .np_udist_bc_grid(state, grid)
  m <- nrow(grid)
  family <- function(y) {
    z <- data.frame(seq(min(y[[1L]]), max(y[[1L]]), length.out = m))
    names(z) <- names(y); z
  }
  if (!identical(sort(grid[[1L]]), family(x)[[1L]]))
    stop("CDS1 requires the complete untrimmed sample-range family", call. = FALSE)
  original <- .np_udist_smoothed_provider(state, grid, B, streams, law, blocklen, tile, chunk)
  if (!requireNamespace("ks", quietly = TRUE))
    stop("smoothed CDF calibration requires ks", call. = FALSE)
  capture <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(force(expr), warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    list(value = value, warnings = messages)
  }
  pilot <- capture(ks::hpi(x[[1L]], nstage = 2, binned = TRUE, deriv.order = 3))
  g <- pilot$value
  if (length(g) != 1L || !is.finite(g) || g <= 0)
    stop("smoothed CDF derivative-three pilot is undefined", call. = FALSE)
  records <- vector("list", 39L)
  .np_plot_progress_outer_run({
  for (j in seq_len(39L)) {
    os <- outer[[j]]; context <- list(index = j, stage = "outer sample", streams = os, g = g)
    records[[j]] <- tryCatch({
      sample <- .np_bias_stream(os$data, {
        donor <- sort(sample.int(state$n, state$n, replace = TRUE))
        y <- data.frame(x[[1L]][donor] + g * stats::rnorm(state$n)); names(y) <- names(x); y
      })
      y <- sample$value; context$data <- y; context$stage <- "outer CDF CV"
      cv <- capture(.np_bias_stream(os$cv, outer.fit(y)))
      fy <- cv$value$value; context$bandwidth <- fy[["bws", exact = TRUE]]
      context$stage <- "outer CDF preparation"
      sy <- .np_udist_smoothed_prepare(fy, y, max.bytes)
      gy <- family(y)
      truth <- vapply(gy[[1L]], function(t) mean(stats::pnorm((t - x[[1L]]) / g)), 0)
      context$grid <- gy; context$stage <- "inner corrected CDF curves and ranks"
      py <- .np_udist_smoothed_provider(sy, gy, B, list(count = os$count, jitter = os$jitter),
                                        law, blocklen, tile, chunk)
      z <- .np_udens_mx1$scan(py, truth, nominal = 1 - alpha, nominal.alpha = alpha)
      list(row = list(index = j, h = sy$h, q.required = z$required), truth = truth,
           support = py$bound, depth = z$target.depth, kS = z$k$simultaneous,
           grid = gy, bandwidth = sy$h, data.end = sample$end, cv.end = cv$value$end,
           draw.end = py$end, warnings = cv$warnings)
    }, error = function(e) {
      attr(e, "cds1.context") <- context
      attr(e, "cds1.completed.outer") <- records[seq_len(j - 1L)]
      stop(e)
    })
    if (!is.null(progress)) progress(j, 39L)
  }
  })
  reduced <- .np_udens_mx1$reduce(original, records, warn = FALSE, nominal = 1 - alpha, nominal.alpha = alpha)
  reduced$calibration <- "CDF derivative-three pilot with CDF support completion"
  messages <- unique(c(reduced$messages, pilot$warnings,
                       unlist(lapply(records, function(z) z$warnings), use.names = FALSE)))
  resolution <- list(B = B, nominal = 1 - alpha, family.size = m,
    order.index = reduced$order.index, nominal.tail.draws = c(pointwise = B * alpha / 2,
    bonferroni = B * alpha / (2 * m)), limited = any(reduced$order.index == 1))
  if (warn && length(messages)) warning(structure(list(message = paste(messages, collapse = "; "),
    call = NULL, diagnostics = resolution), class = c("np_calibrated_percentile_warning", "warning", "condition")))
  list(inference = reduced, pilot = list(g = g, derivative = 3L), outer = records,
       rng = list(start = original$start, end = original$end), grid = grid,
       budget = original$budget, warnings = messages, resolution = resolution,
       settings = list(construction = "cdf-calibrated-unshifted-smoothed-percentile-v1",
         B = B, alpha = alpha, outer = 39L, rank = as.integer(reduced$selection$k),
         law = law, blocklen = blocklen, h = state$h,
         scope = "private q1 fixed Gaussian2 engineering prototype; coverage unassessed"))
}
