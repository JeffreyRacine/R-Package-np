# Same-sample random-radius influence for unconditional ANN estimates.
# This owner is entered only when uncertainty was requested. Point estimates,
# bandwidths, kernel evaluation and computational options are never changed.
.np_ann_se_domain <- function(bws, n, density) {
  if (n < 4L)
    return("the sample is too small for a two-sided NN spacing pilot")
  if (!identical(bws[["ckerbound", exact = TRUE]], "none"))
    return("ANN uncertainty with finite kernel bounds is not yet implemented")
  kernel <- bws[["ckertype", exact = TRUE]]
  if (!kernel %in% c("gaussian", "epanechnikov", "uniform"))
    return("the kernel's ANN support-boundary influence is not yet implemented")
  k <- as.double(unlist(bws[["bandwidth", exact = TRUE]], use.names = FALSE))[bws[["icon", exact = TRUE]]]
  if (!length(k) || any(!is.finite(k)) || any(k != floor(k)) ||
      any(k < 2 | k > n - 3L))
    return("the NN rank has no regular two-sided interior spacing pilot")
  NULL
}

.np_ann_se_geometry <- function(bws, data) {
  x <- as.matrix(data[, bws[["icon", exact = TRUE]], drop = FALSE])
  storage.mode(x) <- "double"
  k <- as.double(unlist(bws[["bandwidth", exact = TRUE]], use.names = FALSE))[bws[["icon", exact = TRUE]]]
  orders <- lapply(seq_len(ncol(x)), function(j) order(x[, j]))
  geometry <- lapply(seq_len(ncol(x)), function(j)
    .Call("C_np_ann_geometry", x[orders[[j]], j], as.integer(k[j]),
          PACKAGE = "np"))
  list(x = x, order = orders, geometry = geometry,
       valid = all(vapply(geometry, function(g) all(g[["status"]] == 0L), logical(1L))))
}

.np_ann_unconditional_se_local <- function(bws, data, evaluation, density,
                                            progress.context = NULL) {
  n <- nrow(data)
  m <- nrow(evaluation)
  if (!m) return(numeric(0L))
  # The public owner already omitted rows and retains their external mapping.
  # Nested npksum must not interpret that marker as another omission request.
  attr(data, "na.action") <- NULL
  attr(evaluation, "na.action") <- NULL
  state <- .np_ann_se_geometry(bws, data)
  if (!state$valid) return(rep.int(NA_real_, m))
  p <- ncol(state$x)
  e <- as.matrix(evaluation[, bws[["icon", exact = TRUE]], drop = FALSE])
  storage.mode(e) <- "double"
  # Reuse the density-family adapter: its normalized ordered Li-Racine
  # spelling differs from the generic kernel-sum spelling. Use raw ranks
  # even when the public bandwidth object displays scale factors.
  kernel.spec <- bws
  kernel.spec[["bw"]] <- as.double(unlist(bws[["bandwidth", exact = TRUE]],
                                        use.names = FALSE))
  kbw <- .np_make_kbandwidth_unconditional(kernel.spec, data)
  uniform <- density && identical(bws[["ckertype", exact = TRUE]], "uniform")
  if (uniform) state$uniform <- .np_ann_uniform_prepare(kbw, data)
  # Bound the entire resident normal + permutation tensor, not each separately.
  tile <- max(1L, min(64L, floor(1048576 / (as.double(n) * (1L + p + as.integer(uniform))))))
  out <- numeric(m)
  for (start in seq.int(1L, m, by = tile)) {
    rows <- seq.int(start, min(m, start + tile - 1L))
    weights <- npksum(bws = kbw, txdat = data,
      exdat = evaluation[rows, , drop = FALSE],
      operator = if (density) "normal" else "integral",
      permutation.operator = if (density) "derivative" else "normal",
      bandwidth.divide = TRUE, return.kernel.weights = TRUE,
      return.derivative.kernel.weights = TRUE,
      .np.internal.bandwidth.divide.weights = TRUE)
    if (uniform) {
      faces <- .np_ann_uniform_faces(state, data,
        evaluation[rows, , drop = FALSE], e[rows, , drop = FALSE])
      out[rows] <- .Call("C_np_ann_variance_faces", state$geometry, state$order,
        state$x, e[rows, , drop = FALSE], weights[["kw", exact = TRUE]],
        weights[["p.kw", exact = TRUE]], density, faces, PACKAGE = "np")
    } else {
      out[rows] <- .Call("C_np_ann_variance", state$geometry, state$order,
        state$x, e[rows, , drop = FALSE], weights[["kw", exact = TRUE]],
        weights[["p.kw", exact = TRUE]], density, PACKAGE = "np")
    }
    weights <- NULL
    if (!is.null(progress.context))
      progress.context$state <- .np_progress_step_at(
        progress.context$state, now = .np_progress_now(), done = max(rows),
        detail = sprintf("asymptotic standard errors (%d/%d rows)", max(rows), m))
  }
  out
}

.np_ann_unconditional_se <- function(bws, data, evaluation, density) {
  reason <- .np_ann_se_domain(bws, nrow(data), density)
  if (!is.null(reason)) {
    .np_warning(paste0("Adaptive-NN standard errors are unavailable: ", reason,
                      ". Point estimates are retained; standard errors are NA."),
                call. = FALSE)
    return(rep.int(NA_real_, nrow(evaluation)))
  }
  progress.context <- new.env(parent = emptyenv())
  progress.context$state <- .np_progress_begin(
    "Adaptive-NN standard errors", total = nrow(evaluation), surface = "fit")
  on.exit(.np_progress_end(progress.context$state), add = TRUE)
  out <- .np_ann_unconditional_se_local(bws, data, evaluation, density, progress.context)
  if (anyNA(out))
    .np_warning(paste0("Adaptive-NN standard errors are unavailable for ",
      sum(is.na(out)), " evaluation row(s): a NN or boundary spacing pilot is ",
      "non-interior or degenerate, or the estimated influence is nonfinite. Point estimates are retained; ",
      "these standard errors are NA."), call. = FALSE)
  out
}
