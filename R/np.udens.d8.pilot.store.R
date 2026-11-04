# Invocation-owned, column-major factor workspace. No persistent numerical cache.
.np_udens_d8_store <- function(m) {
  path <- tempfile("np-d8-factor-")
  connection <- file(path, open = "w+b")
  store <- new.env(parent = emptyenv())
  store$path <- path; store$connection <- connection
  store$row.path <- NULL; store$row.connection <- NULL; store$rows.ready <- FALSE
  store$m <- m; store$rank <- 0L; store$closed <- FALSE
  store
}

.np_udens_d8_store_close <- function(store) {
  if (!store$closed) {
    store$closed <- TRUE
    tryCatch(close(store$connection), finally = {
      unlink(store$path)
      if (!is.null(store$row.connection))
        tryCatch(close(store$row.connection), finally = unlink(store$row.path))
    })
  }
  invisible(NULL)
}

.np_udens_d8_store_append <- function(store, column) {
  if (length(column) != store$m || any(!is.finite(column)))
    stop("invalid D8 factor column", call. = FALSE)
  offset <- 8 * as.double(store$m) * store$rank
  if (!is.finite(offset) || offset + 8 * store$m > 2^53)
    stop("D8 factor workspace exceeds exact file offsets", call. = FALSE)
  seek(store$connection, where = offset, origin = "start", rw = "write")
  writeBin(as.double(column), store$connection, size = 8L)
  store$rank <- store$rank + 1L
  invisible(NULL)
}

.np_udens_d8_store_panel <- function(store, first, last) {
  seek(store$connection, where = 8 * as.double(store$m) * (first - 1),
       origin = "start", rw = "read")
  count <- store$m * (last - first + 1L)
  values <- readBin(store$connection, "double", n = count, size = 8L)
  if (length(values) != count) stop("incomplete D8 factor workspace", call. = FALSE)
  matrix(values, store$m)
}

.np_udens_d8_store_rows <- function(store, first, last) {
  width <- last - first + 1L
  if (isTRUE(store$rows.ready)) {
    seek(store$row.connection, where = 8 * as.double(first-1L) * store$rank,
      origin = "start", rw = "read")
    count <- as.double(width) * store$rank
    values <- readBin(store$row.connection, "double", n = count, size = 8L)
    if (length(values) != count) stop("incomplete D8 row workspace", call. = FALSE)
    return(matrix(values, width, store$rank, byrow = TRUE))
  }
  result <- matrix(0, width, store$rank)
  for (j in seq_len(store$rank)) {
    seek(store$connection, where = 8 * (as.double(store$m) * (j - 1) + first - 1),
         origin = "start", rw = "read")
    values <- readBin(store$connection, "double", n = width, size = 8L)
    if (length(values) != width) stop("incomplete D8 factor workspace", call. = FALSE)
    result[, j] <- values
  }
  result
}

# Transpose only once, in bounded tiles. The factor remains column-major while
# it is being built; all subsequent draw batches use contiguous row reads.
.np_udens_d8_store_for_draws <- function(store, tile) {
  if (store$closed || store$rows.ready || !is.null(store$row.connection))
    stop("invalid D8 workspace transition", call. = FALSE)
  store$row.path <- tempfile("np-d8-rows-")
  store$row.connection <- file(store$row.path, open = "w+b")
  for (first in seq.int(1L, store$m, by = tile)) {
    block <- .np_udens_d8_store_rows(store, first, min(store$m, first+tile-1L))
    writeBin(as.double(t(block)), store$row.connection, size = 8L)
  }
  store$rows.ready <- TRUE
  invisible(NULL)
}

.np_udens_d8_pilot_budget <- function(state, m, B, levels, tile, panel, batch,
                                       processes = 1L) {
  automatic <- is.null(tile)
  tile <- if (automatic) min(m, 32L) else
    min(m, .np_plot_gaussian_integer(tile, "D8 target tile", 256L))
  panel <- min(m, .np_plot_gaussian_integer(panel, "D8 factor panel", 64L))
  batch <- min(B, .np_plot_gaussian_integer(batch, "D8 draw batch", 9999L))
  # Worst rank m. Charges include live inputs, output payload, kernel-expression
  # temporaries, duplicated/serialized factors and actual innovation tasks.
  fixed <- 8 * (6 * state$n * state$q + 40 * m + 24 * m * levels +
                64 * 32 * 256 + 8 * B)
  factor <- 8 * (4 * m * panel + 16 * m)
  draws <- 8 * (8 * m * batch + 6 * tile * m + 6 * tile * batch)
  bytes <- (fixed + max(factor, draws)) * processes
  if (!is.finite(bytes) || bytes > state$max.bytes)
    stop("D8 pilot work exceeds the requested byte budget", call. = FALSE)
  list(bytes = bytes, per.process = bytes/processes, processes = processes,
       tile = tile, panel = panel, batch = batch, automatic.tile = automatic,
       fixed.bytes = fixed, factor.bytes = factor,
       maximum.factor.disk.bytes = 8 * as.double(m)^2,
       maximum.workspace.disk.bytes = 16 * as.double(m)^2,
       representation = "bounded resident panels; invocation-owned factor file")
}

.np_udens_d8_pilot_draw_budget <- function(state, m, rank, budget) {
  if (length(rank) != 1L || !is.numeric(rank) || !is.finite(rank) ||
      rank < 1L || rank > m || rank != floor(rank))
    stop("invalid completed D8 factor rank", call. = FALSE)
  fixed <- budget$fixed.bytes
  base <- fixed + 8 * 8 * as.double(rank) * budget$batch
  per.row <- 8 * (6 * as.double(rank) + 6 * budget$batch)
  room <- floor((state$max.bytes/budget$processes - base)/per.row)
  if (!is.finite(room) || room < 1L)
    stop("D8 completed-factor draws exceed the requested byte budget", call. = FALSE)
  tile <- if (budget$automatic.tile) min(m, room) else budget$tile
  bytes <- max(fixed + budget$factor.bytes, base + tile * per.row)
  if (!is.finite(bytes) || bytes * budget$processes > state$max.bytes)
    stop("D8 completed-factor draws exceed the requested byte budget", call. = FALSE)
  budget$pre.factor.bytes <- budget$bytes
  budget$tile <- as.integer(tile)
  budget$rank <- as.integer(rank)
  budget$per.process <- bytes
  budget$bytes <- bytes * budget$processes
  budget
}

.np_udens_d8_pilot_factor <- function(state, grid, moments, store, panel,
                                      tolerance = 1e-6) {
  m <- nrow(grid); residual <- rep(1, m); pivots <- integer(m)
  while (max(residual) > tolerance && store$rank < m) {
    pivot <- which.max(residual)
    column <- .np_udens_d8_pilot_column(state, grid, moments$mean, pivot) /
      (moments$se * moments$se[pivot])
    if (store$rank) for (first in seq.int(1L, store$rank, by = panel)) {
      block <- .np_udens_d8_store_panel(store, first, min(store$rank, first+panel-1L))
      column <- column - as.vector(block %*% block[pivot, ])
    }
    column <- column/sqrt(residual[pivot])
    updated <- residual - column^2
    if (any(!is.finite(updated)) || min(updated) < -1e-8)
      stop("D8 pilot factor lost positive semidefiniteness", call. = FALSE)
    .np_udens_d8_store_append(store, column)
    pivots[store$rank] <- pivot
    residual <- pmax(updated, 0)
    .np_plot_progress_notify(stage = "Factoring complete D8 pilot covariance",
      done = store$rank, total = m, unit = "columns")
  }
  list(rank = store$rank, pivots = pivots[seq_len(store$rank)], residual = residual,
       tolerance = tolerance, disk.bytes = 8 * as.double(m) * store$rank)
}
