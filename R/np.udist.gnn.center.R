# Exact GNN geometry and center algebra on bounded marginal-coordinate blocks.
.np_udist_gnn_blocks <- function(state, grid, size = 25L) {
  columns <- which(state$numeric)
  categories <- vapply(state$category, `[[`, integer(1), "index")
  groups <- if (length(categories))
    split(seq_len(nrow(grid)), interaction(grid[categories], drop = TRUE)) else
    list(seq_len(nrow(grid)))
  out <- list()
  for (ii in groups) {
    axes <- lapply(columns, function(j) unique(grid[[j]][ii]))
    chunks <- lapply(axes, function(a)
      split(seq_along(a), ceiling(seq_along(a) / size)))
    for (a in chunks[[1L]]) for (b in if (length(axes) == 2L) chunks[[2L]] else list(1L)) {
      coordinates <- list(axes[[1L]][a])
      if (length(axes) == 2L) coordinates[[2L]] <- axes[[2L]][b]
      keep <- grid[[columns[1L]]][ii] %in% coordinates[[1L]]
      if (length(axes) == 2L)
        keep <- keep & grid[[columns[2L]]][ii] %in% coordinates[[2L]]
      jj <- ii[keep]
      if (!length(jj)) next
      # Only requested pairs enter the process, including duplicate targets.
      out[[length(out) + 1L]] <- list(rows = jj, axes = coordinates,
        map = lapply(seq_along(columns), function(k)
          match(grid[[columns[k]]][jj], coordinates[[k]])),
        profile = grid[jj[1L], , drop = FALSE])
    }
  }
  out
}

.np_udist_gnn_block <- function(state, block, corrected = TRUE) {
  columns <- which(state$numeric)
  cache <- lapply(seq_along(columns), function(k) {
    delta <- outer(block$axes[[k]], state$x[[columns[k]]], "-")
    ordering <- t(apply(abs(delta), 1L, order))
    if (nrow(delta) == 1L) ordering <- matrix(ordering, 1L)
    list(delta = delta, order = ordering, k = state$rank$lookup[k])
  })
  categorical <- rep(1, state$n)
  for (category in state$category)
    categorical <- categorical * .np_udist_bc_ordered(
      as.numeric(as.character(block$profile[[category$index]])) - category$training,
      category$lambda, corrected)
  list(cache = cache, category = categorical, map = block$map,
    n = state$n, corrected = corrected)
}

.np_udist_gnn_center <- function(prepared, counts) {
  stopifnot(length(counts) == prepared$n, sum(counts) == prepared$n, all(counts >= 0))
  radii <- list()
  margins <- lapply(prepared$cache, function(cache) {
    radius <- vapply(seq_len(nrow(cache$delta)), function(i) {
      ordering <- cache$order[i, ]
      donor <- ordering[which(cumsum(counts[ordering]) >= cache$k)[1L]]
      abs(cache$delta[i, donor])
    }, numeric(1))
    if (any(!is.finite(radius)) || any(radius <= 0))
      stop("exact GNN center unavailable: nonpositive physical radius", call. = FALSE)
    radii[[length(radii) + 1L]] <<- radius
    z <- cache$delta / radius
    D <- z * dnorm(z); D[is.infinite(z)] <- 0
    list(H = pnorm(z), D = D)
  })
  weights <- counts * prepared$category
  first <- margins[[1L]]
  if (length(margins) == 1L) {
    contribution <- if (prepared$corrected) first$H + .5 * first$D else first$H
    center <- drop(contribution %*% weights) / prepared$n
  } else {
    second <- margins[[2L]]
    weighted.H <- sweep(first$H, 2L, weights, "*")
    weighted.D <- sweep(first$D, 2L, weights, "*")
    center <- (weighted.H %*% t(second$H) + if (prepared$corrected)
      .5 * (weighted.D %*% t(second$H) + weighted.H %*% t(second$D)) else 0) / prepared$n
  }
  list(center = if (length(margins) == 1L) center[prepared$map[[1L]]] else
      center[cbind(prepared$map[[1L]], prepared$map[[2L]])],
    radii = do.call(cbind, lapply(seq_along(radii), function(k)
      radii[[k]][prepared$map[[k]]])))
}

# Preserve the observed center's original recurrence and shifted summation.
# Radius computation here needs no spacing pilot or bootstrap studentizer.
.np_udist_gnn_main <- function(state, grid, corrected = TRUE) {
  n <- state$n; numeric <- which(state$numeric)
  P <- matrix(1, nrow(grid), n); Q <- P * 0
  for (a in seq_along(numeric)) {
    j <- numeric[a]; k <- state$rank$lookup[a]; distinct <- unique(grid[[j]])
    radius <- vapply(distinct, function(t)
      sort.int(abs(state$x[[j]] - t), partial = k)[k], numeric(1))
    h <- radius[match(grid[[j]], distinct)]
    if (any(!is.finite(h)) || any(h <= 0))
      stop("exact GNN center unavailable: nonpositive physical radius", call. = FALSE)
    u <- outer(grid[[j]], state$x[[j]], "-") / h
    H <- pnorm(u); D <- u * dnorm(u); D[is.infinite(u)] <- 0
    Q <- Q * H + P * D; P <- P * H
  }
  original <- P; G <- if (corrected) P + .5 * Q else P
  for (cat in state$category) {
    distance <- outer(as.numeric(as.character(grid[[cat$index]])), cat$training, "-")
    original <- original * .np_udist_bc_ordered(distance, cat$lambda, FALSE)
    G <- G * .np_udist_bc_ordered(distance, cat$lambda, corrected)
  }
  list(center = G[, 1L] + rowMeans(G - G[, 1L]),
    original.center = original[, 1L] + rowMeans(original - original[, 1L]))
}

.np_udist_gnn_budget <- function(state, m, B, processes = 1L) {
  n <- as.double(state$n); p <- ncol(state$x)
  base <- 8 * (16*n*p + 40*n*25 + 8*n*32 + 12*B + 12*m)
  capacity <- floor((state$max.bytes / processes - base) / (6 * 8 * B))
  side <- min(25L, floor(sqrt(max(0, capacity))))
  if (!is.finite(side) || side < 1L)
    stop("exact GNN resident memory forecast exceeds budget", call. = FALSE)
  per.process <- base + 6 * 8 * side^2 * B
  list(bytes = processes * per.process, per.process = per.process,
    processes = processes, marginal.tile = as.integer(side),
    target.tile = as.integer(side^2), draw.chunk = 32L, covariance.cap = 32L)
}
