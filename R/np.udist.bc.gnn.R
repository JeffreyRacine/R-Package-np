# Query-radius CDF influence. Counts are immutable; geometry belongs to a sample.
.np_udist_bc_nn_prepare <- function(bws, x, max.bytes) {
  numeric <- vapply(x, is.numeric, logical(1L))
  if (!identical(bws[["type", exact = TRUE]], "generalized_nn") && any(numeric))
    stop("corrected ANN CDF inference is parked at its population-remainder gate", call. = FALSE)
  requested <- .np_physical_bandwidth(bws)
  package <- getNamespaceName(environment(.np_physical_bandwidth))
  rank <- if (any(numeric)) .Call("C_np_plot_ann_rank_decode", as.double(requested[numeric]),
    as.integer(nrow(x)), PACKAGE = package) else NULL
  if (!is.null(rank) && (any(rank$extended) || any(rank$distance.scale != 1)))
    stop("CDF NN inference requires interior unscaled counts", call. = FALSE)
  spec <- bws; spec[["type"]] <- "fixed"; spec[["scaling"]] <- FALSE
  spec[["bw"]] <- requested; spec[["bw"]][numeric] <- 1
  state <- .np_udist_bc_prepare(spec, x, max.bytes)
  state$type <- bws[["type", exact = TRUE]]; state$bws <- bws
  state$requested <- requested[numeric]; state$rank <- rank
  state$metadata$radius <- "sample-specific query kth distance"
  state$metadata$influence <- "joint direct plus all corrected radius terms"
  state
}

.np_udist_bc_gnn_geometry <- function(state, grid) {
  n <- state$n; q <- state$q; m <- nrow(grid)
  radius <- density <- matrix(0, m, q)
  for (j in seq_len(q)) {
    k <- state$rank$lookup[j]
    w <- min(floor(n^(2/3)), floor(min(k, n-k)/2))
    if (w < 1L || k-w < 1L || k+w > n)
      stop("CDF GNN spacing pilot is non-interior", call. = FALSE)
    x <- state$x[state$numeric][[j]]; targets <- grid[state$numeric][[j]]
    distinct <- unique(targets); r <- d <- numeric(length(distinct))
    for (a in seq_along(distinct)) {
      sorted <- sort.int(abs(x-distinct[a]), partial = c(k-w, k, k+w))
      r[a] <- sorted[k]; d[a] <- (2*w/n)/(sorted[k+w]-sorted[k-w])
    }
    index <- match(targets, distinct); radius[, j] <- r[index]; density[, j] <- d[index]
  }
  if (any(!is.finite(radius)) || any(radius <= 0) || any(!is.finite(density)) || any(density <= 0))
    stop("CDF GNN radius or spacing pilot is unavailable", call. = FALSE)
  list(radius = radius, density.sum = density)
}

.np_udist_bc_gnn_rows <- function(state, grid, corrected) {
  .np_udist_bc_grid(state, grid)
  if (nrow(grid) > 25L) stop("CDF NN rows require bounded target tiles", call. = FALSE)
  .np_udist_bc_budget(state, nrow(grid), 1L, nrow(grid))
  geometry <- .np_udist_bc_gnn_geometry(state, grid)
  x <- state$x[state$numeric]; g <- grid[state$numeric]; q <- state$q
  # Recurrence also evaluates q=0 exactly, without an NN geometry requirement.
  continuous <- function(skip = integer()) {
    P <- matrix(1, nrow(grid), state$n); Q <- P*0
    for (j in setdiff(seq_len(q), skip)) {
      u <- outer(g[[j]], x[[j]], "-")/geometry$radius[, j]
      H <- pnorm(u); moment <- u*dnorm(u); moment[is.infinite(u)] <- 0
      Q <- Q*H + P*moment; P <- P*H
    }
    list(P = P, Q = Q)
  }
  con <- continuous(); original <- con$P
  L <- if (corrected) con$P+.5*con$Q else con$P
  categorical <- matrix(1, nrow(grid), state$n)
  for (cat in state$category) {
    distance <- outer(as.double(toMatrix(grid[cat$index])), cat$training, "-")
    original <- original*.np_udist_bc_ordered(distance, cat$lambda, FALSE)
    categorical <- categorical*.np_udist_bc_ordered(distance, cat$lambda, corrected)
  }
  L <- L*categorical; phi <- L
  for (j in seq_len(q)) {
    h <- geometry$radius[, j]; u <- outer(g[[j]], x[[j]], "-")/h
    other <- continuous(j)
    D <- -u/h*dnorm(u)*if (corrected) ((1.5-.5*u^2)*other$P+.5*other$Q) else other$P
    D <- D*categorical
    inside <- abs(outer(g[[j]], x[[j]], "-")) <= h
    phi <- phi-rowMeans(D)/geometry$density.sum[, j]*inside
  }
  shifted <- phi-phi[, 1L]; F <- (shifted-rowMeans(shifted))/state$n
  if (any(!is.finite(L)) || any(!is.finite(F))) stop("nonfinite CDF GNN influence", call. = FALSE)
  list(original = original, corrected = L,
    original.center = original[, 1L]+rowMeans(original-original[, 1L]),
    center = L[, 1L]+rowMeans(L-L[, 1L]), factor = F, se = sqrt(rowSums(F^2)), geometry = geometry)
}
