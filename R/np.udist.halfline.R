# A fixed Gaussian CDF with one known finite support endpoint. The three
# columns of Q span cubic polynomials that vanish at the transformed endpoint.
.np_udist_halfline_rows <- function(x, g, h, support) {
  left <- is.finite(support[1L])
  endpoint <- if (left) support[1L] else support[2L]
  direction <- if (left) 1 else -1
  distance <- direction * (x - endpoint) / h
  n <- length(x)
  J <- P <- matrix(0, length(g), n)
  index <- 0:3
  for (i in seq_along(g)) {
    if (g[i] == endpoint) {
      if (!left) J[i,] <- P[i,] <- 1
      next
    }
    z <- direction * (g[i] - endpoint) / h
    u <- direction * (g[i] - x) / h
    if (left)
      P[i,] <- (pnorm(u) - pnorm(-distance)) / pnorm(distance)
    else
      P[i,] <- pnorm(u, lower.tail = FALSE) / pnorm(distance)
    density <- dnorm(z)
    M <- c(pnorm(z), -density, numeric(5L))
    for (k in 2:6) {
      edge <- if (density == 0) 0 else z^(k-1L)*density
      M[k+1L] <- (k-1L)*M[k-1L]-edge
    }
    Q <- rbind(-min(1,z)*diag(3L), 0) + rbind(0, diag(3L)/max(1,z))
    H <- matrix(M[1L + outer(index, index, "+")], 4L)
    b <- drop(Q %*% solve(crossprod(Q, H %*% Q), Q[1L,]))
    phi <- dnorm(u)
    uf <- u2f <- numeric(n)
    ok <- phi > 0
    uf[ok] <- u[ok]*phi[ok]
    u2f[ok] <- u[ok]*uf[ok]
    F <- pnorm(u)
    value <- (b[1L]+b[3L])*F - (b[2L]+2*b[4L])*phi - b[3L]*uf - b[4L]*u2f
    J[i,] <- if (left) value else 1-value
  }
  list(original = P, corrected = J)
}

# The minimum and maximum have exact marginal probability laws for a
# continuous iid sample. DKW supplies simultaneous empirical CDF bounds at
# any coordinate at which the signed normal interval excludes all probabilities.
.np_udist_halfline_envelope <- function(answer, widths, x, grid, alpha, structural) {
  n <- length(x)
  sample.range <- range(x)
  at.min <- grid == sample.range[1L]
  at.max <- grid == sample.range[2L]
  below <- grid < sample.range[1L]
  above <- grid > sample.range[2L]
  binding <- vector("list", length(alpha))
  center <- answer$center
  answer$reference.critical <- answer$critical
  for (j in seq_along(alpha)) {
    low <- center - widths$pointwise[, j]
    high <- center + widths$pointwise[, j]
    a <- alpha[j]/2
    normal.low <- low
    normal.high <- high
    # Monotonicity transfers the minimum and maximum rank events to all
    # queries strictly beyond them, including data-selected queries.
    if (any(below)) {
      low[below] <- pmin(low[below], 0)
      high[below] <- pmax(high[below], stats::qbeta(a, 1, n, lower.tail=FALSE))
    }
    if (any(above)) {
      low[above] <- pmin(low[above], stats::qbeta(a, n, 1))
      high[above] <- pmax(high[above], 1)
    }
    for (extreme in c(1L,n)) {
      i <- if (extreme == 1L) at.min else at.max
      if (any(i)) {
        low[i] <- pmin(low[i], stats::qbeta(a, extreme, n+1L-extreme))
        high[i] <- pmax(high[i], stats::qbeta(a, extreme, n+1L-extreme, lower.tail=FALSE))
      }
    }
    rank.low <- low < normal.low
    rank.high <- high > normal.high
    if (length(structural)) rank.low[structural] <- rank.high[structural] <- FALSE
    invalid <- (high < 0 | low > 1) & !below & !above
    binding[[j]] <- list(alpha=alpha[j], rank.lower=which(rank.low),
                         rank.upper=which(rank.high), feasibility=which(invalid),
                         structural=as.integer(structural))
    if (any(invalid)) {
      F <- findInterval(grid[invalid], sort(x))/n
      radius <- sqrt(-log(a)/(2*n))
      low[invalid] <- pmin(low[invalid], pmax(0,F-radius))
      high[invalid] <- pmax(high[invalid], pmin(1,F+radius))
    }
    for (scope in names(answer$intervals)) {
      v <- answer$intervals[[scope]]
      v[, 1L, j] <- pmin(v[, 1L, j], low)
      v[, 2L, j] <- pmax(v[, 2L, j], high)
      if (length(structural)) {
        v[structural, 1L, j] <- center[structural]
        v[structural, 2L, j] <- center[structural]
      }
      answer$intervals[[scope]] <- v
      answer$critical[[scope]][, j] <- pmax(center - v[, 1L, j], v[, 2L, j] - center)
    }
  }
  answer$settings$correction <- "anchored one-sided cubic Gaussian CDF"
  answer$settings$calibration.method <- "boundary-rank"
  answer$settings$coverage.status <- "experimental finite-grid continuous iid CDF"
  answer$settings$endpoint.construction <- "Gaussian process with order-statistic and interior feasibility envelopes"
  answer$settings$critical.components <- binding
  answer$settings$structural.endpoints <- as.integer(structural)
  answer$settings$stochastic.family.size <- length(grid)-length(structural)
  answer$settings$draws <- if (length(answer$process.maximum)) answer$settings$B else 0L
  answer$calibration <- list(method = "sample-extrema rank and feasibility envelope")
  answer
}
