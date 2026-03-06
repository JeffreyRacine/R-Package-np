library(np)

gaussian_row_weights <- function(x, j, bw, extra = rep(1, length(x))) {
  idx <- setdiff(seq_along(x), j)
  u <- (x[idx] - x[[j]]) / bw
  exp(-0.5 * u * u) * extra[idx]
}

degree1_raw_moments <- function(x, y, j, weights) {
  idx <- setdiff(seq_along(x), j)
  xi <- x[idx]
  yi <- y[idx]
  bi <- cbind(1, xi)
  s <- crossprod(bi * weights, bi)
  t <- crossprod(bi * weights, yi)
  list(s = unclass(s), t = as.numeric(t))
}

degree1_centered_moments <- function(xj, raw) {
  s <- raw$s
  t <- raw$t
  k <- nrow(s)
  out_s <- matrix(0, nrow = k, ncol = k)
  out_t <- numeric(k)

  out_t[1] <- t[1]
  out_s[1, 1] <- s[1, 1]

  for (a in 2:k) {
    out_t[a] <- t[a] - xj[a - 1L] * t[1]
    out_s[1, a] <- s[1, a] - xj[a - 1L] * s[1, 1]
    out_s[a, 1] <- out_s[1, a]
  }

  for (a in 2:k) {
    sa0 <- s[a, 1]
    for (b in 2:k) {
      out_s[a, b] <- s[a, b] -
        xj[a - 1L] * s[1, b] -
        xj[b - 1L] * sa0 +
        xj[a - 1L] * xj[b - 1L] * s[1, 1]
    }
  }

  list(s = out_s, t = out_t)
}

solve_degree1_pair <- function(x, y, bw, extra = rep(1, length(x))) {
  n <- length(x)
  raw_fit <- numeric(n)
  centered_fit <- numeric(n)
  raw_cv <- numeric(n)
  centered_cv <- numeric(n)

  for (j in seq_len(n)) {
    w <- gaussian_row_weights(x, j, bw, extra = extra)
    raw <- degree1_raw_moments(x, y, j, w)
    beta_raw <- solve(raw$s, raw$t)
    centered <- degree1_centered_moments(x[j], raw)
    beta_centered <- solve(centered$s, centered$t)

    raw_fit[j] <- c(1, x[j]) %*% beta_raw
    centered_fit[j] <- beta_centered[1]
    raw_cv[j] <- (y[j] - raw_fit[j])^2
    centered_cv[j] <- (y[j] - centered_fit[j])^2
  }

  list(
    raw_fit = raw_fit,
    centered_fit = centered_fit,
    raw_cv = raw_cv,
    centered_cv = centered_cv
  )
}

solve_degree1_one <- function(x, y, j, bw, extra = rep(1, length(x))) {
  w <- gaussian_row_weights(x, j, bw, extra = extra)
  raw <- degree1_raw_moments(x, y, j, w)
  beta_raw <- solve(raw$s, raw$t)
  centered <- degree1_centered_moments(x[j], raw)
  beta_centered <- solve(centered$s, centered$t)

  list(
    raw_fit = as.numeric(c(1, x[j]) %*% beta_raw),
    centered_fit = as.numeric(beta_centered[1]),
    raw_cv = (y[j] - as.numeric(c(1, x[j]) %*% beta_raw))^2,
    centered_cv = (y[j] - as.numeric(beta_centered[1]))^2
  )
}

test_that("degree-1 raw-basis and centered solves are identical for continuous fixed-bandwidth LOO fits", {
  set.seed(20260306)
  x <- sort(runif(40L))
  y <- sin(2 * pi * x) + rnorm(40L, sd = 0.03)

  res <- solve_degree1_pair(x, y, bw = 0.19)

  expect_equal(res$raw_fit, res$centered_fit, tolerance = 1e-12)
  expect_equal(res$raw_cv, res$centered_cv, tolerance = 1e-12)
})

test_that("degree-1 basis-change invariant holds under mixed-data style scalar kernel modifiers", {
  set.seed(20260306)
  x <- sort(runif(45L))
  y <- cos(2 * pi * x) + rnorm(45L, sd = 0.04)
  u <- factor(sample(c("a", "b", "c"), length(x), replace = TRUE))
  o <- ordered(sample(1:4, length(x), replace = TRUE))

  extra <- numeric(length(x))
  for (j in seq_along(x)) {
    same_u <- ifelse(u == u[j], 1.0, 0.35)
    ord_gap <- abs(as.integer(o) - as.integer(o[j]))
    mixed_extra <- same_u * (0.8 ^ ord_gap)
    res_j <- solve_degree1_one(x, y, j, bw = 0.17, extra = mixed_extra)
    extra[j] <- abs(res_j$raw_fit - res_j$centered_fit)
    expect_equal(res_j$raw_fit, res_j$centered_fit, tolerance = 1e-12)
    expect_equal(res_j$raw_cv, res_j$centered_cv, tolerance = 1e-12)
  }

  expect_lte(max(extra), 1e-12)
})
