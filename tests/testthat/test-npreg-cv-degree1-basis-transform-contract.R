library(npRmpi)

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

build_glp_terms_test <- function(degree) {
  k <- length(degree)
  if (k == 0L)
    return(matrix(integer(0), nrow = 1L, ncol = 0L))

  degree <- as.integer(degree)
  z <- as.matrix(do.call(expand.grid, lapply(degree, function(d) 0:d)))
  s <- rowSums(z)
  z <- z[(s > 0L) & (s <= max(degree)), , drop = FALSE]

  if (!all(degree == max(degree))) {
    for (j in seq_along(degree)) {
      d <- degree[j]
      if ((d < max(degree)) && (d > 0L)) {
        s <- rowSums(z)
        dropj <- (s > d) &
          (z[, j, drop = FALSE] == matrix(d, nrow(z), 1L, byrow = TRUE))
        z <- z[!dropj, , drop = FALSE]
      }
    }
  }

  rbind(matrix(0L, nrow = 1L, ncol = k), z)
}

eval_monomial_basis_test <- function(X, terms) {
  X <- as.matrix(X)
  out <- matrix(1, nrow = nrow(X), ncol = nrow(terms))
  for (j in seq_len(nrow(terms))) {
    exponents <- terms[j, ]
    nz <- which(exponents > 0L)
    if (length(nz)) {
      out[, j] <- apply(
        X[, nz, drop = FALSE]^matrix(exponents[nz], nrow(X), length(nz), byrow = TRUE),
        1L,
        prod
      )
    }
  }
  out
}

build_shift_raw_from_center_test <- function(terms, xj) {
  keys <- setNames(seq_len(nrow(terms)), apply(terms, 1L, paste, collapse = ","))
  out <- matrix(0, nrow = nrow(terms), ncol = nrow(terms))

  recurse <- function(alpha, pos = 1L, beta = integer(length(alpha)), acc = list()) {
    if (pos > length(alpha)) {
      acc[[length(acc) + 1L]] <- beta
      return(acc)
    }
    for (v in 0:alpha[pos]) {
      beta[pos] <- v
      acc <- recurse(alpha, pos + 1L, beta, acc)
    }
    acc
  }

  for (r in seq_len(nrow(terms))) {
    alpha <- terms[r, ]
    for (beta in recurse(alpha)) {
      cidx <- keys[[paste(beta, collapse = ",")]]
      out[r, cidx] <- out[r, cidx] + prod(choose(alpha, beta) * xj^(alpha - beta))
    }
  }

  out
}

solve_glp_raw_center_pair <- function(X, y, degree, bw, extra = rep(1, nrow(X))) {
  terms <- build_glp_terms_test(degree)
  B <- eval_monomial_basis_test(X, terms)
  n <- nrow(B)
  raw_fit <- numeric(n)
  centered_fit <- numeric(n)
  transformed_fit <- numeric(n)

  for (j in seq_len(n)) {
    idx <- setdiff(seq_len(n), j)
    u <- sweep(X[idx, , drop = FALSE], 2L, X[j, ], "-")
    w <- exp(-0.5 * rowSums((u / matrix(bw, nrow(u), ncol(u), byrow = TRUE))^2)) * extra[idx]
    Bj <- B[idx, , drop = FALSE]
    sys_raw <- list(
      s = crossprod(Bj, Bj * w),
      t = crossprod(Bj, y[idx] * w)
    )
    beta_raw <- solve(sys_raw$s, sys_raw$t)
    raw_fit[j] <- drop(B[j, ] %*% beta_raw)

    Cj <- eval_monomial_basis_test(u, terms)
    sys_center <- list(
      s = crossprod(Cj, Cj * w),
      t = crossprod(Cj, y[idx] * w)
    )
    beta_center <- solve(sys_center$s, sys_center$t)
    centered_fit[j] <- beta_center[1]

    shift <- build_shift_raw_from_center_test(terms, X[j, ])
    inv_shift <- solve(shift)
    s_transform <- inv_shift %*% sys_raw$s %*% t(inv_shift)
    t_transform <- inv_shift %*% sys_raw$t
    beta_transform <- solve(s_transform, t_transform)
    transformed_fit[j] <- beta_transform[1]
  }

  list(
    raw_fit = raw_fit,
    centered_fit = centered_fit,
    transformed_fit = transformed_fit,
    raw_cv = (y - raw_fit)^2,
    centered_cv = (y - centered_fit)^2,
    transformed_cv = (y - transformed_fit)^2
  )
}

press_global_glp_test <- function(X, y, degree) {
  terms <- build_glp_terms_test(degree)
  B <- eval_monomial_basis_test(X, terms)
  XtX_inv <- solve(crossprod(B))
  beta <- XtX_inv %*% crossprod(B, y)
  fit <- drop(B %*% beta)
  hdiag <- rowSums((B %*% XtX_inv) * B)
  sum(((y - fit) / (1 - hdiag))^2)
}

const_weight_loocv_glp_test <- function(X, y, degree) {
  terms <- build_glp_terms_test(degree)
  B <- eval_monomial_basis_test(X, terms)
  n <- nrow(B)
  fit <- numeric(n)
  for (j in seq_len(n)) {
    idx <- setdiff(seq_len(n), j)
    beta <- solve(crossprod(B[idx, , drop = FALSE]), crossprod(B[idx, , drop = FALSE], y[idx]))
    fit[j] <- drop(B[j, ] %*% beta)
  }
  sum((y - fit)^2)
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

test_that("higher-degree GLP raw, centered, and transformed CVLS fits agree", {
  set.seed(20260306)
  X <- cbind(
    x1 = runif(24L, min = -1, max = 1),
    x2 = runif(24L, min = -1, max = 1),
    x3 = runif(24L, min = -1, max = 1)
  )
  y <- 0.3 + X[, 1] - 0.4 * X[, 2] * X[, 3] + 0.25 * X[, 1]^2 + rnorm(24L, sd = 0.02)

  res <- solve_glp_raw_center_pair(X, y, degree = c(2L, 1L, 1L), bw = c(0.75, 0.75, 0.75))

  expect_equal(res$raw_fit, res$centered_fit, tolerance = 1e-10)
  expect_equal(res$raw_fit, res$transformed_fit, tolerance = 1e-10)
  expect_equal(res$raw_cv, res$centered_cv, tolerance = 1e-10)
  expect_equal(res$raw_cv, res$transformed_cv, tolerance = 1e-10)
})

test_that("higher-degree GLP basis-change invariant survives mixed-data style scalar modifiers", {
  set.seed(20260306)
  X <- cbind(
    x1 = runif(22L, min = -1, max = 1),
    x2 = runif(22L, min = -1, max = 1)
  )
  y <- 0.5 + X[, 1] + 0.2 * X[, 2]^2 + rnorm(22L, sd = 0.03)
  u <- factor(sample(c("a", "b"), nrow(X), replace = TRUE))
  o <- ordered(sample(1:3, nrow(X), replace = TRUE))
  extra <- ifelse(u == "a", 1, 0.45) * (0.8 ^ abs(as.integer(o) - 2L))

  res <- solve_glp_raw_center_pair(X, y, degree = c(2L, 1L), bw = c(0.7, 0.7), extra = extra)

  expect_equal(res$raw_fit, res$centered_fit, tolerance = 1e-10)
  expect_equal(res$raw_fit, res$transformed_fit, tolerance = 1e-10)
})

test_that("large-h GLP CVLS collapses to global polynomial PRESS", {
  set.seed(20260306)
  X <- cbind(
    x1 = runif(20L, min = -1, max = 1),
    x2 = runif(20L, min = -1, max = 1)
  )
  y <- 0.4 + 0.8 * X[, 1] - 0.3 * X[, 2] + 0.25 * X[, 1]^2 + rnorm(20L, sd = 0.01)

  expect_equal(
    const_weight_loocv_glp_test(X, y, degree = c(2L, 2L)),
    press_global_glp_test(X, y, degree = c(2L, 2L)),
    tolerance = 1e-10
  )
})
