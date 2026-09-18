r17_ann_literal_hat <- function(x, k, kernel) {
  n <- length(x)
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    donors <- setdiff(seq_len(n), i)
    h <- vapply(donors, function(j)
      sort(abs(x[-c(i, j)] - x[j]))[k], numeric(1))
    u <- (x[i] - x[donors]) / h
    w <- switch(kernel,
      gaussian = dnorm(u),
      epanechnikov = ifelse(abs(u) < sqrt(5), 3/(4*sqrt(5))*(1-u^2/5), 0),
      uniform = ifelse(abs(u) <= 1, .5, 0)) / h
    H[i, donors] <- w/sum(w)
  }
  H
}

test_that("adaptive compact-kernel LOO hats retain safe tree geometry", {
  withr::local_options(np.messages = FALSE, np.tree = TRUE)
  set.seed(19318)
  x <- data.frame(x = runif(24, -1, 1))
  y <- sin(x$x)
  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    b <- npregbw(xdat = x, ydat = y, bws = 12,
      regtype = "lc", bwtype = "adaptive_nn", ckertype = kernel,
      bandwidth.compute = FALSE)
    options(np.tree = TRUE)
    H <- npreghat(b, txdat = x, leave.one.out = TRUE)
    a <- npreghat(b, txdat = x, y = cbind(y, y^2),
                  output = "apply", leave.one.out = TRUE)
    options(np.tree = FALSE)
    dense <- npreghat(b, txdat = x, leave.one.out = TRUE)
    expect_equal(unname(H), unname(dense), tolerance = 1e-10)
    expect_equal(as.numeric(H), as.numeric(r17_ann_literal_hat(x$x, 12, kernel)),
                 tolerance = 1e-10)
    expect_equal(as.numeric(a), as.numeric(H %*% cbind(y, y^2)), tolerance = 1e-10)
    expect_equal(diag(H), rep(0, nrow(x)), tolerance = 1e-10)
  }
})

test_that("single-row tree repair preserves polynomial and bandwidth siblings", {
  withr::local_options(np.messages = FALSE, np.tree = TRUE)
  set.seed(19319)
  x <- data.frame(x = runif(26, -1, 1), z = runif(26, -1, 1))
  y <- sin(x$x) + x$z
  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn"))
    for (degree in 0:2) for (kernel in c("gaussian", "epanechnikov")) {
      b <- npregbw(xdat = x, ydat = y,
        bws = rep(if (bwtype == "fixed") .8 else 16, 2),
        regtype = "lp", degree = rep(degree, 2), bwtype = bwtype,
        ckertype = kernel, bandwidth.compute = FALSE)
      options(np.tree = FALSE)
      a <- npreghat(b, txdat = x, y = y, output = "apply", leave.one.out = TRUE)
      options(np.tree = TRUE)
      z <- npreghat(b, txdat = x, y = y, output = "apply", leave.one.out = TRUE)
      expect_equal(z, a, tolerance = 1e-10)
    }
})
