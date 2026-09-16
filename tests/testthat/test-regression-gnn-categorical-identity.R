# Tiny dense matrices are test oracles only: production retains streaming rows.
test_that("GNN significance response tiles use the same fitted-row categorical target", {
  skip_on_cran()
  spawn_mpi_slaves()
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(8843)
  n <- 42L
  x <- data.frame(x = runif(n, -1, 1),
    u = factor(rep(letters[1:3], length.out = n)),
    o = ordered(rep(letters[1:3], each = n / 3)))
  y <- sin(x$x) + .4 * (x$u == "b") + rnorm(n, sd = .3)
  response <- cbind(y, rev(y), y + rnorm(n, sd = .2))
  tile <- getFromNamespace(".np_npsig_streamed_iid_tile", "npRmpi")
  joint <- getFromNamespace(".np_npsig_streamed_response_statistic", "npRmpi")
  for (type in c("lc", "ll", "lp")) {
    bw <- npregbw(xdat = x, ydat = y, bws = c(28, .3, .4),
      bandwidth.compute = FALSE, bwtype = "generalized_nn", regtype = type,
      degree = if (type == "lp") 2L else NULL)
    direct <- matrix(0, 3L, 2L)
    for (r in 1:3) {
      fit <- npreg(bw, txdat = x, tydat = response[, r],
                    gradients = TRUE, se = TRUE)
      for (j in 2:3) {
        reference <- if (j == 2L) as.integer(x$u) == 1L else rep(FALSE, n)
        ratio <- numeric(n)
        ratio[!reference] <- fit$grad[!reference, j] / fit$gerr[!reference, j]
        direct[r, j - 1L] <- mean(ratio^2)
      }
    }
    for (j in 2:3)
      expect_identical(tile(bw, x, j, response.matrix = response,
        null.mean = y, residual.pool = y, pivotal = TRUE), direct[, j - 1L])
    expect_equal(joint(bw, x, 2:3, response, pivotal = TRUE), rowMeans(direct),
                  tolerance = 2e-12)
    seed <- .Random.seed
    result <- npsigtest(bw, index = 2:3, B = 9L, random.seed = 883)
    expect_identical(.Random.seed, seed)
    expect_equal(unname(result$In), direct[1L, ], tolerance = 2e-12)
    expect_equal(unname(result$P), colMeans(sweep(result$In.bootstrap, 2L,
      result$In, `>=`)), tolerance = 0)
  }
})

c88_exponents <- function(degree, basis) {
  e <- as.matrix(expand.grid(0:degree, 0:degree))
  if (basis == "glp") e <- e[rowSums(e) <= degree, , drop = FALSE]
  if (basis == "additive") e <- e[rowSums(e > 0L) <= 1L, , drop = FALSE]
  e
}

c88_rows <- function(x, z, k, degree, basis, identity) {
  n <- nrow(x)
  e <- c88_exponents(degree, basis)
  a <- matrix(0, nrow(z), n)
  for (i in seq_len(nrow(z))) {
    w <- rep(1, n)
    for (j in 1:2) {
      distance <- abs(x[[j]] - z[[j]][i])
      # The fitted identity is excluded only while selecting the radius.
      # Its response still enters the weights and the full local solve.
      radius <- sort(distance)[k[j] + as.integer(identity)]
      w <- w * dnorm(distance / radius) / radius
    }
    w <- w * ifelse(x$u == z$u[i], .7, .15) *
      .4^abs(as.integer(x$o) - as.integer(z$o[i]))
    dx <- sweep(as.matrix(x[1:2]), 2L, as.numeric(z[i, 1:2]), "-")
    b <- vapply(seq_len(nrow(e)), function(j)
      dx[, 1L]^e[j, 1L] * dx[, 2L]^e[j, 2L], numeric(n))
    gram <- crossprod(b, w * b)
    if (rcond(gram) < 1e-11)
      stop("CF88 oracle fixture is ill conditioned; do not compare ridge policies")
    a[i, ] <- as.vector(b %*% solve(gram, as.numeric(rowSums(e) == 0))) * w
  }
  a
}

c88_endpoints <- function(z, column) {
  lo <- hi <- z
  v <- z[[column]]
  lev <- levels(v)
  if (is.ordered(v)) {
    code <- as.integer(v)
    lo[[column]] <- ordered(lev[pmax(1L, code - 1L)], levels = lev)
    hi[[column]] <- ordered(lev[ifelse(code == 1L, 2L, code)], levels = lev)
  } else lo[[column]] <- factor(rep(lev[1L], nrow(z)), levels = lev)
  list(lo = lo, hi = hi)
}

test_that("GNN categorical influence uses the requested training or external geometry", {
  skip_on_cran()
  spawn_mpi_slaves()
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(8831)
  n <- 72L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("a", "b", "c"), each = n / 3)))
  y <- sin(2*x$x1) + x$x2^2 + .3*(x$u == "b") +
    .2*as.integer(x$o) + rnorm(n, sd = .3)
  k <- c(48L, 53L)
  cases <- rbind(data.frame(type = c("lc", "ll"), degree = 0:1,
                           basis = "glp", bernstein = FALSE),
    expand.grid(type = "lp", degree = 0:3,
                basis = c("glp", "additive", "tensor"),
                bernstein = c(FALSE, TRUE), stringsAsFactors = FALSE))
  for (case in seq_len(nrow(cases))) {
    cfg <- cases[case, ]
    args <- list(xdat = x, ydat = y, bws = c(k, .3, .4),
      bandwidth.compute = FALSE, bwtype = "generalized_nn",
      regtype = cfg$type, ukertype = "aitchisonaitken", okertype = "liracine")
    if (cfg$type == "lp") args <- c(args, list(degree = rep(cfg$degree, 2L),
      basis = cfg$basis, bernstein.basis = cfg$bernstein))
    bw <- do.call(npregbw, args)
    s <- c88_rows(x, x, k, cfg$degree, cfg$basis, TRUE)
    residual <- y - as.vector(s %*% y)
    variance <- residual^2 / rowSums((diag(n) - s)^2)
    for (identity in c(TRUE, FALSE)) {
      z <- if (identity) x else x[c(2L, 13L, 28L, 51L, 70L), ]
      fit.args <- list(bws = bw, txdat = x, tydat = y,
                       gradients = TRUE, se = TRUE, warn.glp.gradient = FALSE)
      if (!identity) fit.args$exdat <- z
      fit <- do.call(npreg, fit.args)
      fit.args$se <- FALSE
      off <- do.call(npreg, fit.args)
      expect_identical(fit$mean, off$mean)
      expect_identical(fit$grad, off$grad)
      a <- if (identity) s else c88_rows(x, z, k, cfg$degree, cfg$basis, FALSE)
      label <- paste(cfg$type, cfg$degree, cfg$basis, cfg$bernstein, identity)
      expect_equal(fit$mean, as.vector(a %*% y), tolerance = 2e-9, info = label)
      expect_equal(fit$merr, sqrt(as.vector(a^2 %*% variance)),
                   tolerance = 2e-9, info = label)
      for (column in 3:4) {
        endpoints <- c88_endpoints(z, column)
        delta <- c88_rows(x, endpoints$hi, k, cfg$degree, cfg$basis, identity) -
          c88_rows(x, endpoints$lo, k, cfg$degree, cfg$basis, identity)
        expect_equal(fit$grad[, column], as.vector(delta %*% y),
                     tolerance = 2e-9, info = label)
        expect_equal(fit$gerr[, column], sqrt(as.vector(delta^2 %*% variance)),
                     tolerance = 2e-9, info = label)
      }
    }
  }
})

test_that("GNN training categorical effects retain radius identity at the lower bound", {
  skip_on_cran()
  spawn_mpi_slaves()
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 36L),
                  u = factor(rep(c("a", "b"), 18L)))
  x$x[2L] <- x$x[1L]
  y <- sin(x$x) + .1 * (x$u == "b")
  for (type in c("lc", "ll", "lp")) {
    degree <- if (type == "lp") 1L else NULL
    expect_error(npregbw(xdat = x, ydat = y, bws = c(1, .3),
      bandwidth.compute = FALSE, bwtype = "generalized_nn", regtype = type,
      degree = degree),
      "nearest-neighbor bandwidth must be in")
    bw <- npregbw(xdat = x, ydat = y, bws = c(2, .3),
      bandwidth.compute = FALSE, bwtype = "generalized_nn", regtype = type,
      degree = degree)
    fit <- npreg(bw, gradients = TRUE)
    expect_true(all(is.finite(fit$grad[, 2L])))
    expect_error(npreg(bw, exdat = x, gradients = TRUE), "zero.*radius|radius.*zero")
    above <- npregbw(xdat = x, ydat = y, bws = c(3, .3),
      bandwidth.compute = FALSE, bwtype = "generalized_nn", regtype = type,
      degree = degree)
    expect_true(all(is.finite(npreg(above, exdat = x, gradients = TRUE)$grad[, 2L])))
    tied <- x
    tied$x[3L] <- tied$x[1L]
    expect_error(npreg(bw, txdat = tied, tydat = y, gradients = TRUE),
                  "zero.*radius|radius.*zero")
  }
})
