test_that("NN copula errors combine aligned numerator and marginal scores", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  owner <- get(".npcopula_density_se_nn", ns)
  x <- data.frame(x = sin(1:24) + (1:24)/24,
                  y = cos(1:24) + (1:24)/48)
  ex <- data.frame(x = c(.15, .65), y = c(.2, .6))
  n <- nrow(x); k <- 6L
  center <- function(z) z - mean(z)
  for (type in c("generalized_nn", "adaptive_nn")) {
    bw <- npudensbw(dat = x, bws = c(k, k), bandwidth.compute = FALSE,
                    bwtype = type)
    ann <- identical(type, "adaptive_nn")
    radius <- pilot <- matrix(NA_real_, n, 2L)
    membership <- vector("list", 2L)
    if (ann) for (j in 1:2) {
      distance <- abs(outer(x[[j]], x[[j]], "-"))
      width <- min(floor((n - 1)^(2/3)), floor(min(k, n - 1 - k)/2))
      for (i in seq_len(n)) {
        d <- sort(distance[i, -i])
        radius[i, j] <- d[k]
        pilot[i, j] <- (2 * width/(n - 1))/(d[k + width] - d[k - width])
      }
      membership[[j]] <- distance <= radius[, j]
    }
    reference <- vapply(seq_len(nrow(ex)), function(q) {
      h <- if (ann) radius else matrix(rep(vapply(1:2, function(j)
        sort(abs(x[[j]] - ex[q, j]))[k], numeric(1L)), each = n), n, 2L)
      z <- sweep(as.matrix(x), 2L, as.numeric(ex[q, ]), function(a, b) b - a)/h
      w <- dnorm(z)/h
      M <- colMeans(w); D <- prod(M); K <- w[, 1L] * w[, 2L]
      J <- mean(K); phi <- K; marginal <- w
      if (ann) for (j in 1:2) {
        phi <- phi + as.numeric(crossprod(
          K * (1 - z[, j]^2)/(h[, j] * pilot[, j]), membership[[j]]))/n
        marginal[, j] <- marginal[, j] + as.numeric(crossprod(
          w[, j] * (1 - z[, j]^2)/(h[, j] * pilot[, j]), membership[[j]]))/n
      }
      score <- center(phi)/D - (J/D) *
        (center(marginal[, 1L])/M[1L] + center(marginal[, 2L])/M[2L])
      sqrt(sum(center(score)^2)) * if (ann) 1/sqrt(n * (n - 1)) else 1/n
    }, numeric(1L))
    expect_equal(owner(bw, x, ex), reference, tolerance = 2e-10)
    single <- data.frame(x = ordered("a"))
    single.bw <- npudensbw(dat = single, bws = 0, bwtype = type,
                           bandwidth.compute = FALSE)
    expect_identical(owner(single.bw, single, single), 0)
    b1 <- npudensbw(dat = x[1L], bws = k, bandwidth.compute = FALSE, bwtype = type)
    expect_identical(owner(b1, x[1L], ex[1L]), c(0, 0))
    off <- npcopula(bws = bw, data = x)
    on <- npcopula(bws = bw, data = x, se = TRUE)
    expect_identical(as.data.frame(off), as.data.frame(on))
    expect_identical(off$bws, on$bws)
    expect_null(off$copulaerr)
    expect_identical(predict(on, se.fit = TRUE)$se.fit, se(on))
  }
})

test_that("NN ratio norm is stable without erasing small positive uncertainty", {
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  norm <- get(".npcopula_density_se_nn_norm", ns)
  expect_equal(norm(c(-1e200, 1e200)), sqrt(2) * 1e200)
  expect_equal(norm(c(-1e-200, 1e-200)), sqrt(2) * 1e-200)
  expect_gt(norm(c(-1e-200, 1e-200)), 0)
  expect_identical(norm(c(4, 4)), 0)
  expect_identical(norm(c(NA_real_, 1)), NA_real_)
})
