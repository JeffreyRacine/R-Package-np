test_that("single-index pairs SDs match an explicit Gaussian resample oracle", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 24L
  for (p in 1:2) {
    x <- data.frame(x1 = seq(-1, 1, length.out = n))
    if (p == 2L) x$x2 <- cos(seq_len(n))
    beta <- c(1, .6)[seq_len(p)]
    index <- as.vector(as.matrix(x) %*% beta)
    y <- sin(index) + .1*cos(seq_len(n))
    bw <- suppressWarnings(npindexbw(xdat = x, ydat = y,
      bws = c(beta, .7), method = "ichimura", regtype = "lc",
      bandwidth.compute = FALSE))
    set.seed(619L)
    plan <- boot::boot(data.frame(x, y), function(data, indices) 0, R = 3L)
    indices <- boot::boot.array(plan, indices = TRUE)
    expected.rng <- .Random.seed
    means <- vapply(seq_len(3L), function(b) {
      take <- indices[b, ]
      projected <- as.vector(as.matrix(x[take, , drop = FALSE]) %*% beta)
      weight <- dnorm(outer(index, projected, "-")/.7)
      as.vector(weight %*% y[take])/rowSums(weight)
    }, numeric(n))
    set.seed(619L)
    fit <- npindex(bws = bw, txdat = x, tydat = y, se = TRUE, se.type = "bootstrap", B = 3L)
    expect_equal(as.vector(se(fit)), apply(means, 1L, sd), tolerance = 1e-12)
    expect_identical(.Random.seed, expected.rng)
    expect_true(all(is.finite(se(fit))))
  }
})
