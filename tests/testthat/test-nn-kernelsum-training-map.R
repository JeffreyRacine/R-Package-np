test_that("mapped generalized-NN sums preserve training occurrence geometry", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.extendednn = TRUE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(17L)
  z <- data.frame(z = sin(i * sqrt(2)) + i / 17)
  tensor <- cbind(1, z$z, cos(i))
  idx <- c(17L, 2L, 9L, 2L)
  evaluate <- function(args) {
    args$bws <- kbandwidth(bw = args$bws, nobs = nrow(args$txdat),
      xdati = untangle(args$txdat), xnames = names(args$txdat),
      bwtype = args$bwtype)
    args$bwtype <- NULL
    .npRmpi_with_local_regression(do.call(npksum, args))
  }
  for(tree in c(FALSE, TRUE)) for(k in c(1L, 8L, 16L, 32L)) {
    options(np.tree = tree)
    for(divide in c(FALSE, TRUE)) {
      args <- list(txdat = z, bws = k, bwtype = "generalized_nn",
        tydat = tensor, weights = tensor, bandwidth.divide = divide,
        return.kernel.weights = TRUE,
        .np.internal.bandwidth.divide.weights = divide)
      direct <- evaluate(c(args, list(leave.one.out = TRUE)))
      mapped <- evaluate(c(args, list(exdat = z[idx, , drop = FALSE],
        .np.internal.eval.train.index = idx)))
      radius <- vapply(idx, function(j)
        sort(abs(z$z[-j] - z$z[j]))[min(k, 16L)] * max(1, k / 16), numeric(1))
      literal <- vapply(seq_along(idx), function(j) {
        w <- dnorm((z$z - z$z[idx[j]]) / radius[j])
        if(divide) w <- w / radius[j]
        w
      }, numeric(17L))
      expect_equal(mapped$kw, literal, tolerance = 2e-13)
      for(j in seq_along(idx))
        mapped$ksum[, , j] <- mapped$ksum[, , j] -
          mapped$kw[idx[j], j] * tcrossprod(tensor[idx[j], ])
      expect_equal(mapped$ksum, direct$ksum[, , idx, drop = FALSE],
                   tolerance = 2e-12)
    }
  }
  options(np.tree = FALSE)
  base <- list(txdat = z, exdat = z[idx, , drop = FALSE], bws = 8,
    bwtype = "generalized_nn")
  for(bad in list(as.double(idx), c(0L, idx[-1L]), c(NA_integer_, idx[-1L]),
                  idx[-1L]))
    expect_error(evaluate(c(base, list(.np.internal.eval.train.index = bad))),
                 "internal mapped NN", fixed = TRUE)
  tied <- data.frame(z = rep(0:2, each = 6L))
  expect_error(evaluate(list(txdat = tied, exdat = tied[1L, , drop = FALSE],
    bws = 1, bwtype = "generalized_nn", .np.internal.eval.train.index = 1L)),
    "zero", fixed = TRUE)
  expect_true(all(is.finite(evaluate(c(base, list(
    .np.internal.eval.train.index = idx)))$ksum)))

  # Protect the consumer's ANN normalization as well as the GNN mapping.
  check.objective <- function(type, degree) {
    i <- seq_len(24L)
    x <- data.frame(w = cos(i / 5) + i / 24)
    z <- data.frame(z = sin(i * sqrt(2)) + i / 24)
    y <- (1 + z$z / 4) * x$w + sin(i * sqrt(5)) / 20
    bw <- .npscoefbw_build_scbandwidth(x, y, z, 12, FALSE,
      list(regtype = "lp", degree = degree, bwtype = type))
    bw <- .npscoefbw_normalize_nomad_scbw(bw, z, 12)
    ctx <- .npscoefbw_nomad_context_prepare(x, y, z)
    pool <- NULL
    on.exit({
      if(!is.null(pool)) .npscoefbw_nomad_pool_stop(pool)
      .npscoefbw_nomad_context_cleanup(ctx)
    })
    pool <- .npscoefbw_nomad_pool_start(ctx)
    direct <- .npscoefbw_nomad_eval_direct(ctx, bw)
    pooled <- .npscoefbw_eval_pool(ctx, bw, pool)
    expect_lt(direct$objective, 1)
    expect_equal(pooled$objective, direct$objective, tolerance = 1e-11)
  }
  for(type in c("generalized_nn", "adaptive_nn")) for(degree in 0:1)
    check.objective(type, degree)
})
