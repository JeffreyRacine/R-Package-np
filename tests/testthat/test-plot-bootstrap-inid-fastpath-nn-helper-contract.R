test_that("generalized nearest-neighbor regression helper matches explicit refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(32905)
  n <- 12
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.5 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(3), , drop = FALSE]
  B <- 1L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  reg.fast <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")

  for (bt in "generalized_nn") {
    bw.val <- c(2, 2)
    rbw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      bws = bw.val,
      bwtype = bt,
      bandwidth.compute = FALSE
    )
    reg.out <- reg.fast(
      xdat = tx,
      exdat = ex,
      bws = rbw,
      ydat = y,
      B = B,
      counts = counts
    )
    reg.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      reg.explicit[b, ] <- npreg(
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = ex,
        bws = rbw,
        gradients = FALSE,
        warn.glp.gradient = FALSE
      )$mean
    }
    expect_equal(reg.out$t, reg.explicit, tolerance = 1e-10, info = paste(bt, "npreg"))
    expect_equal(
      reg.out$t0,
      npreg(txdat = tx, tydat = y, exdat = ex, bws = rbw, gradients = FALSE, warn.glp.gradient = FALSE)$mean,
      tolerance = 1e-10,
      info = paste(bt, "npreg-t0")
    )
  }
})
