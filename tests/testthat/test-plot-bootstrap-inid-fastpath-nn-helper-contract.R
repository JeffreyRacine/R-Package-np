test_that("nearest-neighbor regression and ksum helpers match explicit resample refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(32905)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.5 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(12), , drop = FALSE]
  ty <- data.frame(y = runif(n))
  ey <- data.frame(y = seq(0.1, 0.9, length.out = 12))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  reg.fast <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")
  fast.u <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "npRmpi")
  fast.c <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "npRmpi")

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) c(5, 5) else c(2, 2)
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

    ubw <- npudensbw(dat = data.frame(x = x1), bws = bw.val[1L], bwtype = bt, bandwidth.compute = FALSE)
    uout <- fast.u(
      xdat = data.frame(x = x1),
      exdat = data.frame(x = ex$x1),
      bws = ubw,
      B = B,
      operator = "normal",
      counts = counts
    )
    uexplicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      uexplicit[b, ] <- npudens(
        tdat = data.frame(x = x1)[idx, , drop = FALSE],
        edat = data.frame(x = ex$x1),
        bws = ubw
      )$dens
    }
    expect_equal(uout$t, uexplicit, tolerance = 1e-10, info = paste(bt, "npudens"))

    cbw <- npcdensbw(
      xdat = data.frame(x = x1),
      ydat = ty,
      bws = rep.int(bw.val[1L], 2L),
      bwtype = bt,
      bandwidth.compute = FALSE
    )
    cout <- fast.c(
      xdat = data.frame(x = x1),
      ydat = ty,
      exdat = data.frame(x = ex$x1),
      eydat = ey,
      bws = cbw,
      B = B,
      cdf = FALSE,
      counts = counts
    )
    cexplicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      cexplicit[b, ] <- npcdens(
        txdat = data.frame(x = x1)[idx, , drop = FALSE],
        tydat = ty[idx, , drop = FALSE],
        exdat = data.frame(x = ex$x1),
        eydat = ey,
        bws = cbw
      )$condens
    }
    expect_equal(cout$t, cexplicit, tolerance = 1e-10, info = paste(bt, "npcdens"))
  }
})

test_that("nearest-neighbor frozen hat helpers match frozen operator application", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(32907)
  n <- 36
  tx <- data.frame(x = runif(n))
  ty <- data.frame(y = runif(n))
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 9))
  ey <- data.frame(y = seq(0.15, 0.85, length.out = 9))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  frozen.u <- getFromNamespace(".np_inid_boot_from_hat_unconditional_frozen", "npRmpi")
  frozen.c <- getFromNamespace(".np_inid_boot_from_hat_conditional_frozen", "npRmpi")
  udens.hat <- getFromNamespace("npudenshat", "npRmpi")
  udist.hat <- getFromNamespace("npudisthat", "npRmpi")
  cdens.hat <- getFromNamespace("npcdenshat", "npRmpi")
  cdist.hat <- getFromNamespace("npcdisthat", "npRmpi")
  make_xkbw <- getFromNamespace(".npcdhat_make_xkbw", "npRmpi")
  make_ybw <- getFromNamespace(".npcdhat_make_ybw", "npRmpi")
  make_kmat <- getFromNamespace(".npcdhat_make_kernel_matrix", "npRmpi")

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) 5 else 2
    cbw.val <- rep.int(bw.val, 2L)

    u.dens.bw <- npudensbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    u.dist.bw <- npudistbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dens.bw <- npcdensbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dist.bw <- npcdistbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)

    H.udens <- unclass(udens.hat(bws = u.dens.bw, tdat = tx, edat = ex, output = "matrix"))
    H.udist <- unclass(udist.hat(bws = u.dist.bw, tdat = tx, edat = ex, output = "matrix"))
    H.cdens <- unclass(cdens.hat(bws = c.dens.bw, txdat = tx, tydat = ty, exdat = ex, eydat = ey, output = "matrix"))
    H.cdist <- unclass(cdist.hat(bws = c.dist.bw, txdat = tx, tydat = ty, exdat = ex, eydat = ey, output = "matrix"))
    xkbw.dens <- make_xkbw(bws = c.dens.bw, txdat = tx)
    xkbw.dist <- make_xkbw(bws = c.dist.bw, txdat = tx)
    ybw.dens <- make_ybw(bws = c.dens.bw, tydat = ty)
    ybw.dist <- make_ybw(bws = c.dist.bw, tydat = ty)
    Kx.dens <- make_kmat(
      kbw = xkbw.dens,
      txdat = tx,
      exdat = ex,
      operator = rep.int("normal", ncol(tx))
    )
    Ky.dens <- make_kmat(
      kbw = ybw.dens,
      txdat = ty,
      exdat = ey,
      operator = rep.int("normal", ncol(ty))
    )
    Kx.dist <- make_kmat(
      kbw = xkbw.dist,
      txdat = tx,
      exdat = ex,
      operator = rep.int("normal", ncol(tx))
    )
    Ky.dist <- make_kmat(
      kbw = ybw.dist,
      txdat = ty,
      exdat = ey,
      operator = rep.int("integral", ncol(ty))
    )
    den.op.dens <- Kx.dens / n
    num.op.dens <- (Kx.dens * Ky.dens) / n
    den.op.dist <- Kx.dist / n
    num.op.dist <- (Kx.dist * Ky.dist) / n

    u.dens.frozen <- frozen.u(
      xdat = tx,
      exdat = ex,
      bws = u.dens.bw,
      B = B,
      operator = "normal",
      counts = counts
    )
    u.dist.frozen <- frozen.u(
      xdat = tx,
      exdat = ex,
      bws = u.dist.bw,
      B = B,
      operator = "integral",
      counts = counts
    )
    c.dens.frozen <- frozen.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dens.bw,
      B = B,
      cdf = FALSE,
      counts = counts
    )
    c.dist.frozen <- frozen.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dist.bw,
      B = B,
      cdf = TRUE,
      counts = counts
    )

    expect_equal(u.dens.frozen$t, crossprod(counts, t(H.udens)), tolerance = 1e-10, info = paste(bt, "udens"))
    expect_equal(u.dist.frozen$t, crossprod(counts, t(H.udist)), tolerance = 1e-10, info = paste(bt, "udist"))
    expect_equal(
      c.dens.frozen$t,
      crossprod(counts, t(num.op.dens)) / pmax(crossprod(counts, t(den.op.dens)), .Machine$double.eps),
      tolerance = 1e-10,
      info = paste(bt, "npcdens")
    )
    expect_equal(
      c.dist.frozen$t,
      crossprod(counts, t(num.op.dist)) / pmax(crossprod(counts, t(den.op.dist)), .Machine$double.eps),
      tolerance = 1e-10,
      info = paste(bt, "npcdist")
    )

    expect_equal(u.dens.frozen$t0, rowSums(H.udens), tolerance = 1e-12, info = paste(bt, "udens-t0"))
    expect_equal(u.dist.frozen$t0, rowSums(H.udist), tolerance = 1e-12, info = paste(bt, "udist-t0"))
    expect_equal(
      c.dens.frozen$t0,
      rowSums(num.op.dens) / pmax(rowSums(den.op.dens), .Machine$double.eps),
      tolerance = 1e-12,
      info = paste(bt, "npcdens-t0")
    )
    expect_equal(
      c.dist.frozen$t0,
      rowSums(num.op.dist) / pmax(rowSums(den.op.dist), .Machine$double.eps),
      tolerance = 1e-12,
      info = paste(bt, "npcdist-t0")
    )
  }
})
