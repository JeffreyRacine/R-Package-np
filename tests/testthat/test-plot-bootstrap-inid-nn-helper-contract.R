test_that("nearest-neighbor regression helper matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(32905)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.6 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(15), , drop = FALSE]
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))
  fast.fun <- getFromNamespace(".np_inid_boot_from_reghat_exact", "np")

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL),
    list(regtype = "ll", basis = NULL, degree = NULL),
    list(regtype = "lp", basis = "glp", degree = c(2L, 2L))
  )

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) c(5, 5) else c(2, 2)

    for (cfg in cfgs) {
      bw.args <- list(
        xdat = tx,
        ydat = y,
        regtype = cfg$regtype,
        bws = bw.val,
        bwtype = bt,
        bandwidth.compute = FALSE
      )
      if (!is.null(cfg$basis)) {
        bw.args$basis <- cfg$basis
        bw.args$degree <- cfg$degree
      }
      bw <- do.call(npregbw, bw.args)

      fast.out <- fast.fun(
        xdat = tx,
        exdat = ex,
        bws = bw,
        ydat = y,
        B = B,
        counts = counts
      )

      explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
      for (b in seq_len(B)) {
        idx <- rep.int(seq_len(n), counts[, b])
        explicit.t[b, ] <- npreg(
          txdat = tx[idx, , drop = FALSE],
          tydat = y[idx],
          exdat = ex,
          bws = bw,
          gradients = FALSE,
          warn.glp.gradient = FALSE
        )$mean
      }

      expect_equal(
        fast.out$t,
        explicit.t,
        tolerance = 1e-10,
        info = paste(bt, cfg$regtype)
      )
      expect_equal(
        fast.out$t0,
        npreg(txdat = tx, tydat = y, exdat = ex, bws = bw, gradients = FALSE, warn.glp.gradient = FALSE)$mean,
        tolerance = 1e-10,
        info = paste(bt, cfg$regtype)
      )
    }

    bw.grad <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      bws = bw.val,
      bwtype = bt,
      bandwidth.compute = FALSE
    )
    fast.grad <- fast.fun(
      xdat = tx,
      exdat = ex,
      bws = bw.grad,
      ydat = y,
      B = B,
      counts = counts,
      gradients = TRUE,
      gradient.order = 1L,
      slice.index = 1L
    )

    explicit.grad <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      explicit.grad[b, ] <- npreg(
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = ex,
        bws = bw.grad,
        gradients = TRUE,
        gradient.order = 1L,
        warn.glp.gradient = FALSE
      )$grad[, 1L]
    }

    expect_equal(
      fast.grad$t,
      explicit.grad,
      tolerance = 1e-10,
      info = paste(bt, "ll-grad")
    )
    expect_equal(
      fast.grad$t0,
      npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.grad, gradients = TRUE, gradient.order = 1L, warn.glp.gradient = FALSE)$grad[, 1L],
      tolerance = 1e-10,
      info = paste(bt, "ll-grad")
    )
  }
})

test_that("nearest-neighbor ksum helpers match explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(32906)
  n <- 40
  tx <- data.frame(x = runif(n))
  ty <- data.frame(y = runif(n))
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 11))
  ey <- data.frame(y = seq(0.1, 0.9, length.out = 11))
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.u <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "np")
  fast.c <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) 5 else 2
    cbw.val <- rep.int(bw.val, 2L)

    u.dens.bw <- npudensbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    u.dist.bw <- npudistbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dens.bw <- npcdensbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dist.bw <- npcdistbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)

    u.dens.fast <- fast.u(
      xdat = tx,
      exdat = ex,
      bws = u.dens.bw,
      B = B,
      operator = "normal",
      counts = counts
    )
    u.dist.fast <- fast.u(
      xdat = tx,
      exdat = ex,
      bws = u.dist.bw,
      B = B,
      operator = "integral",
      counts = counts
    )
    c.dens.fast <- fast.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dens.bw,
      B = B,
      cdf = FALSE,
      counts = counts
    )
    c.dist.fast <- fast.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dist.bw,
      B = B,
      cdf = TRUE,
      counts = counts
    )

    u.dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    u.dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    c.dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    c.dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      u.dens.explicit[b, ] <- npudens(tdat = tx[idx, , drop = FALSE], edat = ex, bws = u.dens.bw)$dens
      u.dist.explicit[b, ] <- npudist(tdat = tx[idx, , drop = FALSE], edat = ex, bws = u.dist.bw)$dist
      c.dens.explicit[b, ] <- npcdens(
        txdat = tx[idx, , drop = FALSE],
        tydat = ty[idx, , drop = FALSE],
        exdat = ex,
        eydat = ey,
        bws = c.dens.bw
      )$condens
      c.dist.explicit[b, ] <- npcdist(
        txdat = tx[idx, , drop = FALSE],
        tydat = ty[idx, , drop = FALSE],
        exdat = ex,
        eydat = ey,
        bws = c.dist.bw
      )$condist
    }

    expect_equal(u.dens.fast$t, u.dens.explicit, tolerance = 1e-10, info = paste(bt, "udens"))
    expect_equal(u.dist.fast$t, u.dist.explicit, tolerance = 1e-10, info = paste(bt, "udist"))
    expect_equal(c.dens.fast$t, c.dens.explicit, tolerance = 1e-10, info = paste(bt, "npcdens"))
    expect_equal(c.dist.fast$t, c.dist.explicit, tolerance = 1e-10, info = paste(bt, "npcdist"))
    expect_equal(u.dens.fast$t0, npudens(tdat = tx, edat = ex, bws = u.dens.bw)$dens, tolerance = 1e-12, info = paste(bt, "udens-t0"))
    expect_equal(u.dist.fast$t0, npudist(tdat = tx, edat = ex, bws = u.dist.bw)$dist, tolerance = 1e-12, info = paste(bt, "udist-t0"))
    expect_equal(c.dens.fast$t0, npcdens(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = c.dens.bw)$condens, tolerance = 1e-12, info = paste(bt, "npcdens-t0"))
    expect_equal(c.dist.fast$t0, npcdist(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = c.dist.bw)$condist, tolerance = 1e-12, info = paste(bt, "npcdist-t0"))
  }
})

test_that("nearest-neighbor frozen hat helpers match frozen operator application", {
  skip_if_not_installed("np")

  set.seed(32907)
  n <- 36
  tx <- data.frame(x = runif(n))
  ty <- data.frame(y = runif(n))
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 9))
  ey <- data.frame(y = seq(0.15, 0.85, length.out = 9))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  frozen.u <- getFromNamespace(".np_inid_boot_from_hat_unconditional_frozen", "np")
  frozen.c <- getFromNamespace(".np_inid_boot_from_hat_conditional_frozen", "np")
  udens.hat <- getFromNamespace("npudenshat", "np")
  udist.hat <- getFromNamespace("npudisthat", "np")
  cdens.hat <- getFromNamespace("npcdenshat", "np")
  cdist.hat <- getFromNamespace("npcdisthat", "np")
  make_xkbw <- getFromNamespace(".npcdhat_make_xkbw", "np")
  make_ybw <- getFromNamespace(".npcdhat_make_ybw", "np")
  make_kmat <- getFromNamespace(".npcdhat_make_kernel_matrix", "np")

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

    expect_equal(u.dens.frozen$t, t(H.udens %*% counts), tolerance = 1e-10, info = paste(bt, "udens"))
    expect_equal(u.dist.frozen$t, t(H.udist %*% counts), tolerance = 1e-10, info = paste(bt, "udist"))
    expect_equal(
      c.dens.frozen$t,
      t(num.op.dens %*% counts) / pmax(t(den.op.dens %*% counts), .Machine$double.eps),
      tolerance = 1e-10,
      info = paste(bt, "npcdens")
    )
    expect_equal(
      c.dist.frozen$t,
      t(num.op.dist %*% counts) / pmax(t(den.op.dist %*% counts), .Machine$double.eps),
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

test_that("npindex nearest-neighbor inid helper matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(3292)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))
  fast.fun <- getFromNamespace(".np_inid_boot_from_index", "np")

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL, h = 5L),
    list(regtype = "ll", basis = NULL, degree = NULL, h = 5L),
    list(regtype = "lp", basis = "tensor", degree = 2L, h = 5L)
  )

  for (bt in c("generalized_nn", "adaptive_nn")) {
    for (cfg in cfgs) {
      bw.args <- list(
        xdat = tx,
        ydat = y,
        bws = c(1, 1, cfg$h),
        bandwidth.compute = FALSE,
        regtype = cfg$regtype,
        bwtype = bt
      )
      if (!is.null(cfg$basis)) {
        bw.args$basis <- cfg$basis
        bw.args$degree <- cfg$degree
      }
      bw <- do.call(npindexbw, bw.args)

      fast.out <- fast.fun(
        xdat = tx,
        ydat = y,
        bws = bw,
        B = B,
        counts = counts,
        gradients = FALSE
      )

      explicit.t <- matrix(NA_real_, nrow = B, ncol = n)
      for (b in seq_len(B)) {
        idx <- rep.int(seq_len(n), counts[, b])
        explicit.t[b, ] <- npindex(
          txdat = tx[idx, , drop = FALSE],
          tydat = y[idx],
          exdat = tx,
          bws = bw,
          gradients = FALSE
        )$mean
      }

      expect_equal(
        fast.out$t,
        explicit.t,
        tolerance = 1e-8,
        info = paste(bt, cfg$regtype)
      )
      expect_equal(
        fast.out$t0,
        npindex(txdat = tx, tydat = y, exdat = tx, bws = bw, gradients = FALSE)$mean,
        tolerance = 1e-8,
        info = paste(bt, cfg$regtype)
      )
    }
  }
})
