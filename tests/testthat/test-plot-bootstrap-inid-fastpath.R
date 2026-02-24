test_that("inid lc fast path matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(321)
  n <- 40
  x <- runif(n)
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 21))

  bw <- npregbw(
    y ~ x,
    regtype = "lc",
    bws = c(0.2),
    bandwidth.compute = FALSE
  )

  H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  B <- 13L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  fast.out <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    fit.b <- npreg(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = ex,
      bws = bw,
      gradients = FALSE,
      warn.glp.gradient = FALSE
    )
    explicit.t[b, ] <- fit.b$mean
  }

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-10)
  expect_equal(fast.out$t0, as.vector(H %*% y), tolerance = 1e-12)
})

test_that("inid lc fast path toggle preserves plot bootstrap contract", {
  skip_if_not_installed("np")

  set.seed(322)
  n <- 60
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.15)
  bw <- npregbw(y ~ x, regtype = "lc", nmulti = 1)

  run_plot <- function(disable) {
    old <- getOption("np.plot.inid.fastpath.disable")
    on.exit(options(np.plot.inid.fastpath.disable = old), add = TRUE)
    options(np.plot.inid.fastpath.disable = disable)

    suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 9
      )
    )
  }

  out.fast <- run_plot(disable = FALSE)
  out.legacy <- run_plot(disable = TRUE)

  expect_type(out.fast, "list")
  expect_true(length(out.fast) > 0)
  expect_type(out.legacy, "list")
  expect_true(length(out.legacy) > 0)
})

test_that("plot bootstrap accepts wild alias", {
  skip_if_not_installed("np")

  set.seed(323)
  n <- 80
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  bw <- npregbw(y ~ x, bws = 0.2, bandwidth.compute = FALSE)

  out <- suppressWarnings(
    plot(
      bw,
      xdat = data.frame(x = x),
      ydat = y,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 9
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("inid lc chunked generation matches explicit chunked draws", {
  skip_if_not_installed("np")

  set.seed(324)
  n <- 30
  x <- runif(n)
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 17))
  bw <- npregbw(y ~ x, regtype = "lc", bws = c(0.2), bandwidth.compute = FALSE)
  H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  B <- 11L

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  old.chunk <- getOption("np.plot.inid.chunk.size")
  on.exit(options(np.plot.inid.chunk.size = old.chunk), add = TRUE)
  options(np.plot.inid.chunk.size = 3L)

  set.seed(7001)
  auto <- fast.fun(H = H, ydat = y, B = B, counts = NULL)

  set.seed(7001)
  counts <- matrix(NA_integer_, nrow = n, ncol = B)
  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + 3L - 1L)
    bsz <- stopi - start + 1L
    counts[, start:stopi] <- rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    start <- stopi + 1L
  }
  explicit <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  expect_equal(auto$t, explicit$t, tolerance = 1e-12)
  expect_equal(auto$t0, explicit$t0, tolerance = 1e-12)
})

test_that("inid ksum fast path matches explicit resample refits for npudens/npudist", {
  skip_if_not_installed("np")

  set.seed(325)
  n <- 70
  x <- rnorm(n)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 19))
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  dens.bw <- npudensbw(dat = tx, bws = 0.35, bandwidth.compute = FALSE)
  dist.bw <- npudistbw(dat = tx, bws = 0.35, bandwidth.compute = FALSE)
  fast.fun <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "np")

  dens.fast <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = dens.bw,
    B = B,
    operator = "normal",
    counts = counts
  )
  dist.fast <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = dist.bw,
    B = B,
    operator = "integral",
    counts = counts
  )

  dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    dens.explicit[b, ] <- npudens(tdat = tx[idx, , drop = FALSE], edat = ex, bws = dens.bw)$dens
    dist.explicit[b, ] <- npudist(tdat = tx[idx, , drop = FALSE], edat = ex, bws = dist.bw)$dist
  }

  expect_equal(dens.fast$t, dens.explicit, tolerance = 1e-10)
  expect_equal(dist.fast$t, dist.explicit, tolerance = 1e-10)
  expect_equal(dens.fast$t0, npudens(tdat = tx, edat = ex, bws = dens.bw)$dens, tolerance = 1e-12)
  expect_equal(dist.fast$t0, npudist(tdat = tx, edat = ex, bws = dist.bw)$dist, tolerance = 1e-12)
})

test_that("inid ksum fast path matches explicit resample refits for npcdens/npcdist", {
  skip_if_not_installed("np")

  set.seed(326)
  n <- 80
  x <- rnorm(n)
  y <- rnorm(n)
  tx <- data.frame(x = x)
  ty <- data.frame(y = y)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 17))
  ey <- data.frame(y = seq(min(y), max(y), length.out = 17))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  dens.bw <- npcdensbw(xdat = tx, ydat = ty, bws = c(0.45, 0.45), bandwidth.compute = FALSE)
  dist.bw <- npcdistbw(xdat = tx, ydat = ty, bws = c(0.45, 0.45), bandwidth.compute = FALSE)
  fast.fun <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")

  dens.fast <- fast.fun(
    xdat = tx,
    ydat = ty,
    exdat = ex,
    eydat = ey,
    bws = dens.bw,
    B = B,
    cdf = FALSE,
    counts = counts
  )
  dist.fast <- fast.fun(
    xdat = tx,
    ydat = ty,
    exdat = ex,
    eydat = ey,
    bws = dist.bw,
    B = B,
    cdf = TRUE,
    counts = counts
  )

  expect_true(is.list(dens.fast))
  expect_true(is.list(dist.fast))

  dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    dens.explicit[b, ] <- npcdens(
      txdat = tx[idx, , drop = FALSE],
      tydat = ty[idx, , drop = FALSE],
      exdat = ex,
      eydat = ey,
      bws = dens.bw
    )$condens
    dist.explicit[b, ] <- npcdist(
      txdat = tx[idx, , drop = FALSE],
      tydat = ty[idx, , drop = FALSE],
      exdat = ex,
      eydat = ey,
      bws = dist.bw
    )$condist
  }

  expect_equal(dens.fast$t, dens.explicit, tolerance = 1e-10)
  expect_equal(dist.fast$t, dist.explicit, tolerance = 1e-10)
  expect_equal(dens.fast$t0, npcdens(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = dens.bw)$condens, tolerance = 1e-12)
  expect_equal(dist.fast$t0, npcdist(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = dist.bw)$condist, tolerance = 1e-12)
})

test_that("npreg plot bootstrap inid supports lp basis variants", {
  skip_if_not_installed("np")

  set.seed(327)
  n <- 70
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)

  for (basis in c("additive", "tensor", "glp")) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = 2,
      basis = basis,
      bws = 0.2,
      bandwidth.compute = FALSE
    )

    out <- suppressWarnings(
      plot(
        bw,
        xdat = tx,
        ydat = y,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 7
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0)
  }
})
