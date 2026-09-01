test_that("streamed IID tiles reproduce incumbent categorical and HC0 statistics", {
  set.seed(20260831)
  n <- 48L
  xdat <- data.frame(
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high", "top"), each = 12L),
                levels = c("low", "mid", "high", "top")),
    x = seq(-1.25, 1.25, length.out = n)
  )
  ydat <- 0.7 + c(a = 0, b = 0.4, c = -0.25)[xdat$u] +
    c(low = -0.2, mid = 0, high = 0.3, top = 0.55)[xdat$o] +
    0.8 * xdat$x + rnorm(n, sd = 0.15)
  bw <- npregbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.2, 0.2, 0.6), bandwidth.compute = FALSE,
    regtype = "ll", bwtype = "fixed"
  )
  unrestricted <- npreg(
    txdat = xdat, tydat = ydat, bws = bw,
    gradients = FALSE, se = FALSE
  )
  unrestricted.scaled <- scale(residuals(unrestricted))
  unrestricted.scale <- attr(unrestricted.scaled, "scaled:scale")
  unrestricted.center <- attr(unrestricted.scaled, "scaled:center")

  set.seed(771)
  donor <- vapply(
    seq_len(8L),
    function(unused) sample.int(n, replace = TRUE),
    integer(n)
  )

  for (tested.index in seq_len(ncol(xdat))) {
    null.frame <- xdat
    xq <- uocquantile(xdat[[tested.index]], 0.5)
    null.frame[[tested.index]] <- if (is.factor(xdat[[tested.index]])) {
      np:::cast(xq, xdat[[tested.index]], same.levels = TRUE)
    } else xq
    null.mean <- npreg(
      txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
      gradients = FALSE, se = FALSE
    )$mean
    residual.pool <- as.numeric(
      scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
    )
    residual.pool <- residual.pool - mean(residual.pool)

    streamed <- np:::.np_npsig_streamed_iid_tile(
      bws = bw,
      xdat = xdat,
      tested.index = tested.index,
      donor.index = donor,
      null.mean = null.mean,
      residual.pool = residual.pool
    )
    incumbent <- vapply(seq_len(ncol(donor)), function(replication) {
      ystar <- null.mean + residual.pool[donor[, replication]]
      fit <- npreg(
        txdat = xdat, tydat = ystar, bws = bw,
        gradients = TRUE,
        se = !(is.factor(xdat[[tested.index]]) ||
                 is.ordered(xdat[[tested.index]]))
      )
      np:::.np_npsig_statistic(
        fit,
        index = tested.index,
        pivot = !(is.factor(xdat[[tested.index]]) ||
                    is.ordered(xdat[[tested.index]]))
      )
    }, numeric(1L))

    expect_equal(streamed, incumbent, tolerance = 2e-12)
  }
})

test_that("streamed IID tiles reproduce mixed local-constant statistics", {
  set.seed(20260901)
  n <- 48L
  xdat <- data.frame(
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high", "top"), each = 12L),
                levels = c("low", "mid", "high", "top")),
    x1 = seq(-1.25, 1.25, length.out = n),
    x2 = sin(seq(-1.1, 1.3, length.out = n))
  )
  ydat <- 0.6 + c(a = 0, b = 0.35, c = -0.2)[xdat$u] +
    c(low = -0.25, mid = 0, high = 0.3, top = 0.5)[xdat$o] +
    0.8 * xdat$x1 - 0.35 * xdat$x2 + rnorm(n, sd = 0.14)
  bw <- npregbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.18, 0.16, 0.58, 0.64), bandwidth.compute = FALSE,
    regtype = "lc", bwtype = "fixed"
  )

  expect_true(np:::.np_npsig_streamed_iid_eligible(
    bw, xdat, seq_len(ncol(xdat)), FALSE, "I", "iid", NULL, list()
  ))

  unrestricted <- npreg(
    txdat = xdat, tydat = ydat, bws = bw,
    gradients = FALSE, se = FALSE
  )
  unrestricted.scaled <- scale(residuals(unrestricted))
  unrestricted.scale <- attr(unrestricted.scaled, "scaled:scale")
  unrestricted.center <- attr(unrestricted.scaled, "scaled:center")

  set.seed(771)
  donor <- vapply(
    seq_len(8L),
    function(unused) sample.int(n, replace = TRUE),
    integer(n)
  )

  for (tested.index in seq_len(ncol(xdat))) {
    null.frame <- xdat
    xq <- uocquantile(xdat[[tested.index]], 0.5)
    null.frame[[tested.index]] <- if (is.factor(xdat[[tested.index]])) {
      np:::cast(xq, xdat[[tested.index]], same.levels = TRUE)
    } else xq
    null.mean <- npreg(
      txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
      gradients = FALSE, se = FALSE
    )$mean
    residual.pool <- as.numeric(
      scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
    )
    residual.pool <- residual.pool - mean(residual.pool)

    streamed <- np:::.np_npsig_streamed_iid_tile(
      bws = bw,
      xdat = xdat,
      tested.index = tested.index,
      donor.index = donor,
      null.mean = null.mean,
      residual.pool = residual.pool
    )
    incumbent <- vapply(seq_len(ncol(donor)), function(replication) {
      ystar <- null.mean + residual.pool[donor[, replication]]
      categorical <- is.factor(xdat[[tested.index]]) ||
        is.ordered(xdat[[tested.index]])
      fit <- npreg(
        txdat = xdat, tydat = ystar, bws = bw,
        gradients = TRUE, se = !categorical
      )
      np:::.np_npsig_statistic(
        fit, index = tested.index, pivot = !categorical
      )
    }, numeric(1L))

    expect_equal(streamed, incumbent, tolerance = 2e-12)
  }
})

test_that("streamed IID capability is deterministic and semantics-exact", {
  xdat <- data.frame(
    u = factor(rep(c("a", "b"), 12L)),
    x = seq(-1, 1, length.out = 24L)
  )
  ydat <- seq_len(24L) / 24
  bw <- npregbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.2, 0.55), bandwidth.compute = FALSE,
    regtype = "ll", bwtype = "fixed"
  )
  eligible <- function(...) np:::.np_npsig_streamed_iid_eligible(
    bws = bw, xdat = xdat, index = 1:2, joint = FALSE,
    boot.type = "I", boot.method = "iid", pivot = NULL,
    extra.args = list(), ...
  )

  expect_true(eligible())
  expect_false(np:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, FALSE, "I", "wild", NULL, list()
  ))
  expect_false(np:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 2L, FALSE, "I", "iid", FALSE, list()
  ))
  expect_false(np:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, TRUE, "I", "iid", NULL, list()
  ))

  set.seed(31415)
  before <- .Random.seed
  result <- npsigtest(
    bw, xdat = xdat, ydat = ydat,
    B = 9L, random.seed = 81L
  )
  expect_identical(.Random.seed, before)
  expect_identical(dim(result$In.bootstrap), c(9L, 2L))
})

test_that("streamed IID categorical-only tiles reproduce scalar-LP fits", {
  set.seed(906)
  n <- 48L
  xdat <- data.frame(
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high", "top"), each = 12L),
                levels = c("low", "mid", "high", "top"))
  )
  ydat <- 0.5 + c(a = 0, b = 0.45, c = -0.3)[xdat$u] +
    c(low = -0.25, mid = 0, high = 0.35, top = 0.6)[xdat$o] +
    rnorm(n, sd = 0.12)
  bw <- npregbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.2, 0.2), bandwidth.compute = FALSE,
    regtype = "lc", bwtype = "fixed"
  )
  unrestricted <- npreg(
    txdat = xdat, tydat = ydat, bws = bw,
    gradients = FALSE, se = FALSE
  )
  unrestricted.scaled <- scale(residuals(unrestricted))
  unrestricted.scale <- attr(unrestricted.scaled, "scaled:scale")
  unrestricted.center <- attr(unrestricted.scaled, "scaled:center")

  set.seed(118)
  donor <- vapply(
    seq_len(8L),
    function(unused) sample.int(n, replace = TRUE),
    integer(n)
  )

  for (tested.index in seq_len(ncol(xdat))) {
    null.frame <- xdat
    xq <- uocquantile(xdat[[tested.index]], 0.5)
    null.frame[[tested.index]] <- np:::cast(
      xq, xdat[[tested.index]], same.levels = TRUE
    )
    null.mean <- npreg(
      txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
      gradients = FALSE, se = FALSE
    )$mean
    residual.pool <- as.numeric(
      scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
    )
    residual.pool <- residual.pool - mean(residual.pool)

    streamed <- np:::.np_npsig_streamed_iid_tile(
      bws = bw,
      xdat = xdat,
      tested.index = tested.index,
      donor.index = donor,
      null.mean = null.mean,
      residual.pool = residual.pool
    )
    incumbent <- vapply(seq_len(ncol(donor)), function(replication) {
      ystar <- null.mean + residual.pool[donor[, replication]]
      fit <- npreg(
        txdat = xdat, tydat = ystar, bws = bw,
        gradients = TRUE, se = FALSE
      )
      np:::.np_npsig_statistic(
        fit, index = tested.index, pivot = FALSE
      )
    }, numeric(1L))

    expect_equal(streamed, incumbent, tolerance = 2e-12)
  }
})
