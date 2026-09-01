test_that("streamed IID tiles reproduce incumbent categorical and HC0 statistics", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

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
      npRmpi:::cast(xq, xdat[[tested.index]], same.levels = TRUE)
    } else xq
    null.mean <- npreg(
      txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
      gradients = FALSE, se = FALSE
    )$mean
    residual.pool <- as.numeric(
      scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
    )
    residual.pool <- residual.pool - mean(residual.pool)

    streamed <- npRmpi:::.np_npsig_streamed_iid_tile(
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
      npRmpi:::.np_npsig_statistic(
        fit,
        index = tested.index,
        pivot = !(is.factor(xdat[[tested.index]]) ||
                    is.ordered(xdat[[tested.index]]))
      )
    }, numeric(1L))

    expect_equal(streamed, incumbent, tolerance = 2e-12)
  }
})

test_that("streamed IID tiles preserve general LP basis and degree semantics", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260902)
  n <- 42L
  xdat <- data.frame(
    z = factor(rep(c("a", "b", "c"), length.out = n)),
    x1 = seq(-1.1, 1.2, length.out = n),
    x2 = sin(seq(-1.0, 1.3, length.out = n))
  )
  ydat <- 0.4 + c(a = 0, b = 0.3, c = -0.2)[xdat$z] +
    0.7 * xdat$x1 - 0.25 * xdat$x2^2 + rnorm(n, sd = 0.12)
  specifications <- list(
    list(basis = "glp", degree = c(1L, 2L), shifted = FALSE),
    list(basis = "glp", degree = c(1L, 2L), shifted = TRUE),
    list(basis = "additive", degree = c(2L, 2L), shifted = TRUE),
    list(basis = "tensor", degree = c(1L, 1L), shifted = TRUE)
  )
  set.seed(771)
  donor <- vapply(
    seq_len(2L),
    function(unused) sample.int(n, replace = TRUE),
    integer(n)
  )

  for (specification in specifications) {
    bw <- npregbw(
      xdat = xdat, ydat = ydat,
      bws = c(0.2, 0.55, 0.6), bandwidth.compute = FALSE,
      regtype = "lp", bwtype = "fixed",
      basis = specification$basis,
      degree = specification$degree,
      bernstein.basis = specification$shifted
    )
    expect_true(npRmpi:::.np_npsig_streamed_iid_eligible(
      bw, xdat, seq_len(ncol(xdat)), FALSE, "I", "iid", NULL, list()
    ))
    unrestricted <- npreg(
      txdat = xdat, tydat = ydat, bws = bw,
      gradients = FALSE, se = FALSE
    )
    unrestricted.scaled <- scale(residuals(unrestricted))
    unrestricted.scale <- attr(unrestricted.scaled, "scaled:scale")
    unrestricted.center <- attr(unrestricted.scaled, "scaled:center")

    for (tested.index in seq_len(ncol(xdat))) {
      null.frame <- xdat
      xq <- uocquantile(xdat[[tested.index]], 0.5)
      null.frame[[tested.index]] <- if (is.factor(xdat[[tested.index]])) {
        npRmpi:::cast(xq, xdat[[tested.index]], same.levels = TRUE)
      } else xq
      null.mean <- npreg(
        txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
        gradients = FALSE, se = FALSE
      )$mean
      residual.pool <- as.numeric(
        scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
      )
      residual.pool <- residual.pool - mean(residual.pool)
      categorical <- is.factor(xdat[[tested.index]]) ||
        is.ordered(xdat[[tested.index]])
      streamed <- npRmpi:::.np_npsig_streamed_iid_tile(
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
          gradients = TRUE, se = !categorical
        )
        npRmpi:::.np_npsig_statistic(
          fit, index = tested.index, pivot = !categorical
        )
      }, numeric(1L))
      expect_equal(streamed, incumbent, tolerance = 1e-9)
    }
  }
})

test_that("streamed IID tiles preserve generalized and adaptive NN geometry", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260903)
  n <- 40L
  xdat <- data.frame(
    z = factor(rep(c("a", "b"), length.out = n)),
    x1 = seq(-1.15, 1.25, length.out = n),
    x2 = cos(seq(-0.9, 1.4, length.out = n))
  )
  ydat <- 0.5 + 0.35 * (xdat$z == "b") +
    0.8 * xdat$x1 - 0.3 * xdat$x2 + rnorm(n, sd = 0.1)
  set.seed(773)
  donor <- vapply(
    seq_len(2L),
    function(unused) sample.int(n, replace = TRUE),
    integer(n)
  )

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    bw <- npregbw(
      xdat = xdat, ydat = ydat,
      bws = c(0.2, 9, 10), bandwidth.compute = FALSE,
      regtype = "ll", bwtype = bwtype
    )
    expect_true(npRmpi:::.np_npsig_streamed_iid_eligible(
      bw, xdat, seq_len(ncol(xdat)), FALSE, "I", "iid", NULL, list()
    ))
    unrestricted <- npreg(
      txdat = xdat, tydat = ydat, bws = bw,
      gradients = FALSE, se = FALSE
    )
    unrestricted.scaled <- scale(residuals(unrestricted))
    unrestricted.scale <- attr(unrestricted.scaled, "scaled:scale")
    unrestricted.center <- attr(unrestricted.scaled, "scaled:center")

    for (tested.index in seq_len(ncol(xdat))) {
      null.frame <- xdat
      xq <- uocquantile(xdat[[tested.index]], 0.5)
      null.frame[[tested.index]] <- if (is.factor(xdat[[tested.index]])) {
        npRmpi:::cast(xq, xdat[[tested.index]], same.levels = TRUE)
      } else xq
      null.mean <- npreg(
        txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
        gradients = FALSE, se = FALSE
      )$mean
      residual.pool <- as.numeric(
        scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
      )
      residual.pool <- residual.pool - mean(residual.pool)
      categorical <- is.factor(xdat[[tested.index]]) ||
        is.ordered(xdat[[tested.index]])
      streamed <- npRmpi:::.np_npsig_streamed_iid_tile(
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
          gradients = TRUE, se = !categorical
        )
        npRmpi:::.np_npsig_statistic(
          fit, index = tested.index, pivot = !categorical
        )
      }, numeric(1L))
      expect_equal(
        streamed,
        incumbent,
        tolerance = if (identical(bwtype, "adaptive_nn") && !categorical)
          5e-7 else 2e-10
      )
    }
  }
})

test_that("streamed IID tiles reproduce mixed local-constant statistics", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

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

  expect_true(npRmpi:::.np_npsig_streamed_iid_eligible(
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
      npRmpi:::cast(xq, xdat[[tested.index]], same.levels = TRUE)
    } else xq
    null.mean <- npreg(
      txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
      gradients = FALSE, se = FALSE
    )$mean
    residual.pool <- as.numeric(
      scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
    )
    residual.pool <- residual.pool - mean(residual.pool)

    streamed <- npRmpi:::.np_npsig_streamed_iid_tile(
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
      npRmpi:::.np_npsig_statistic(
        fit, index = tested.index, pivot = !categorical
      )
    }, numeric(1L))

    expect_equal(streamed, incumbent, tolerance = 2e-12)
  }
})

test_that("streamed response tiles reproduce wild-bootstrap statistics", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260904)
  n <- 42L
  xdat <- data.frame(
    z = factor(rep(c("a", "b", "c"), length.out = n)),
    x1 = seq(-1.1, 1.2, length.out = n),
    x2 = sin(seq(-0.8, 1.3, length.out = n))
  )
  ydat <- 0.4 * (xdat$z == "b") + 0.8 * xdat$x1 -
    0.2 * xdat$x2 + rnorm(n, sd = 0.08)
  bw <- npregbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.2, 0.6, 0.55), bandwidth.compute = FALSE,
    regtype = "ll", bwtype = "fixed"
  )
  unrestricted <- npreg(
    txdat = xdat, tydat = ydat, bws = bw,
    gradients = FALSE, se = FALSE
  )
  unrestricted.scaled <- scale(residuals(unrestricted))
  unrestricted.scale <- attr(unrestricted.scaled, "scaled:scale")
  unrestricted.center <- attr(unrestricted.scaled, "scaled:center")

  for (tested.index in seq_len(ncol(xdat))) {
    null.frame <- xdat
    xq <- uocquantile(xdat[[tested.index]], 0.5)
    null.frame[[tested.index]] <- if (is.factor(xdat[[tested.index]])) {
      npRmpi:::cast(xq, xdat[[tested.index]], same.levels = TRUE)
    } else xq
    null.mean <- npreg(
      txdat = xdat, tydat = ydat, exdat = null.frame, bws = bw,
      gradients = FALSE, se = FALSE
    )$mean
    residual.pool <- as.numeric(
      scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
    )
    residual.pool <- residual.pool - mean(residual.pool)
    set.seed(700L + tested.index)
    response <- vapply(seq_len(3L), function(unused) {
      multiplier <- ifelse(runif(n) <= 0.72360679774998,
                           -0.6180339887499, 1.6180339887499)
      null.mean + residual.pool * multiplier
    }, numeric(n))
    streamed <- npRmpi:::.np_npsig_streamed_iid_tile(
      bws = bw,
      xdat = xdat,
      tested.index = tested.index,
      response.matrix = response,
      null.mean = null.mean,
      residual.pool = residual.pool
    )
    categorical <- is.factor(xdat[[tested.index]]) ||
      is.ordered(xdat[[tested.index]])
    incumbent <- vapply(seq_len(ncol(response)), function(replication) {
      fit <- npreg(
        txdat = xdat, tydat = response[, replication], bws = bw,
        gradients = TRUE, se = !categorical
      )
      npRmpi:::.np_npsig_statistic(
        fit, index = tested.index, pivot = !categorical
      )
    }, numeric(1L))
    expect_equal(streamed, incumbent, tolerance = 2e-10)
    if (!categorical) {
      streamed.raw <- npRmpi:::.np_npsig_streamed_iid_tile(
        bws = bw,
        xdat = xdat,
        tested.index = tested.index,
        response.matrix = response,
        null.mean = null.mean,
        residual.pool = residual.pool,
        pivotal = FALSE
      )
      incumbent.raw <- vapply(
        seq_len(ncol(response)),
        function(replication) {
          fit <- npreg(
            txdat = xdat, tydat = response[, replication], bws = bw,
            gradients = TRUE, se = FALSE
          )
          npRmpi:::.np_npsig_statistic(
            fit, index = tested.index, pivot = FALSE
          )
        },
        numeric(1L)
      )
      expect_equal(streamed.raw, incumbent.raw, tolerance = 2e-10)
    }
  }
})

test_that("streamed IID capability is deterministic and semantics-exact", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

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
  eligible <- function(...) npRmpi:::.np_npsig_streamed_iid_eligible(
    bws = bw, xdat = xdat, index = 1:2, joint = FALSE,
    boot.type = "I", boot.method = "iid", pivot = NULL,
    extra.args = list(), ...
  )

  expect_true(eligible())
  expect_true(npRmpi:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, FALSE, "I", "wild", NULL, list()
  ))
  expect_true(npRmpi:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, FALSE, "I", "wild-rademacher", NULL, list()
  ))
  expect_false(npRmpi:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, FALSE, "I", "pairwise", NULL, list()
  ))
  expect_false(npRmpi:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 2L, FALSE, "I", "iid", FALSE, list()
  ))
  expect_true(npRmpi:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, TRUE, "I", "iid", NULL, list()
  ))
  expect_true(npRmpi:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, TRUE, "I", "wild", FALSE, list()
  ))
  expect_false(npRmpi:::.np_npsig_streamed_iid_eligible(
    bw, xdat, 1:2, TRUE, "I", "iid", TRUE, list()
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

test_that("direct individual pairwise statistics reproduce public fits", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260905)
  n <- 30L
  xdat <- data.frame(
    z = factor(rep(c("a", "b"), length.out = n)),
    x = seq(-1, 1, length.out = n)
  )
  ydat <- 0.4 * (xdat$z == "b") + 0.7 * xdat$x + rnorm(n, sd = 0.1)
  bw <- npregbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.2, 0.55), bandwidth.compute = FALSE,
    regtype = "ll", bwtype = "fixed"
  )
  seed <- 867L
  result <- npsigtest(
    bw, xdat = xdat, ydat = ydat,
    B = 9L, boot.method = "pairwise", random.seed = seed
  )

  set.seed(seed)
  oracle <- matrix(0, nrow = 9L, ncol = ncol(xdat))
  for (tested.index in seq_len(ncol(xdat))) {
    categorical <- is.factor(xdat[[tested.index]]) ||
      is.ordered(xdat[[tested.index]])
    for (replication in seq_len(nrow(oracle))) {
      donor <- sample.int(n, replace = TRUE)
      xstar <- xdat
      xstar[, -tested.index] <- xdat[donor, -tested.index, drop = FALSE]
      fit <- npreg(
        txdat = xstar, tydat = ydat[donor], bws = bw,
        gradients = TRUE, se = !categorical
      )
      oracle[replication, tested.index] <- npRmpi:::.np_npsig_statistic(
        fit, index = tested.index, pivot = !categorical
      )
    }
  }
  expect_equal(result$In.bootstrap, oracle, tolerance = 2e-10)
})

test_that("streamed IID categorical-only tiles reproduce scalar-LP fits", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

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
    null.frame[[tested.index]] <- npRmpi:::cast(
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

    streamed <- npRmpi:::.np_npsig_streamed_iid_tile(
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
      npRmpi:::.np_npsig_statistic(
        fit, index = tested.index, pivot = FALSE
      )
    }, numeric(1L))

    expect_equal(streamed, incumbent, tolerance = 2e-12)
  }
})
