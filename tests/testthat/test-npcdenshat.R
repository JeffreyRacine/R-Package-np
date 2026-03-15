test_that("npcdenshat matches npcdens and preserves matrix/apply parity across supported bwtype/regtype cells", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npcdenshat <- getFromNamespace("npcdenshat", "npRmpi")

  set.seed(20260310)
  n <- 18
  x <- data.frame(x = sort(runif(n)))
  y <- data.frame(y = rnorm(n))
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 8))
  ey <- data.frame(y = seq(min(y$y), max(y$y), length.out = 8))
  iota <- rep.int(1, n)

  cfgs <- expand.grid(
    bwtype = c("fixed", "generalized_nn", "adaptive_nn"),
    regtype = c("lc", "ll", "lp"),
    stringsAsFactors = FALSE
  )

  for (ii in seq_len(nrow(cfgs))) {
    cfg <- cfgs[ii, ]
    bw.args <- list(
      xdat = x,
      ydat = y,
      regtype = cfg$regtype,
      bwtype = cfg$bwtype,
      bws = if (identical(cfg$bwtype, "fixed")) c(0.24, 0.28) else c(7, 7),
      bandwidth.compute = FALSE
    )
    if (identical(cfg$regtype, "lp")) {
      bw.args$basis <- "glp"
      bw.args$degree <- 2L
    }

    bw <- do.call(npcdensbw, bw.args)
    fit <- npcdens(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey)
    H <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, output = "matrix")
    hy <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, y = iota, output = "apply")

    expect_s3_class(H, "npcdenshat")
    expect_equal(
      as.vector(H %*% iota),
      as.vector(fit$condens),
      tolerance = 1e-10,
      info = paste(cfg$bwtype, cfg$regtype)
    )
    expect_equal(
      as.vector(hy),
      as.vector(H %*% iota),
      tolerance = 1e-12,
      info = paste(cfg$bwtype, cfg$regtype)
    )
  }
})

test_that("npcdenshat apply mode matches matrix RHS multiplication", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npcdenshat <- getFromNamespace("npcdenshat", "npRmpi")

  set.seed(20260311)
  n <- 14
  x <- data.frame(x = sort(runif(n)))
  y <- data.frame(y = rnorm(n))
  ex <- data.frame(x = seq(0.15, 0.85, length.out = 6))
  ey <- data.frame(y = seq(min(y$y), max(y$y), length.out = 6))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bwtype = "generalized_nn",
    bws = c(7, 7),
    bandwidth.compute = FALSE
  )

  H <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, output = "matrix")
  rhs <- cbind(rep.int(1, n), seq_len(n) / n)
  hy <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, y = rhs, output = "apply")

  expect_equal(hy, H %*% rhs, tolerance = 1e-12)
})

test_that("npcdenshat matches npcdens on multivariate and mixed common-use evaluation cells", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npcdenshat <- getFromNamespace("npcdenshat", "npRmpi")

  set.seed(20260315)
  n <- 20
  x.cont <- data.frame(
    x1 = sort(runif(n)),
    x2 = sort(runif(n))
  )
  y.cont <- data.frame(y = rnorm(n))
  ex.cont <- data.frame(
    x1 = seq(0.1, 0.9, length.out = 7),
    x2 = seq(0.2, 0.8, length.out = 7)
  )
  ey.cont <- data.frame(y = seq(min(y.cont$y), max(y.cont$y), length.out = 7))

  x.mixed <- data.frame(
    x1 = runif(n),
    x2 = runif(n),
    uf = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    of = ordered(sample(c("low", "mid", "high"), n, replace = TRUE),
      levels = c("low", "mid", "high")
    )
  )
  y.mixed <- data.frame(y = rnorm(n))
  ex.mixed <- data.frame(
    x1 = seq(0.1, 0.9, length.out = 7),
    x2 = seq(0.15, 0.85, length.out = 7),
    uf = factor(rep(c("a", "b"), length.out = 7), levels = levels(x.mixed$uf)),
    of = ordered(rep(c("low", "high"), length.out = 7), levels = levels(x.mixed$of))
  )
  ey.mixed <- data.frame(y = seq(min(y.mixed$y), max(y.mixed$y), length.out = 7))

  cases <- list(
    list(
      label = "cont2 generalized ll",
      bw = npcdensbw(
        xdat = x.cont,
        ydat = y.cont,
        regtype = "ll",
        bwtype = "generalized_nn",
        bws = c(6, 7, 7),
        bandwidth.compute = FALSE
      ),
      txdat = x.cont,
      tydat = y.cont,
      exdat = ex.cont,
      eydat = ey.cont
    ),
    list(
      label = "mixed generalized lp2",
      bw = npcdensbw(
        xdat = x.mixed,
        ydat = y.mixed,
        regtype = "lp",
        basis = "glp",
        degree = c(2L, 2L),
        bwtype = "generalized_nn",
        bws = c(6, 7, 7, 0.45, 0.55),
        bandwidth.compute = FALSE
      ),
      txdat = x.mixed,
      tydat = y.mixed,
      exdat = ex.mixed,
      eydat = ey.mixed
    )
  )

  for (case in cases) {
    fit <- npcdens(
      bws = case$bw,
      txdat = case$txdat,
      tydat = case$tydat,
      exdat = case$exdat,
      eydat = case$eydat
    )
    H <- npcdenshat(
      bws = case$bw,
      txdat = case$txdat,
      tydat = case$tydat,
      exdat = case$exdat,
      eydat = case$eydat,
      output = "matrix"
    )
    rhs <- cbind(rep.int(1, nrow(case$txdat)), seq_len(nrow(case$txdat)) / nrow(case$txdat))
    hy <- npcdenshat(
      bws = case$bw,
      txdat = case$txdat,
      tydat = case$tydat,
      exdat = case$exdat,
      eydat = case$eydat,
      y = rhs,
      output = "apply"
    )

    expect_equal(
      as.vector(H %*% rhs[, 1L]),
      as.vector(fit$condens),
      tolerance = 1e-10,
      info = case$label
    )
    expect_equal(hy, H %*% rhs, tolerance = 1e-12, info = case$label)
  }
})

test_that("npcdenshat preserves bounded gaussian manual-bandwidth semantics", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  npcdenshat <- getFromNamespace("npcdenshat", "npRmpi")

  x <- data.frame(x = sort(c(5, 11, 21, 31, 46, 75, 98, 122, 145, 165, 195, 224)))
  y <- data.frame(y = sort(c(7, 15, 22, 35, 48, 59, 77, 95, 121, 141, 166, 188)))
  iota <- rep.int(1, nrow(x))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bws = c(100, 100),
    bandwidth.compute = FALSE,
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerbound = "range",
    cykerbound = "range"
  )

  fit <- npcdens(bws = bw, txdat = x, tydat = y)
  H <- npcdenshat(bws = bw, txdat = x, tydat = y, output = "matrix")
  hy <- npcdenshat(bws = bw, txdat = x, tydat = y, y = iota, output = "apply")

  expect_equal(as.vector(H %*% iota), as.vector(fit$condens), tolerance = 1e-8)
  expect_equal(as.vector(hy), as.vector(fit$condens), tolerance = 1e-8)
})
