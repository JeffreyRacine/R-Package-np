library(np)

test_that("npcdenshat matches npcdens and preserves matrix/apply parity across fixed-bandwidth x-side cells", {
  npcdenshat <- getFromNamespace("npcdenshat", "np")

  set.seed(20260310)
  n <- 18
  x <- data.frame(x = sort(runif(n)))
  y <- data.frame(y = rnorm(n))
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 8))
  ey <- data.frame(y = seq(min(y$y), max(y$y), length.out = 8))
  iota <- rep.int(1, n)

  cfgs <- list(
    list(regtype = "lc"),
    list(regtype = "ll"),
    list(regtype = "lp", basis = "glp", degree = 2L)
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = x,
      ydat = y,
      regtype = cfg$regtype,
      bwtype = "fixed",
      bws = c(0.24, 0.28),
      bandwidth.compute = FALSE
    )
    if (!is.null(cfg$basis))
      bw.args$basis <- cfg$basis
    if (!is.null(cfg$degree))
      bw.args$degree <- cfg$degree

    bw <- do.call(npcdensbw, bw.args)
    fit <- npcdens(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey)
    H <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, output = "matrix")
    hy <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, y = iota, output = "apply")

    expect_s3_class(H, "npcdenshat")
    expect_equal(as.vector(H %*% iota), as.vector(fit$condens), tolerance = 1e-10,
                 info = cfg$regtype)
    expect_equal(as.vector(hy), as.vector(H %*% iota), tolerance = 1e-12,
                 info = cfg$regtype)
  }
})

test_that("npcdenshat apply mode matches matrix RHS multiplication", {
  npcdenshat <- getFromNamespace("npcdenshat", "np")

  set.seed(20260311)
  n <- 14
  x <- data.frame(x = sort(runif(n)))
  y <- data.frame(y = rnorm(n))
  ex <- data.frame(x = seq(0.15, 0.85, length.out = 6))
  ey <- data.frame(y = seq(min(y$y), max(y$y), length.out = 6))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "fixed",
    bws = c(0.22, 0.3),
    bandwidth.compute = FALSE
  )

  H <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, output = "matrix")
  rhs <- cbind(rep.int(1, n), seq_len(n) / n)
  hy <- npcdenshat(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey, y = rhs, output = "apply")

  expect_equal(hy, H %*% rhs, tolerance = 1e-12)
})

test_that("npcdenshat preserves bounded gaussian manual-bandwidth semantics", {
  npcdenshat <- getFromNamespace("npcdenshat", "np")

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

test_that("npcdenshat fails fast on non-fixed bandwidth types", {
  npcdenshat <- getFromNamespace("npcdenshat", "np")

  set.seed(20260312)
  x <- data.frame(x = sort(runif(12)))
  y <- data.frame(y = rnorm(12))

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      regtype = "ll",
      bwtype = bwtype,
      bws = c(7, 7),
      bandwidth.compute = FALSE
    )

    expect_error(
      npcdenshat(bws = bw, txdat = x, tydat = y, output = "matrix"),
      "supports bwtype = 'fixed' only"
    )
  }
})
