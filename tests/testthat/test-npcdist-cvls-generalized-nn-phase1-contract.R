library(npRmpi)

phase1_npcdist_cvls_gnn_fixture <- function() {
  set.seed(20260309)
  n <- 24L
  x <- data.frame(
    x1 = runif(n),
    x2 = runif(n)
  )
  y <- data.frame(
    y1 = x$x1^2 - 0.35 * x$x2 + 0.2 * sin(2 * pi * x$x1) + rnorm(n, sd = 0.08)
  )
  list(x = x, y = y)
}

# Independent literal deleted-sample external-grid reference for LC/raw degree1.
phase1_npcdist_cvls_gnn_oracle <- function(dat, bw, degree) {
  x <- as.matrix(dat$x)
  y <- dat$y$y1
  n <- nrow(x)
  grid <- as.numeric(stats::quantile(y, probs = seq(0, 1, length.out = 100L)))
  stopifnot(degree %in% 0:1, length(bw$xbw) == ncol(x),
            all(bw$xbw == floor(bw$xbw)), all(bw$xbw >= 1),
            all(bw$xbw <= n - 1), bw$ybw == floor(bw$ybw),
            bw$ybw >= 1, bw$ybw <= n - 1)
  basis <- if (degree == 0L) matrix(1, n, 1L) else cbind(1, x)
  loss <- 0
  for (i in seq_len(n)) {
    keep <- setdiff(seq_len(n), i)
    weight <- rep(1, n - 1L)
    for (d in seq_len(ncol(x))) {
      radius <- sort(abs(x[keep, d] - x[i, d]))[[bw$xbw[[d]]]]
      stopifnot(is.finite(radius), radius > 0)
      weight <- weight * stats::dnorm((x[keep, d] - x[i, d]) / radius)
    }
    decomposition <- qr(sqrt(weight) * basis[keep, , drop = FALSE])
    stopifnot(decomposition$rank == ncol(basis))
    for (z in grid) {
      radius <- sort(abs(y[keep] - z))[[bw$ybw]]
      stopifnot(is.finite(radius), radius > 0)
      response <- stats::pnorm((z - y[keep]) / radius)
      coefficients <- qr.coef(decomposition, sqrt(weight) * response)
      prediction <- sum(basis[i, ] * coefficients)
      loss <- loss + (as.numeric(y[i] <= z) - prediction)^2
    }
  }
  loss / (n * length(grid))
}

phase1_npcdist_cvls_gnn_cases <- local({
  cache <- NULL

  function() {
    if (!is.null(cache))
      return(cache)
    if (!spawn_mpi_slaves())
      skip("Could not spawn MPI slaves")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

    dat <- phase1_npcdist_cvls_gnn_fixture()
    degree1 <- rep.int(1L, ncol(dat$x))
    degree2 <- rep.int(2L, ncol(dat$x))

    cache <<- list(
      dat = dat,
      degree1 = degree1,
      degree2 = degree2,
      bw.lc = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "lc",
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      ),
      bw.lp0 = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "lp",
        basis = "glp",
        degree = rep.int(0L, ncol(dat$x)),
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      ),
      bw.ll = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "ll",
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      ),
      bw.lp = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "lp",
        basis = "glp",
        degree = degree1,
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      ),
      bw.d2 = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "lp",
        basis = "glp",
        degree = degree2,
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      )
    )

    cache
  }
})

test_that("phase1 npcdistbw cv.ls generalized-nn lc matches literal external-grid CV", {
  cases <- phase1_npcdist_cvls_gnn_cases()
  bw.lc <- cases$bw.lc
  bw.lp0 <- cases$bw.lp0

  expect_true(is.finite(bw.lc$fval))
  expect_equal(bw.lc$fval, phase1_npcdist_cvls_gnn_oracle(cases$dat, bw.lc, 0L),
               tolerance = 1e-10)
  expect_identical(bw.lc[["fval"]], bw.lp0[["fval"]])
  expect_identical(bw.lc[["xbw"]], bw.lp0[["xbw"]])
  expect_identical(bw.lc[["ybw"]], bw.lp0[["ybw"]])
})

test_that("phase1 npcdistbw cv.ls generalized-nn keeps ll on canonical lp degree-1 glp", {
  cases <- phase1_npcdist_cvls_gnn_cases()
  degree <- cases$degree1
  bw.ll <- cases$bw.ll
  bw.lp <- cases$bw.lp

  expect_identical(bw.ll$regtype.engine, "lp")
  expect_identical(bw.ll$basis.engine, "glp")
  expect_identical(as.integer(bw.ll$degree.engine), degree)
  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, phase1_npcdist_cvls_gnn_oracle(cases$dat, bw.ll, 1L),
               tolerance = 1e-10)
  expect_equal(bw.lp$fval, phase1_npcdist_cvls_gnn_oracle(cases$dat, bw.lp, 1L),
               tolerance = 1e-10)
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-10)
})

test_that("phase1 npcdistbw cv.ls generalized-nn lp degree-2 succeeds on a higher-order fixture", {
  cases <- phase1_npcdist_cvls_gnn_cases()
  bw.d1 <- cases$bw.lp
  bw.d2 <- cases$bw.d2

  expect_identical(as.integer(bw.d2$degree.engine), cases$degree2)
  expect_true(is.finite(bw.d2$fval))
  expect_gte(bw.d2$fval, 0)
  expect_true(all(is.finite(c(bw.d2$xbw, bw.d2$ybw))))
  expect_gt(abs(bw.d2$fval - bw.d1$fval), 1e-6)
})

test_that("phase1 npcdistbw cv.ls generalized-nn avoids search-boundary collapse", {
  if (!spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260309)
  n <- 30L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = rbeta(n, shape1 = 1 + 2 * x$x1, shape2 = 2))

  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1L,
    itmax = 1L
  )
  fit <- npcdist(
    bws = bw,
    exdat = data.frame(x1 = rep(0.5, 4L)),
    eydat = data.frame(y1 = c(0, 0.02, 0.98, 1))
  )
  pred <- fitted(fit)

  expect_true(is.finite(bw$fval))
  expect_true(all(is.finite(pred)))
  expect_true(all(pred >= -1e-6))
  expect_true(all(pred <= 1 + 1e-6))
  expect_gt(pred[3] - pred[1], 0.5)
  expect_true(all(bw$xbw > 1))
  expect_true(all(bw$ybw > 1))
  expect_true(all(bw$xbw < n))
  expect_true(all(bw$ybw < n))
})
