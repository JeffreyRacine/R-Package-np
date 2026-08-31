test_that("npindexbw raw certification excludes mapped and terminal sentinels", {
  valid <- getFromNamespace(".npindexbw_raw_objective_valid", "np")
  result <- getFromNamespace(".npindexbw_objective_result", "np")

  expect_true(valid(0))
  expect_true(valid(-1))
  expect_false(valid(.Machine$double.xmax))
  expect_false(valid(Inf))
  expect_false(valid(NA_real_))
  expect_false(valid(numeric()))
  expect_null(result(.Machine$double.xmax, 0L)$raw.valid)
  expect_true(result(1, 0L, certify = TRUE)$raw.valid)
  expect_false(result(.Machine$double.xmax, 0L, certify = TRUE)$raw.valid)
})

make_npindex_invalid_legacy_fixture <- function(method) {
  x <- data.frame(
    x1 = seq(-1, 1, length.out = 24L),
    x2 = seq(1, -1, length.out = 24L)^2
  )
  y <- if (identical(method, "kleinspady")) {
    rep(c(0, 1), 12L)
  } else {
    sin(3 * x$x1) + seq_len(24L) / 1000
  }
  bws <- npindexbw(
    xdat = x,
    ydat = y,
    method = method,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    bandwidth.compute = FALSE
  )
  bws$beta <- c(1, 0.5)
  bws$bw <- 1e-9
  bws$bandwidth[[1L]] <- 1e-9
  list(xdat = x, ydat = y, bws = bws)
}

test_that("npindexbw legacy owners reject an invalid selected candidate", {
  for (method in c("ichimura", "kleinspady")) {
    fixture <- make_npindex_invalid_legacy_fixture(method)
    expect_error(
      npindexbw(
        xdat = fixture$xdat,
        ydat = fixture$ydat,
        bws = fixture$bws,
        only.optimize.beta = TRUE,
        nmulti = 1L,
        optim.maxit = 5L,
        scale.factor.search.lower = 0
      ),
      "npindexbw search did not return a raw-valid selected candidate",
      fixed = TRUE,
      info = method
    )
  }
})

make_npindex_invalid_nomad_fixture <- function(engine) {
  set.seed(941L)
  x <- data.frame(
    x1 = seq(-1, 1, length.out = 24L),
    x2 = runif(24L, -1, 1)
  )
  list(
    xdat = x,
    ydat = sin(3 * x$x1) + 0.2 * x$x2 + seq_len(24L) / 1000,
    bws = c(1, 0.5, 1e-9),
    method = "ichimura",
    regtype = "lp",
    nomad = TRUE,
    search.engine = engine,
    degree.min = 1L,
    degree.max = 2L,
    degree.start = 1L,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    nmulti = 1L,
    only.optimize.beta = identical(engine, "nomad+powell"),
    optim.maxit = 5L,
    powell.remin = FALSE,
    scale.factor.search.lower = 0,
    scale.factor.init.lower = 1e-12,
    nomad.opts = list(MAX_BB_EVAL = 1L)
  )
}

test_that("npindexbw NOMAD boundaries reject an invalid selected candidate", {
  for (engine in c("nomad", "nomad+powell")) {
    expect_error(
      do.call(npindexbw, make_npindex_invalid_nomad_fixture(engine)),
      "npindexbw search did not return a raw-valid selected candidate",
      fixed = TRUE,
      info = engine
    )
  }
})
