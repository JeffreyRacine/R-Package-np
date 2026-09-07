test_that("npRmpi npindexbw raw certification excludes mapped and terminal sentinels", {
  valid <- getFromNamespace(".npindexbw_raw_objective_valid", "npRmpi")
  result <- getFromNamespace(".npindexbw_objective_result", "npRmpi")

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

make_nprmpi_npindex_invalid_legacy_fixture <- function(method) {
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

test_that("npRmpi npindexbw legacy owners reject a raw-invalid held start", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  for (method in c("ichimura", "kleinspady")) {
    fixture <- make_nprmpi_npindex_invalid_legacy_fixture(method)
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
      "raw-invalid held bandwidth; restoration is disabled",
      fixed = TRUE,
      info = method
    )
  }
})

make_nprmpi_npindex_invalid_nomad_fixture <- function(engine) {
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

test_that("npRmpi npindexbw NOMAD boundaries reject invalid starts or endpoints", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  for (engine in c("nomad", "nomad+powell")) {
    expect_error(
      do.call(npindexbw, make_nprmpi_npindex_invalid_nomad_fixture(engine)),
      if (engine == "nomad+powell")
        "raw-invalid held bandwidth; restoration is disabled" else
        "npindexbw search did not return a raw-valid selected candidate",
      fixed = TRUE,
      info = engine
    )
  }
})

test_that("npRmpi npindexbw final certificate preserves the point and rejects raw invalidity", {
  certify <- getFromNamespace(".npindexbw_certify_selected_candidate", "npRmpi")
  result <- getFromNamespace(".npindexbw_objective_result", "npRmpi")
  point <- c(0.5, 0.75)
  seen <- NULL
  raw <- 0.25
  testthat::local_mocked_bindings(
    # MPI has a dedicated raw evaluator, rather than the serial certify flag.
    .npindexbw_eval_objective_raw = function(param, xmat, ydat, bws, spec,
                                            leaf.descriptor) {
      seen <<- param
      result(raw, 1L, TRUE)
    },
    .package = "npRmpi"
  )
  args <- list(param = point, xmat = matrix(1, 2L, 2L), ydat = c(0, 1),
               bws = list(), spec = list())
  expect_identical(do.call(certify, args), result(raw, 1L, TRUE))
  expect_identical(seen, point)
  raw <- .Machine$double.xmax
  seen <- NULL
  expect_error(do.call(certify, args),
               "npindexbw search did not return a raw-valid selected candidate",
               fixed = TRUE)
  expect_identical(seen, point)
})
