auto_beta_arguments <- function(dat, degree.start, bws = NULL) {
  arguments <- list(
    formula = y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    degree.start = degree.start,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    ckertype = "beta",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 10L)
  )
  if (!is.null(bws))
    arguments[["bws"]] <- bws
  arguments
}

auto_beta_replay <- function(dat, bandwidth) {
  getFromNamespace(".npregbw_eval_only", "np")(
    xdat = dat["x"],
    ydat = dat[["y"]],
    bws = bandwidth,
    invalid.penalty = "baseline"
  )[["objective"]][[1L]]
}

test_that("automatic beta regression starts with search-owned bandwidths", {
  skip_if_not_installed("crs")

  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.04, 0.96, length.out = 32L)
  dat <- data.frame(y = sin(2 * pi * x) + 0.04 * cos(7 * pi * x), x = x)
  bandwidth <- do.call(npregbw, auto_beta_arguments(dat, degree.start = 1L))
  replay <- auto_beta_replay(dat, bandwidth)

  expect_true(all(is.finite(bandwidth[["bw"]])))
  expect_true(all(bandwidth[["bw"]] > 0))
  expect_identical(as.double(bandwidth[["fval"]][[1L]]), as.double(replay))
})

test_that("automatic beta degree zero uses the canonical scalar objective", {
  skip_if_not_installed("crs")

  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.04, 0.96, length.out = 32L)
  dat <- data.frame(y = sin(2 * pi * x) + 0.04 * cos(7 * pi * x), x = x)
  automatic <- do.call(
    npregbw,
    auto_beta_arguments(dat, degree.start = 0L, bws = 0.17)
  )
  restart <- automatic[["degree.search"]][["restart.results"]][[1L]]

  fixed <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 0L,
    bandwidth.compute = FALSE,
    bws = 0.17,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    ckertype = "beta",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )
  fixed.objective <- auto_beta_replay(dat, fixed)
  replay <- auto_beta_replay(dat, automatic)

  expect_identical(as.integer(restart[["first_degree"]]), 0L)
  expect_identical(
    as.double(restart[["first_objective"]][[1L]]),
    as.double(fixed.objective)
  )
  expect_identical(as.double(automatic[["fval"]][[1L]]), as.double(replay))
})
