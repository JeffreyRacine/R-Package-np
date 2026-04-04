npscoef_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

npscoef_fit_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

npscoef_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

make_npscoef_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 24L
  dat <- data.frame(
    x = runif(n),
    z = sort(runif(n))
  )
  dat$y <- (1 + dat$z^2) * dat$x + rnorm(n, sd = 0.08)

  list(
    dat = dat,
    tx = dat["x"],
    tz = dat["z"],
    y = dat$y,
    bw = npscoefbw(
      xdat = dat["x"],
      zdat = dat["z"],
      ydat = dat$y,
      bws = 0.20,
      bandwidth.compute = FALSE
    )
  )
}

test_that("npscoef fit progress uses the bandwidth single-line surface", {
  begin <- getFromNamespace(".np_scoef_fit_progress_begin", "np")

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  state <- with_np_progress_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = function() 0
    ),
    begin()
  )

  expect_identical(state$surface, "bandwidth")
  expect_identical(state$renderer, "single_line")
  expect_null(state$total)
})

test_that("npscoef direct bws fit emits fit progress", {
  fixture <- make_npscoef_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npscoef(
      bws = fixture$bw,
      txdat = fixture$tx,
      tzdat = fixture$tz,
      tydat = fixture$y,
      errors = FALSE,
      iterate = FALSE
    ),
    force_renderer = "single_line",
    now = npscoef_fit_progress_time_counter()
  )

  lines <- npscoef_fit_progress_lines(actual)

  expect_s3_class(actual$value, "smoothcoefficient")
  expect_true(any(grepl(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\.$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\. iteration 1, elapsed [0-9]+\\.[0-9]s: solving coefficient rows$",
    lines
  )))
})

test_that("npscoef direct bws fit stays silent below start grace without handoff", {
  fixture <- make_npscoef_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.unknown.sec = 0.75,
    np.progress.interval.unknown.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npscoef(
      bws = fixture$bw,
      txdat = fixture$tx,
      tzdat = fixture$tz,
      tydat = fixture$y,
      errors = FALSE,
      iterate = FALSE
    ),
    force_renderer = "single_line",
    now = npscoef_fit_progress_time_values(c(0, 0.2, 0.4, 0.6))
  )

  expect_length(actual$trace, 0L)
})

test_that("npscoef formula bw to fit route hands off immediately into fit progress", {
  fixture <- make_npscoef_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npscoef(
      y ~ x | z,
      data = fixture$dat,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = npscoef_fit_progress_time_counter()
  )

  lines <- npscoef_fit_progress_lines(actual)
  fit.start.pos <- grep(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\. elapsed 0\\.0s: building moments$",
    lines
  )
  fit.solve.pos <- grep(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\. iteration 1, elapsed [0-9]+\\.[0-9]s: solving coefficient rows$",
    lines
  )
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)

  expect_s3_class(actual$value, "smoothcoefficient")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.solve.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.solve.pos[[1L]])
})

test_that("npscoef nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")

  fixture <- make_npscoef_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npscoef(
      y ~ x | z,
      data = fixture$dat,
      regtype = "lp",
      degree.select = "coordinate",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L,
      errors = FALSE,
      iterate = FALSE
    ),
    force_renderer = "single_line",
    now = npscoef_fit_progress_time_counter()
  )

  lines <- npscoef_fit_progress_lines(actual)
  fit.start.pos <- grep(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\. elapsed 0\\.0s: building moments$",
    lines
  )
  fit.solve.pos <- grep(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\. iteration 1, elapsed [0-9]+\\.[0-9]s: solving coefficient rows$",
    lines
  )
  powell.pos <- grep("^\\[np\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[np\\] Selecting degree and bandwidth", lines)

  expect_s3_class(actual$value, "smoothcoefficient")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.solve.pos) >= 1L)
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.solve.pos[[1L]])
})

test_that("predict.smoothcoefficient re-entry emits fit progress", {
  fixture <- make_npscoef_fit_progress_fixture()

  fit <- npscoef(
    bws = fixture$bw,
    txdat = fixture$tx,
    tzdat = fixture$tz,
    tydat = fixture$y,
    errors = FALSE,
    iterate = FALSE
  )

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    predict(fit, newdata = fixture$dat[c(2L, 7L), c("x", "z")]),
    force_renderer = "single_line",
    now = npscoef_fit_progress_time_counter()
  )

  lines <- npscoef_fit_progress_lines(actual)

  expect_true(any(grepl(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\. iteration 1, elapsed [0-9]+\\.[0-9]s: solving coefficient rows$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[np\\] Fitting smooth coefficient model\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s: estimating standard errors$",
    lines
  )))
})
