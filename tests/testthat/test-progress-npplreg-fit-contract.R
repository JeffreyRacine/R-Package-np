npplreg_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

npplreg_fit_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

npplreg_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

make_npplreg_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 32L
  dat <- data.frame(
    x = seq(-0.9, 0.9, length.out = n),
    z = seq(0.2, 1.8, length.out = n)
  )
  dat$y <- 1.5 * dat$x + sin(2 * pi * dat$z) + 0.1 * cos(pi * dat$x)

  list(
    dat = dat,
    tx = dat["x"],
    tz = dat["z"],
    y = dat$y,
    bw = npplregbw(
      xdat = dat["x"],
      zdat = dat["z"],
      ydat = dat$y,
      bws = matrix(c(0.35, 0.35), nrow = 2L),
      bandwidth.compute = FALSE
    )
  )
}

test_that("npplreg fit progress uses the bandwidth single-line surface", {
  begin <- getFromNamespace(".np_plreg_fit_progress_begin", "np")

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  state <- with_np_progress_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = function() 0
    ),
    begin(xnames = "x")
  )

  expect_identical(state$surface, "bandwidth")
  expect_identical(state$renderer, "single_line")
  expect_equal(state$total, 2L)
})

test_that("npplreg direct bws fit emits known-total fit progress", {
  fixture <- make_npplreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npplreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tzdat = fixture$tz,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npplreg_fit_progress_time_counter()
  )

  lines <- npplreg_fit_progress_lines(actual)

  expect_s3_class(actual$value, "plregression")
  expect_true(any(grepl(
    "^\\[np\\] Fitting partially linear regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): y~z$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[np\\] Fitting partially linear regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\): x~z$",
    lines
  )))
})

test_that("npplreg direct bws fit stays silent below start grace without handoff", {
  fixture <- make_npplreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0.75,
    np.progress.interval.known.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npplreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tzdat = fixture$tz,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npplreg_fit_progress_time_values(c(0, 0.2, 0.4, 0.6))
  )

  expect_length(actual$trace, 0L)
})

test_that("npplreg formula bw to fit route hands off immediately into fit progress", {
  fixture <- make_npplreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npplreg(
      y ~ x | z,
      data = fixture$dat,
      nmulti = 1
    ),
    force_renderer = "single_line",
    now = npplreg_fit_progress_time_counter()
  )

  lines <- npplreg_fit_progress_lines(actual)
  fit.zero.pos <- grep(
    "^\\[np\\] Fitting partially linear regression 0/2 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting y~z$",
    lines
  )
  fit.one.pos <- grep(
    "^\\[np\\] Fitting partially linear regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): y~z$",
    lines
  )
  fit.two.pos <- grep(
    "^\\[np\\] Fitting partially linear regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\): x~z$",
    lines
  )
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)

  expect_s3_class(actual$value, "plregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.zero.pos) == 1L)
  expect_true(length(fit.one.pos) >= 1L)
  expect_true(length(fit.two.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.zero.pos[[1L]])
  expect_lt(fit.zero.pos[[1L]], fit.one.pos[[1L]])
  expect_lt(fit.one.pos[[1L]], fit.two.pos[[1L]])
})

test_that("npplreg nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")

  set.seed(20260404)
  n <- 18L
  dat <- data.frame(x = rnorm(n), z = sort(runif(n)))
  dat$y <- 1 + 0.5 * dat$x + sin(2 * pi * dat$z) + rnorm(n, sd = 0.05)

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npplreg(
      y ~ x | z,
      data = dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = npplreg_fit_progress_time_counter()
  )

  lines <- npplreg_fit_progress_lines(actual)
  fit.zero.pos <- grep(
    "^\\[np\\] Fitting partially linear regression 0/2 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting y~z$",
    lines
  )
  fit.one.pos <- grep(
    "^\\[np\\] Fitting partially linear regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): y~z$",
    lines
  )
  fit.two.pos <- grep(
    "^\\[np\\] Fitting partially linear regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\): x~z$",
    lines
  )
  powell.pos <- grep("^\\[np\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[np\\] Selecting degree and bandwidth", lines)

  expect_s3_class(actual$value, "plregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.zero.pos) == 1L)
  expect_true(length(fit.one.pos) >= 1L)
  expect_true(length(fit.two.pos) >= 1L)
  expect_lt(max(powell.pos), fit.zero.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.zero.pos[[1L]])
  expect_lt(fit.zero.pos[[1L]], fit.one.pos[[1L]])
  expect_lt(fit.one.pos[[1L]], fit.two.pos[[1L]])
})
