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

npplreg_fit_progress_targets <- function(xnames = "x") {
  getFromNamespace(".np_plreg_fit_progress_targets", "np")(xnames)
}

npplreg_fit_progress_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

npplreg_fit_progress_pos <- function(lines, done, total, detail) {
  pct <- sprintf("%.1f%%", 100 * done / total)
  grep(
    sprintf(
      "^\\[np\\] Fitting partially linear regression %d/%d \\(%s, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): %s$",
      done,
      total,
      pct,
      npplreg_fit_progress_escape(detail)
    ),
    lines
  )
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

  targets <- npplreg_fit_progress_targets("x")

  expect_identical(state$surface, "bandwidth")
  expect_identical(state$renderer, "single_line")
  expect_equal(state$total, length(targets))
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
  targets <- npplreg_fit_progress_targets("x")
  total <- length(targets)

  expect_s3_class(actual$value, "plregression")
  expect_true(length(npplreg_fit_progress_pos(lines, 1L, total, targets[1L])) >= 1L)
  expect_true(length(npplreg_fit_progress_pos(lines, total, total, targets[total])) >= 1L)
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
  targets <- npplreg_fit_progress_targets("x")
  total <- length(targets)
  fit.zero.pos <- grep(
    sprintf(
      "^\\[np\\] Fitting partially linear regression 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting %s$",
      total,
      npplreg_fit_progress_escape(targets[1L])
    ),
    lines
  )
  fit.one.pos <- npplreg_fit_progress_pos(lines, 1L, total, targets[1L])
  fit.final.pos <- npplreg_fit_progress_pos(lines, total, total, targets[total])
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)

  expect_s3_class(actual$value, "plregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.zero.pos) == 1L)
  expect_true(length(fit.one.pos) >= 1L)
  expect_true(length(fit.final.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.zero.pos[[1L]])
  expect_lt(fit.zero.pos[[1L]], fit.one.pos[[1L]])
  expect_lt(fit.one.pos[[1L]], fit.final.pos[[1L]])
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
  targets <- npplreg_fit_progress_targets("x")
  total <- length(targets)
  fit.zero.pos <- grep(
    sprintf(
      "^\\[np\\] Fitting partially linear regression 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting %s$",
      total,
      npplreg_fit_progress_escape(targets[1L])
    ),
    lines
  )
  fit.one.pos <- npplreg_fit_progress_pos(lines, 1L, total, targets[1L])
  fit.final.pos <- npplreg_fit_progress_pos(lines, total, total, targets[total])
  powell.pos <- grep("^\\[np\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[np\\] (Selecting degree and bandwidth|NOMAD degree/bw|Exhaustive degree/bw|Auto:NOMAD degree/bw|Auto:exhaustive degree/bw)", lines)

  expect_s3_class(actual$value, "plregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.zero.pos) == 1L)
  expect_true(length(fit.one.pos) >= 1L)
  expect_true(length(fit.final.pos) >= 1L)
  if (length(powell.pos) > 0L)
    expect_lt(max(powell.pos), fit.zero.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.zero.pos[[1L]])
  expect_lt(fit.zero.pos[[1L]], fit.one.pos[[1L]])
  expect_lt(fit.one.pos[[1L]], fit.final.pos[[1L]])
})
