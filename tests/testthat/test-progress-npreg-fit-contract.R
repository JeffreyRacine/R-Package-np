npreg_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

npreg_fit_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

npreg_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

expect_np_npreg_powell_progress_surface <- function(lines) {
  powell.lines <- lines[grepl("^\\[np\\] Refining bandwidth \\(", lines)]
  info <- paste(lines, collapse = "\n")

  expect_true(length(powell.lines) > 0L, info = info)
  expect_false(any(grepl(
    "^\\[np\\] Refining NOMAD solution with one Powell hot start at degree ",
    lines
  )), info = info)
  expect_false(any(grepl("best (", powell.lines, fixed = TRUE)), info = info)
  expect_true(any(grepl(
    "^\\[np\\] Refining bandwidth \\(elapsed [0-9]+\\.[0-9]s, degree \\([^)]*\\), iter [0-9]+\\)$",
    powell.lines
  )), info = info)
}

make_npreg_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 24L
  dat <- data.frame(x = seq(-0.9, 0.9, length.out = n))
  dat$y <- sin(2 * pi * dat$x) + 0.35 * dat$x

  list(
    dat = dat,
    tx = dat["x"],
    y = dat$y,
    bw = npregbw(
      xdat = dat["x"],
      ydat = dat$y,
      bws = 0.22,
      bandwidth.compute = FALSE
    )
  )
}

test_that("npreg direct bws fit emits known-total fit progress", {
  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)

  expect_s3_class(actual$value, "npregression")
  expect_true(any(grepl(
    "^\\[np\\] Fitting regression 1/24 \\([0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[np\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})

test_that("npreg direct bws fit stays silent below start grace without handoff", {
  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0.75,
    np.progress.interval.known.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npreg_fit_progress_time_values(c(0, 0.2, 0.4, 0.6))
  )

  expect_length(actual$trace, 0L)
})

test_that("npreg bw to fit route hands off immediately into single-line fit progress", {
  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npreg(
      y ~ x,
      data = fixture$dat,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    "^\\[np\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )
  fit.finish.pos <- grep(
    "^\\[np\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )

  expect_s3_class(actual$value, "npregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npreg nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")

  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npreg(
      y ~ x,
      data = fixture$dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[np\\] Selecting degree and bandwidth", lines)
  powell.pos <- grep("^\\[np\\] Refining bandwidth \\(", lines)
  fit.start.pos <- grep(
    "^\\[np\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )
  fit.finish.pos <- grep(
    "^\\[np\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )

  expect_s3_class(actual$value, "npregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_np_npreg_powell_progress_surface(lines)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("predict.npregression re-entry emits fit progress", {
  fixture <- make_npreg_fit_progress_fixture()

  fit <- npreg(
    bws = fixture$bw,
    txdat = fixture$tx,
    tydat = fixture$y
  )

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    predict(fit, newdata = fixture$dat[c(2L, 7L), "x", drop = FALSE]),
    force_renderer = "single_line",
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)

  expect_true(any(grepl(
    "^\\[np\\] Fitting regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[np\\] Fitting regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})
