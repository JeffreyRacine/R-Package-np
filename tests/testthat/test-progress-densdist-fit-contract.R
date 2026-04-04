densdist_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

densdist_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

make_densdist_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 24L
  dat <- data.frame(x = seq(-0.9, 0.9, length.out = n))

  list(
    dat = dat,
    tdat = dat["x"],
    n = n,
    dens.bw = npudensbw(
      dat = dat["x"],
      bws = 0.24,
      bandwidth.compute = FALSE
    ),
    dist.bw = npudistbw(
      dat = dat["x"],
      bws = 0.24,
      bandwidth.compute = FALSE
    )
  )
}

test_that("npudens direct bws fit emits single-line fit progress", {
  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npudens(
      bws = fixture$dens.bw,
      tdat = fixture$tdat
    ),
    force_renderer = "single_line",
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "npdensity")
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting density 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", fixture$n),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )))
})

test_that("npudens bw to fit route hands off immediately into single-line fit progress", {
  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npudens(
      tdat = fixture$tdat,
      nmulti = 2L
    ),
    force_renderer = "single_line",
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[np\\] Fitting density 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[np\\] Fitting density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )

  expect_s3_class(actual$value, "npdensity")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npudist direct bws fit emits single-line fit progress", {
  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npudist(
      bws = fixture$dist.bw,
      tdat = fixture$tdat
    ),
    force_renderer = "single_line",
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "npdistribution")
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting distribution 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", fixture$n),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )))
})

test_that("npudist bw to fit route hands off immediately into single-line fit progress", {
  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npudist(
      tdat = fixture$tdat,
      nmulti = 2L
    ),
    force_renderer = "single_line",
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[np\\] Fitting distribution 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[np\\] Fitting distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )

  expect_s3_class(actual$value, "npdistribution")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})
