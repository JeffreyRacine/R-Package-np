condensdist_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

condensdist_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

expect_condensdist_clean_powell_surface <- function(lines, pkg_pattern = "np", info = NULL) {
  powell.lines <- grep(sprintf("^\\[%s\\] Refining bandwidth \\(", pkg_pattern), lines, value = TRUE)
  detail.info <- if (!is.null(info)) info else paste(lines, collapse = "\n")

  expect_true(length(powell.lines) > 0L, info = detail.info)
  expect_true(
    any(grepl(
      sprintf("^\\[%s\\] Refining bandwidth \\(elapsed [0-9]+\\.[0-9]s, degree \\([0-9]+\\)(, iter [0-9]+)?\\)$", pkg_pattern),
      powell.lines
    )),
    info = paste(c(detail.info, powell.lines), collapse = "\n")
  )
  expect_false(any(grepl("best \\(", powell.lines)), info = paste(c(detail.info, powell.lines), collapse = "\n"))
  expect_false(
    any(grepl(sprintf("^\\[%s\\] Bandwidth selection \\(Refining NOMAD solution", pkg_pattern), lines)),
    info = detail.info
  )
}

make_condensdist_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 18L
  dat <- data.frame(
    x = seq(-0.9, 0.9, length.out = n),
    y = 0.4 * seq(-0.9, 0.9, length.out = n) + sin(seq(-0.9, 0.9, length.out = n))
  )

  list(
    dat = dat,
    tx = dat["x"],
    ty = dat["y"],
    n = n,
    cdens.bw = npcdensbw(
      xdat = dat["x"],
      ydat = dat["y"],
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE
    ),
    cdist.bw = npcdistbw(
      xdat = dat["x"],
      ydat = dat["y"],
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE
    )
  )
}

test_that("npcdens direct bws fit emits single-line fit progress", {
  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npcdens(
      bws = fixture$cdens.bw,
      txdat = fixture$tx,
      tydat = fixture$ty
    ),
    force_renderer = "single_line",
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "condensity")
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting conditional density 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", total),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting conditional density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )))
})

test_that("npcdens bw to fit route hands off immediately into single-line fit progress", {
  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npcdens(
      txdat = fixture$tx,
      tydat = fixture$ty,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional density 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", total),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )

  expect_s3_class(actual$value, "condensity")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npcdens nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")

  fixture <- make_condensdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npcdens(
      y ~ x,
      data = fixture$dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  fit.start.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional density 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )
  powell.pos <- grep("^\\[np\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[np\\] Selecting degree and bandwidth", lines)

  expect_s3_class(actual$value, "condensity")
  expect_condensdist_clean_powell_surface(lines, pkg_pattern = "np")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npcdist direct bws fit emits single-line fit progress", {
  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npcdist(
      bws = fixture$cdist.bw,
      txdat = fixture$tx,
      tydat = fixture$ty
    ),
    force_renderer = "single_line",
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "condistribution")
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting conditional distribution 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", total),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[np\\] Fitting conditional distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )))
})

test_that("npcdist bw to fit route hands off immediately into single-line fit progress", {
  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npcdist(
      txdat = fixture$tx,
      tydat = fixture$ty,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[np\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional distribution 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", total),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )

  expect_s3_class(actual$value, "condistribution")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npcdist nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")

  fixture <- make_condensdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npcdist(
      y ~ x,
      data = fixture$dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  fit.start.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional distribution 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[np\\] Fitting conditional distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )
  powell.pos <- grep("^\\[np\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[np\\] Selecting degree and bandwidth", lines)

  expect_s3_class(actual$value, "condistribution")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})
