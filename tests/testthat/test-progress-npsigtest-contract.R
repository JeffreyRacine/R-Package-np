progress_time_counter <- function(start = 0, by = 0.6) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_npsigtest_signature <- function(shadow) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  events <- vapply(shadow$trace, `[[`, character(1L), "event")
  keep <- grepl("^\\[np\\] Testing ", lines)

  data.frame(
    id = vapply(shadow$trace, `[[`, character(1L), "id")[keep],
    event = events[keep],
    line = lines[keep],
    stringsAsFactors = FALSE
  )
}

shadow_lines <- function(shadow) {
  shadow_npsigtest_signature(shadow)$line
}

npsigtest_fun <- function(...) {
  getFromNamespace("npsigtest", "np")(...)
}

make_sigtest_fixture <- function(seed = 42, n = 30) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.1)
  bw <- getFromNamespace("npregbw", "np")(
    y ~ x1 + x2,
    bws = c(0.2, 0.4),
    bandwidth.compute = FALSE
  )
  list(bw = bw)
}

test_that("npsigtest joint progress delays bootstrap ETA until work completes", {
  fixture <- make_sigtest_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, B = 9, joint = TRUE, index = 1),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  single_line <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, B = 9, joint = TRUE, index = 1),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)
  legacy.signature <- shadow_npsigtest_signature(legacy)
  single.signature <- shadow_npsigtest_signature(single_line)

  expect_s3_class(single_line$value, "sigtest")
  expect_equal(
    single.signature[single.signature$event != "finish", ],
    legacy.signature[legacy.signature$event != "finish", ]
  )
  expect_match(lines[[1L]], "^\\[np\\] Testing joint significance\\.\\.\\. elapsed 0\\.0s$")
  expect_false(grepl("eta", lines[[1L]], fixed = TRUE))
  expect_true(any(grepl("^\\[np\\] Testing joint significance 1/9 ", lines)))
  expect_true(any(grepl("^\\[np\\] Testing joint significance 9/9 ", lines)))
  expect_true(any(grepl("eta [0-9]+\\.[0-9]s", lines)))
  expect_length(unique(single.signature$id), 1L)
  expect_identical(tail(single.signature$event, 1L), "finish")
})

test_that("npsigtest individual progress uses completed predictors for ETA", {
  fixture <- make_sigtest_fixture(seed = 99)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, B = 9, joint = FALSE, index = c(1, 2)),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  single_line <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, B = 9, joint = FALSE, index = c(1, 2)),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)
  legacy.signature <- shadow_npsigtest_signature(legacy)
  single.signature <- shadow_npsigtest_signature(single_line)

  expect_s3_class(single_line$value, "sigtest")
  expect_equal(
    single.signature[single.signature$event != "finish", ],
    legacy.signature[legacy.signature$event != "finish", ]
  )
  expect_match(lines[[1L]], "^\\[np\\] Testing x1\\.\\.\\. elapsed 0\\.0s$")
  expect_false(grepl("eta", lines[[1L]], fixed = TRUE))
  expect_true(any(grepl("^\\[np\\] Testing x2 1/2 ", lines)))
  expect_true(any(grepl("^\\[np\\] Testing x2 2/2 ", lines)))
  expect_false(any(grepl("/[9] ", lines)))
  expect_false(any(grepl("of \\(1,2\\)", lines)))
  expect_length(unique(single.signature$id), 1L)
  expect_identical(tail(single.signature$event, 1L), "finish")
})

test_that("npsigtest progress respects np.messages FALSE", {
  fixture <- make_sigtest_fixture(seed = 17)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, B = 9, joint = TRUE, index = 1),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})

test_that("npsigtest progress respects suppressMessages", {
  fixture <- make_sigtest_fixture(seed = 23)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    suppressMessages(npsigtest_fun(bws = fixture$bw, B = 9, joint = TRUE, index = 1)),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})
