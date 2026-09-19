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
  expect_match(lines[[1L]], "rep 0/9, elapsed 0\\.0s, eta estimating")
  expect_true(any(grepl("rep 9/9", lines, fixed = TRUE)))
  expect_false(any(grepl("target", lines, fixed = TRUE)))
  expect_true(any(grepl("eta 0s", lines, fixed = TRUE)))
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
  expect_match(lines[[1L]], "Testing x1 \\(target 1/2, rep 0/9")
  expect_match(lines[[1L]], "eta estimating")
  expect_true(any(grepl("Testing x2 (target 2/2, rep 0/9", lines, fixed = TRUE)))
  expect_true(any(grepl("Testing x2 (target 2/2, rep 9/9", lines, fixed = TRUE)))
  expect_false(any(grepl("Testing x2 (target 1/2", lines, fixed = TRUE)))
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
# np active-target/replication state boundaries.

test_that("npsigtest reordered targets keep their requested ordinal", {
  fixture <- make_sigtest_fixture(seed = 109)
  old <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old), add = TRUE)
  result <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, B = 9, index = c(2, 1)),
    now = progress_time_counter())
  lines <- shadow_lines(result)
  expect_match(lines[[1L]], "Testing x2 \\(target 1/2, rep 0/9")
  expect_true(any(grepl("Testing x1 (target 2/2, rep 0/9", lines, fixed = TRUE)))
})

test_that("npsigtest progress separates active ordinal and fractional work", {
  fields <- getFromNamespace(".np_npsig_progress_fields", "np")
  state <- new.env(parent = emptyenv())
  state$started <- 0
  state$npsig.completed <- 0L
  state$npsig.target <- 1L
  state$npsig.total <- 3L
  state$npsig.B <- 10L
  state$npsig.joint <- FALSE
  state$npsig.skipped <- FALSE
  expect_identical(fields(state, 0L, NULL, 0),
    c("target 1/3", "rep 0/10", "elapsed 0.0s", "eta estimating"))
  expect_identical(fields(state, 5L, NULL, 2),
    c("target 1/3", "rep 5/10", "elapsed 2.0s", "eta 10.0s"))
  state$npsig.completed <- 1L
  state$npsig.target <- 2L
  expect_identical(fields(state, 0L, NULL, 4),
    c("target 2/3", "rep 0/10", "elapsed 4.0s", "eta 8.0s"))
  state$npsig.completed <- 2L
  state$npsig.skipped <- TRUE
  expect_identical(fields(state, 0L, NULL, 4),
    c("target 2/3", "rep skipped", "elapsed 4.0s", "eta 2.0s"))
  state$npsig.completed <- 3L
  state$npsig.target <- 3L
  expect_identical(tail(fields(state, 0L, NULL, 4), 1), "eta 0s")
  compact <- getFromNamespace(".np_progress_fit_single_line", "np")
  line <- "[np] Testing nonwhite (target 3/6, rep 145/999, elapsed 69.8s, eta 139.6s)"
  for (width in c(40L, 60L, 80L)) {
    out <- compact(line, width)
    expect_lte(nchar(out, type = "width"), width)
    expect_match(out, "(target|tgt) 3/6")
  }
})

test_that("analytic non-rejection does not report simulated draws", {
  set.seed(133)
  d <- data.frame(z = factor(rep(0:1, 15)))
  y <- rnorm(30)
  bw <- getFromNamespace("npregbw", "np")(
    d, y, bws = .5, bandwidth.compute = FALSE)
  old <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old), add = TRUE)
  result <- capture_progress_shadow_trace(
    npsigtest_fun(bws = bw, B = 9), now = progress_time_counter())
  lines <- shadow_lines(result)
  expect_true(any(grepl("rep skipped", lines, fixed = TRUE)))
  expect_false(any(grepl("rep 9/9", lines, fixed = TRUE)))
  expect_identical(unname(result$value$P), 1)
  expect_true(all(is.na(result$value$In.bootstrap)))
  expect_identical(tail(shadow_npsigtest_signature(result)$event, 1L), "finish")
})

test_that("npsigtest restores the native heartbeat owner after an error", {
  fixture <- make_sigtest_fixture(seed = 219)
  runtime <- getFromNamespace(".np_progress_runtime", "np")
  previous <- runtime$fit_forward
  old <- options(np.messages = TRUE)
  on.exit(options(old), add = TRUE)
  result <- with_np_progress_bindings(
    list(.npreg_complete = function(...) stop("progress cleanup sentinel")),
    capture_progress_shadow_trace(
      tryCatch(npsigtest_fun(bws = fixture$bw, B = 9), error = identity),
      now = progress_time_counter()))
  expect_match(conditionMessage(result$value), "progress cleanup sentinel")
  expect_identical(runtime$fit_forward, previous)
  expect_length(unique(vapply(result$trace, `[[`, character(1), "id")), 1L)
  expect_identical(tail(shadow_npsigtest_signature(result)$event, 1L), "abort")
})
