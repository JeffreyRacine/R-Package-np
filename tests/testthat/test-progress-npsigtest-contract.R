progress_time_counter <- function(start = 0, by = 0.6) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_bootstrap_signature <- function(shadow) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  events <- vapply(shadow$trace, `[[`, character(1L), "event")
  keep <- grepl("^\\[npRmpi\\] Bootstrap replications", lines)

  data.frame(
    event = events[keep],
    line = lines[keep],
    stringsAsFactors = FALSE
  )
}

shadow_lines <- function(shadow) {
  shadow_bootstrap_signature(shadow)$line
}

skip_live_route_slice <- function() {
  skip_if_not(
    identical(Sys.getenv("NP_RMPI_PROGRESS_LIVE_ROUTE_TESTS", ""), "true"),
    "live npRmpi route slice is gated to manual session/attach/profile proof artifacts"
  )
}

npsigtest_fun <- function(...) {
  getFromNamespace("npsigtest", "npRmpi")(...)
}

make_sigtest_fixture <- function(seed = 42, n = 30) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.1)
  bw <- getFromNamespace("npregbw", "npRmpi")(
    y ~ x1 + x2,
    bws = c(0.2, 0.4),
    bandwidth.compute = FALSE
  )
  list(bw = bw)
}

test_that("npsigtest joint single-line bootstrap progress matches legacy semantics", {
  skip_on_cran()
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_sigtest_fixture()

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = TRUE, index = 1),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  single_line <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = TRUE, index = 1),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)

  expect_s3_class(single_line$value, "sigtest")
  expect_equal(shadow_bootstrap_signature(single_line), shadow_bootstrap_signature(legacy))
  expect_true(any(grepl("^\\[npRmpi\\] Bootstrap replications [0-9]+/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("npsigtest individual single-line bootstrap progress matches legacy semantics", {
  skip_on_cran()
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_sigtest_fixture(seed = 99)

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = FALSE, index = c(1, 2)),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  single_line <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = FALSE, index = c(1, 2)),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)

  expect_s3_class(single_line$value, "sigtest")
  expect_equal(shadow_bootstrap_signature(single_line), shadow_bootstrap_signature(legacy))
  expect_true(sum(grepl("^\\[npRmpi\\] Bootstrap replications [0-9]+/9 ", lines)) >= 2L)
  expect_true(sum(grepl("^\\[npRmpi\\] Bootstrap replications 9/9 ", lines)) >= 2L)
})

test_that("npsigtest progress respects np.messages FALSE", {
  skip_on_cran()
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_sigtest_fixture()

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = TRUE, index = 1),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})

test_that("npsigtest progress respects suppressMessages", {
  skip_on_cran()
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_sigtest_fixture()

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    suppressMessages(npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = TRUE, index = 1)),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})

test_that("npsigtest source routes use canonical bootstrap surface tags", {
  src <- paste(
    installed_function_text("npsigtest.rbandwidth"),
    installed_function_text("npsigtest.default"),
    sep = "\n"
  )
  expect_true(grepl('surface = "bootstrap"', src, fixed = TRUE))
})
