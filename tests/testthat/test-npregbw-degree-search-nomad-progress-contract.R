test_that("npregbw NOMAD plus Powell progress keeps stage text on managed lines", {
  skip_if_not_installed("crs")

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260321)
  dat <- data.frame(x = sort(runif(16)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  msgs <- with_np_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 40, by = 0.25))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("np"), inherits = FALSE)(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "nomad+powell",
        degree.min = 0L,
        degree.max = 1L,
        degree.verify = FALSE,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L,
        max.bb.eval = 8
      )
    )
  )

  expect_false(any(grepl("^\\[np\\] Starting NOMAD automatic polynomial degree search at degree ", msgs)))
  expect_false(any(grepl("^\\[np\\] Refining NOMAD solution with one Powell hot start at degree ", msgs)))
  expect_true(any(grepl("starting at degree \\(", msgs)))
})
