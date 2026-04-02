expect_nomad_unknown_bound_progress <- function(msgs, pkg_pattern, info = NULL) {
  expect_false(any(grepl("nomad\\+powell", msgs, ignore.case = TRUE)), info = info)
  expect_false(any(grepl("eval [0-9]+", msgs)), info = info)
  expect_false(any(grepl("fval=", msgs, fixed = TRUE)), info = info)
  expect_false(any(grepl("%|eta ", msgs)), info = info)
  expect_true(any(grepl(sprintf("^\\[%s\\] Selecting degree and bandwidth \\(", pkg_pattern), msgs)), info = info)
  expect_true(any(grepl(sprintf("^\\[%s\\] Refining bandwidth \\(", pkg_pattern), msgs)), info = info)
  expect_true(any(grepl("multistart [12]/2", msgs)), info = info)
  expect_true(any(grepl("iteration [0-9]+", msgs)), info = info)
  expect_true(any(grepl("deg \\(", msgs)), info = info)
  expect_true(any(grepl("best \\(", msgs)), info = info)
}

test_that("remaining serial NOMAD families use unknown-bound restart progress", {
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

  set.seed(20260401)
  n <- 18L
  x <- data.frame(x = sort(runif(n)))
  z <- data.frame(z = sort(runif(n)))
  y_reg <- sin(2 * pi * x$x) + rnorm(n, sd = 0.05)
  y_beta <- rbeta(n, 1 + x$x, 2 - x$x)
  y_pl <- 1 + 0.5 * x$x + cos(2 * pi * z$z) + rnorm(n, sd = 0.05)
  y_sc <- (1 + z$z^2) * x$x + rnorm(n, sd = 0.05)
  y_idx <- sin(x$x + 0.5 * z$z) + rnorm(n, sd = 0.05)

  cases <- list(
    npcdistbw = function() np::npcdistbw(
      y_reg ~ x,
      data = data.frame(y_reg = y_reg, x = x$x),
      regtype = "lp",
      degree.select = "coordinate",
      search.engine = "nomad+powell",
      degree.min = 0L,
      degree.max = 1L,
      degree.verify = FALSE,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 2L,
      ngrid = 30L,
      max.bb.eval = 8L
    ),
    npplregbw = function() np::npplregbw(
      xdat = x,
      zdat = z,
      ydat = y_pl,
      regtype = "lp",
      bernstein.basis = TRUE,
      degree.select = "coordinate",
      search.engine = "nomad+powell",
      degree.min = 0L,
      degree.max = 1L,
      degree.verify = FALSE,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 2L,
      max.bb.eval = 8L
    ),
    npscoefbw = function() np::npscoefbw(
      xdat = x,
      zdat = z,
      ydat = y_sc,
      regtype = "lp",
      bernstein.basis = TRUE,
      degree.select = "coordinate",
      search.engine = "nomad+powell",
      degree.min = 0L,
      degree.max = 1L,
      degree.verify = FALSE,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 2L,
      max.bb.eval = 8L
    ),
    npindexbw = function() np::npindexbw(
      xdat = data.frame(x1 = x$x, x2 = z$z),
      ydat = y_idx,
      bws = c(1, 0.5, 0.35),
      method = "ichimura",
      regtype = "lp",
      bernstein.basis = TRUE,
      degree.select = "coordinate",
      search.engine = "nomad+powell",
      degree.min = 0L,
      degree.max = 1L,
      degree.verify = FALSE,
      bwtype = "fixed",
      nmulti = 2L,
      max.bb.eval = 8L
    )
  )

  for (case in names(cases)) {
    msgs <- with_np_degree_bindings(
      list(
        .np_progress_is_interactive = function() TRUE,
        .np_progress_renderer_for_surface = function(surface, capability) "legacy",
        .np_progress_now = degree_progress_time_values(seq(0, 60, by = 0.25))
      ),
      capture_degree_messages_only(cases[[case]]())
    )
    expect_nomad_unknown_bound_progress(msgs, pkg_pattern = "np", info = case)
  }
})
