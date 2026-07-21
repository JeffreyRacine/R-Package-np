capture_progress_shadow_with_conditions <- function(expr, force_renderer = NULL, now = function() 0) {
  messages <- character()
  warnings <- character()

  shadow <- withCallingHandlers(
    capture_progress_shadow_trace(expr, force_renderer = force_renderer, now = now),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    value = shadow$value,
    trace = shadow$trace,
    final_line = shadow$final_line,
    messages = sub("\n$", "", messages),
    warnings = warnings
  )
}

skip_live_route_slice <- function() {
  skip_on_cran()
}

progress_time_counter <- function(start = 0, by = 2.1) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

make_iv_data <- function(n = 40) {
  set.seed(42)
  z <- runif(n)
  w <- z + rnorm(n, sd = 0.1)
  y <- z^2 + rnorm(n, sd = 0.1)
  list(y = y, z = z, w = w)
}

shadow_lines_matching <- function(shadow, pattern) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  lines[grepl(pattern, lines)]
}

shadow_signature <- function(shadow, pattern) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  events <- vapply(shadow$trace, `[[`, character(1L), "event")
  keep <- grepl(pattern, lines) & events == "render"

  data.frame(
    event = events[keep],
    line = lines[keep],
    stringsAsFactors = FALSE
  )
}

test_that("Landweber npregiv single-line progress reports object labels with outer iterations", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  single_line <- capture_progress_shadow_with_conditions(
    npregiv(y = dat$y, z = dat$z, w = dat$w, method = "Landweber-Fridman", iterate.max = 4),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- vapply(single_line$trace, `[[`, character(1L), "line")

  expect_s3_class(single_line$value, "npregiv")
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[y\\|w\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[[^)]*\\|z\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[y-phi\\(z\\)\\|w\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[y-phi\\(z\\)\\|w\\]\\|z\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_false(any(grepl("Iterating Landweber-Fridman solve", lines, fixed = TRUE)))
  expect_false(any(grepl("%|eta ", lines)))
})

test_that("Tikhonov npregiv single-line progress restores historical object labels", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  single_line <- capture_progress_shadow_with_conditions(
    npregiv(
      y = dat$y,
      z = dat$z,
      w = dat$w,
      method = "Tikhonov",
      iterate.Tikhonov = TRUE,
      iterate.Tikhonov.num = 3
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- vapply(single_line$trace, `[[`, character(1L), "line")

  expect_s3_class(single_line$value, "npregiv")
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[y\\|w\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[y\\|w\\]\\|z\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(alpha, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[phi\\(z\\)\\|w\\], iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[phi\\(z\\)\\|w\\]\\|z\\], iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(phi\\(z\\), iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_false(any(grepl("Iterating Tikhonov solve", lines, fixed = TRUE)))
  expect_false(any(grepl("%|eta ", lines)))
})

test_that("Landweber npregiv without residual smoothing reports the alternate historical objects", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  single_line <- capture_progress_shadow_with_conditions(
    npregiv(
      y = dat$y,
      z = dat$z,
      w = dat$w,
      method = "Landweber-Fridman",
      smooth.residuals = FALSE,
      iterate.max = 4
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- vapply(single_line$trace, `[[`, character(1L), "line")

  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[phi\\(z\\)\\|w\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[y\\|w\\]-E\\[phi\\(z\\)\\|w\\]\\|z\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("npregivderiv single-line progress matches legacy semantics", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_with_conditions(
    npregivderiv(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  dat <- make_iv_data()
  single_line <- capture_progress_shadow_with_conditions(
    npregivderiv(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  pattern <- "^\\[npRmpi\\] IV derivative"
  lines <- shadow_lines_matching(single_line, pattern)

  expect_s3_class(single_line$value, "npregivderiv")
  expect_equal(shadow_signature(single_line, pattern), shadow_signature(legacy, pattern))
  expect_true(any(grepl("^\\[npRmpi\\] Preparing f_Z\\(z\\), S_Z\\(z\\)$", single_line$messages)))
  expect_true(any(grepl("^\\[npRmpi\\] IV derivative\\.\\.\\. elapsed [0-9]+\\.[0-9]s: E\\[y\\|w\\]$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV derivative\\.\\.\\. elapsed [0-9]+\\.[0-9]s: d/dz E\\[y\\|z\\]$", lines)))
  expect_true(any(grepl("E\\[y-phi_0\\(z\\)\\|w\\]$", lines)))
  expect_true(any(grepl("T\\*\\{E\\[y-phi_0\\(z\\)\\|w\\]\\}$", lines)))
  expect_true(any(grepl("iteration 1, elapsed [0-9]+\\.[0-9]s: E\\[y-phi\\(z\\)\\|w\\]$", lines)))
  expect_true(any(grepl("iteration 1, elapsed [0-9]+\\.[0-9]s: T\\*\\{E\\[y-phi\\(z\\)\\|w\\]\\}$", lines)))
  expect_false(any(grepl("iteration 0", lines, fixed=TRUE)))
  expect_false(any(grepl("iteration 3, .*T\\*", lines)))
  expect_false(any(grepl("Iterating Landweber-Fridman derivative solve|updating phi|E\\(mu\\|w\\)", lines)))
})

test_that("npregivderiv MPI progress labels alternate initialization routes truthfully", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages=TRUE)
  on.exit(options(old_opts), add=TRUE)

  capture_route <- function(...) {
    capture_progress_shadow_with_conditions(
      npregivderiv(y=dat$y, z=dat$z, w=dat$w,
                   nmulti=1L, iterate.max=2L, ...),
      force_renderer="single_line",
      now=progress_time_counter()
    )
  }
  route_lines <- function(x) {
    lines <- vapply(x$trace, `[[`, character(1L), "line")
    lines[grepl("^\\[npRmpi\\] IV derivative", lines)]
  }

  component <- route_lines(capture_route(smooth.residuals=FALSE))
  expect_true(any(grepl("E\\[phi_0\\(z\\)\\|w\\]$", component)))
  expect_true(any(grepl("iteration 1, .*E\\[phi\\(z\\)\\|w\\]$", component)))
  expect_true(any(grepl("T\\*\\{E\\[y-phi\\(z\\)\\|w\\]\\}$", component)))

  nested <- route_lines(capture_route(start.from="EEywz"))
  expect_true(any(grepl("d/dz E\\[E\\[y\\|w\\]\\|z\\]$", nested)))

  supplied <- route_lines(capture_route(starting.values=rep(0, length(dat$y))))
  expect_false(any(grepl("d/dz", supplied, fixed=TRUE)))
})

test_that("npregivderiv MPI mathematical labels fit ordinary terminal budgets", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "npRmpi")
  label <- "[npRmpi] IV derivative... iteration 1000, elapsed 1234.5s: T*{E[y-phi(z)|w]}"

  for(width in c(76L, 80L, 116L, 120L)) {
    fitted <- fit(label, max_width=width)
    expect_lte(nchar(fitted, type="width"), width)
    expect_match(fitted, "iteration 1000", fixed=TRUE)
    expect_match(fitted, "elapsed 1234.5s", fixed=TRUE)
    expect_match(fitted, "T*{E[y-phi(z)|w]}", fixed=TRUE)
  }
})

test_that("npregivderiv MPI progress preserves disabled and absent message options", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)

  dat <- make_iv_data()
  original <- options()
  had_option <- "np.messages" %in% names(original)
  original_value <- original[["np.messages"]]
  on.exit({
    if(had_option) options(np.messages=original_value) else options(np.messages=NULL)
  }, add=TRUE)

  options(np.messages=FALSE)
  disabled <- capture_progress_shadow_with_conditions(
    npregivderiv(y=dat$y, z=dat$z, w=dat$w,
                 nmulti=1L, iterate.max=2L),
    now=progress_time_counter()
  )
  expect_length(disabled$trace, 0L)
  expect_identical(getOption("np.messages"), FALSE)

  options(np.messages=NULL)
  absent <- suppressWarnings(suppressMessages(npregivderiv(
    y=dat$y, z=dat$z, w=dat$w, nmulti=1L, iterate.max=2L
  )))
  expect_s3_class(absent, "npregivderiv")
  expect_false("np.messages" %in% names(options()))
})

test_that("npregivderiv MPI releases its owner after a rank-symmetric error", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages=TRUE)
  on.exit(options(old_opts), add=TRUE)
  registry <- getFromNamespace(".np_progress_registry", "npRmpi")

  expect_error(
    capture_progress_shadow_with_conditions(
      npregivderiv(y=dat$y, z=dat$z, w=dat$w,
                   bwmethod="invalid", nmulti=1L, iterate.max=2L),
      force_renderer="single_line",
      now=progress_time_counter()
    ),
    "invalid|arg"
  )
  expect_null(registry$active_id)

  reuse <- suppressWarnings(suppressMessages(npregivderiv(
    y=dat$y, z=dat$z, w=dat$w, nmulti=1L, iterate.max=2L
  )))
  expect_s3_class(reuse, "npregivderiv")
  expect_null(registry$active_id)
})

test_that("npregiv proof-slice progress respects np.messages FALSE", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_with_conditions(
    npregiv(y = dat$y, z = dat$z, w = dat$w, method = "Landweber-Fridman", iterate.max = 4),
    now = progress_time_counter()
  )

  expect_length(res$messages, 0)
  expect_length(res$trace, 0)
})

test_that("npregivderiv proof-slice progress respects suppressMessages", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_with_conditions(
    suppressMessages(npregivderiv(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3)),
    now = progress_time_counter()
  )

  expect_length(res$messages, 0)
  expect_length(res$trace, 0)
})

test_that("IV source routes use the canonical single-line surface tag", {
  src <- installed_function_text("npregiv.default")
  expect_true(grepl(".np_progress_select_iv(", src, fixed = TRUE))

  src_deriv <- installed_function_text("npregivderiv.default")
  expect_true(grepl('surface = "iv_solve"', src_deriv, fixed = TRUE))
})
