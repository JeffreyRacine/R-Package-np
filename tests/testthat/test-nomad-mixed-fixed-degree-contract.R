test_that("NOMAD mixed fixed/free degree geometry is isolated", {
  pkg <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  geometry_for <- getFromNamespace(".np_nomad_fixed_degree_geometry", pkg)
  add_option <- getFromNamespace(".np_nomad_add_fixed_degree_option", pkg)

  original_options <- list(MAX_BB_EVAL = 40L)
  all_free <- geometry_for(
    degree_spec = list(lower = c(0L, 0L), upper = c(1L, 1L)),
    ncoordinate = 4L,
    bbin = c(0L, 0L, 1L, 1L),
    lb = c(0.1, 0.1, 0, 0),
    ub = c(1, 1, 1, 1),
    start_matrix = rbind(c(0.5, 0.5, 0, 0))
  )
  expect_null(all_free)
  expect_identical(add_option(original_options, all_free), original_options)

  mixed <- geometry_for(
    degree_spec = list(lower = c(0L, 2L), upper = c(1L, 2L)),
    ncoordinate = 4L,
    bbin = c(0L, 0L, 1L, 1L),
    lb = c(0.1, 0.1, 0, 2),
    ub = c(1, 1, 1, 2),
    start_matrix = rbind(c(0.5, 0.5, 0, 2))
  )
  expect_identical(mixed$fixed.index, 4L)
  expect_identical(mixed$fixed.value, 2L)
  expect_equal(mixed$lower, c(0.1, 0.1, 0, 1), tolerance = 0)
  expect_equal(mixed$upper, c(1, 1, 1, 2), tolerance = 0)
  expect_identical(mixed$option.value, "( - - - 2 )")
  expect_identical(
    add_option(original_options, mixed)$FIXED_VARIABLE,
    mixed$option.value
  )
  assert_solution <- getFromNamespace(
    ".np_nomad_assert_fixed_degree_solution", pkg
  )
  expect_true(assert_solution(c(0.5, 0.5, 0, 2), mixed))
  expect_error(
    assert_solution(c(0.5, 0.5, 0, 1), mixed, "synthetic NOMAD"),
    "coordinate 4: requested 2, returned 1",
    fixed = TRUE
  )
  expect_error(
    add_option(list(` fixed_variable ` = "( - - - 2 )"), mixed),
    "'FIXED_VARIABLE' is reserved internally",
    fixed = TRUE
  )
})

test_that("R-bridge NOMAD search asserts mixed fixed degrees", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")

  result <- .np_nomad_search(
    engine = "nomad",
    baseline_record = NULL,
    start_degree = c(1L, 0L),
    x0 = c(1, 0),
    bbin = c(1L, 1L),
    lb = c(1, 0),
    ub = c(1, 1),
    eval_fun = function(point) list(
      objective = (point[2L] - 1)^2,
      degree = as.integer(round(point)),
      num.feval = 1L
    ),
    build_payload = function(point, best_record, solution, interrupted)
      list(payload = list(point = point)),
    nmulti = 1L,
    degree_spec = list(
      initial = c(1L, 0L), lower = c(1L, 0L), upper = c(1L, 1L),
      basis = "glp", nobs = 20L, user_supplied = TRUE
    ),
    nomad.opts = list(MAX_BB_EVAL = 12L),
    native.r.bridge = TRUE
  )

  expect_identical(result$best_point[1L], 1)
  expect_true(all(vapply(result$restart.results, function(x)
    identical(as.numeric(x$solution[1L]), 1), logical(1L))))
})

test_that("prepared NOMAD search honors a positive fixed degree", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  i <- seq_len(32L)
  dat <- data.frame(
    x1 = 2 * ((i * sqrt(2)) %% 1) - 1,
    x2 = 2 * ((i * sqrt(3)) %% 1) - 1
  )
  dat$y <- sin(dat$x1) + 0.4 * dat$x2 + cos(i) / 30

  bw <- npregbw(
    y ~ x1 + x2,
    data = dat,
    nomad = TRUE,
    search.engine = "nomad",
    degree.min = c(1L, 0L),
    degree.max = c(1L, 1L),
    degree.start = c(1L, 0L),
    bwtype = "fixed",
    nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 40L)
  )

  trace_degree <- do.call(
    rbind,
    strsplit(as.character(bw$degree.search$trace$degree), ",", fixed = TRUE)
  )
  expect_s3_class(bw, "rbandwidth")
  expect_identical(as.integer(bw$degree[1L]), 1L)
  expect_true(all(trace_degree[, 1L] == "1"))
  expect_true(is.finite(bw$fval))
})
