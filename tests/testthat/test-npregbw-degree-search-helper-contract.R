test_that("internal degree search returns best-so-far on interrupt", {
  degree_search <- getFromNamespace(".np_degree_search", "npRmpi")

  eval_count <- 0L
  result <- degree_search(
    method = "exhaustive",
    candidates = list(0:1),
    baseline_degree = 0L,
    start_degree = 0L,
    eval_fun = function(degree) {
      eval_count <<- eval_count + 1L
      if (eval_count >= 2L) {
        stop(structure(list(message = "interrupt"), class = c("interrupt", "condition")))
      }
      list(
        objective = as.numeric(sum(degree)),
        payload = list(degree = as.integer(degree)),
        num.feval = 1L
      )
    },
    direction = "min",
    trace_level = "full"
  )

  expect_true(isTRUE(result$interrupted))
  expect_false(isTRUE(result$completed))
  expect_identical(as.integer(result$best$degree), 0L)
  expect_identical(as.integer(result$best_payload$degree), 0L)
})

test_that("restart starts are deterministic and RNG-independent", {
  restart_starts <- getFromNamespace(".np_degree_restart_starts", "npRmpi")

  set.seed(1)
  junk1 <- runif(10)
  starts1 <- restart_starts(
    candidates = list(0:2, 0:2),
    restarts = 4L,
    exclude = list(c(0L, 0L))
  )

  set.seed(999)
  junk2 <- runif(10)
  starts2 <- restart_starts(
    candidates = list(0:2, 0:2),
    restarts = 4L,
    exclude = list(c(0L, 0L))
  )

  expect_false(identical(junk1, junk2))
  expect_identical(starts1, starts2)
})

test_that("NOMAD LP degree starts are deterministic, safe, and prefix-stable", {
  build_degree_starts <- getFromNamespace(".np_lp_nomad_build_degree_starts", "npRmpi")

  starts3 <- build_degree_starts(
    initial = c(1L, 1L, 1L, 1L),
    lower = rep(0L, 4L),
    upper = rep(10L, 4L),
    basis = "tensor",
    nobs = 100L,
    nmulti = 3L,
    random.seed = 42L
  )
  starts5 <- build_degree_starts(
    initial = c(1L, 1L, 1L, 1L),
    lower = rep(0L, 4L),
    upper = rep(10L, 4L),
    basis = "tensor",
    nobs = 100L,
    nmulti = 5L,
    random.seed = 42L
  )

  expect_identical(starts5[seq_len(3L), , drop = FALSE], starts3)
  expect_identical(as.integer(starts3[1L, ]), c(1L, 1L, 1L, 1L))
  expect_true(all(apply(starts5, 1L, function(d) npRmpi:::dim_basis(basis = "tensor", degree = d) <= floor(0.25 * (100L - 1L)))))
})

test_that("NOMAD mixed starts preserve user start 1 and expose prefix-stable restart points", {
  build_starts <- getFromNamespace(".np_nomad_build_starts", "npRmpi")

  x0 <- c(1.5, 0.25, 1, 1, 1)
  spec <- list(
    initial = c(1L, 1L, 1L),
    lower = c(0L, 0L, 0L),
    upper = c(10L, 10L, 10L),
    basis = "glp",
    nobs = 80L,
    user_supplied = TRUE
  )

  starts2 <- build_starts(
    x0 = x0,
    bbin = c(0L, 0L, 1L, 1L, 1L),
    lb = c(1e-2, 0, 0, 0, 0),
    ub = c(1e6, 1, 10, 10, 10),
    nmulti = 2L,
    random.seed = 99L,
    degree_spec = spec
  )
  starts4 <- build_starts(
    x0 = x0,
    bbin = c(0L, 0L, 1L, 1L, 1L),
    lb = c(1e-2, 0, 0, 0, 0),
    ub = c(1e6, 1, 10, 10, 10),
    nmulti = 4L,
    random.seed = 99L,
    degree_spec = spec
  )

  expect_equal(starts2[1L, ], x0)
  expect_identical(starts4[seq_len(2L), , drop = FALSE], starts2)
})

test_that("NOMAD Powell hot-start helpers never emit zero public multistarts", {
  hot_nmulti <- getFromNamespace(".np_nomad_powell_hotstart_nmulti", "npRmpi")

  expect_identical(hot_nmulti("disable_multistart"), 1L)
  expect_identical(hot_nmulti("single_iteration"), 1L)
})

test_that("coordinate search skips incumbent cell revisits within a sweep", {
  degree_search <- getFromNamespace(".np_degree_search", "npRmpi")

  result <- degree_search(
    method = "coordinate",
    candidates = list(0:2, 0:2),
    baseline_degree = c(0L, 0L),
    start_degree = c(0L, 0L),
    restarts = 0L,
    max_cycles = 1L,
    eval_fun = function(degree) {
      list(
        objective = as.numeric(sum(degree)),
        payload = list(degree = as.integer(degree)),
        num.feval = 1L
      )
    },
    direction = "min",
    trace_level = "full"
  )

  expect_identical(result$n.unique, 5L)
  expect_identical(result$n.visits, 5L)
  expect_identical(result$n.cached, 0L)
  expect_identical(nrow(result$trace), result$n.unique)
  expect_false(any(result$trace$cached))
})

test_that("automatic exhaustive search emits a safety warning on large grids", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    npRmpi.autodispatch = TRUE,
    np.degree.search.warn.grid = 3L
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(20),
    x2 = runif(20)
  )
  dat$y <- dat$x1 + dat$x2 + rnorm(nrow(dat), sd = 0.05)

  expect_warning(
    npregbw(
      y ~ x1 + x2,
      data = dat,
      regtype = "lp",
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "exhaustive degree search will evaluate 4 degree combinations"
  )
})

test_that("automatic exhaustive search honors internal hard grid limits", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    npRmpi.autodispatch = TRUE,
    np.degree.search.warn.grid = Inf,
    np.degree.search.max.grid = 3L
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(20),
    x2 = runif(20)
  )
  dat$y <- dat$x1 + dat$x2 + rnorm(nrow(dat), sd = 0.05)

  expect_error(
    npregbw(
      y ~ x1 + x2,
      data = dat,
      regtype = "lp",
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "exceeding the configured limit of 3"
  )
})

test_that("automatic degree search emits staged progress output", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = TRUE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(20)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  msgs <- with_npRmpi_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 20, by = 0.5))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("npRmpi"), inherits = FALSE)(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "exhaustive",
        search.engine = "cell",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    )
  )

  expect_true(any(grepl("Automatic polynomial degree search baseline \\(0\\)", msgs)))
  expect_true(any(grepl("Selecting degree and bandwidth", msgs, fixed = TRUE)))
  expect_true(any(grepl("exhaustive", msgs)))
  expect_true(any(grepl("best \\(", msgs)))

  coord_msgs <- with_npRmpi_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 40, by = 0.5))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("npRmpi"), inherits = FALSE)(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "cell",
        degree.min = 0L,
        degree.max = 1L,
        degree.verify = TRUE,
        degree.max.cycles = 2L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    )
  )

  expect_true(any(grepl("Coordinate automatic polynomial degree search over 0:1", coord_msgs)))
  expect_true(any(grepl("max 3 search evaluations", coord_msgs)))
  expect_true(any(grepl("Exhaustively certifying automatic polynomial degree search over 2 degree combinations \\(re-optimizing bandwidths\\)", coord_msgs)))
})
