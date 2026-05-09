test_that("internal degree search returns best-so-far on interrupt", {
  degree_search <- getFromNamespace(".np_degree_search", "np")

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
  restart_starts <- getFromNamespace(".np_degree_restart_starts", "np")

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
  build_degree_starts <- getFromNamespace(".np_lp_nomad_build_degree_starts", "np")

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
  expect_true(all(apply(starts5, 1L, function(d) np:::dim_basis(basis = "tensor", degree = d) <= floor(0.25 * (100L - 1L)))))
})

test_that("NOMAD mixed starts preserve user start 1 and expose prefix-stable restart points", {
  build_starts <- getFromNamespace(".np_nomad_build_starts", "np")

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
  hot_nmulti <- getFromNamespace(".np_nomad_powell_hotstart_nmulti", "np")

  expect_identical(hot_nmulti("disable_multistart"), 1L)
  expect_identical(hot_nmulti("single_iteration"), 1L)
})

test_that("NOMAD Powell hot-start option helper controls internal remin explicitly", {
  hot_opt_args <- getFromNamespace(".np_nomad_powell_hotstart_opt_args", "np")

  opt.args <- list(
    nmulti = 7L,
    powell.remin = TRUE,
    ftol = 1e-6,
    tol = 1e-5,
    custom = "preserve-me"
  )

  out <- hot_opt_args(opt.args, strategy = "disable_multistart")

  expect_identical(out$nmulti, 1L)
  expect_false(out$powell.remin)
  expect_identical(out$ftol, opt.args$ftol)
  expect_identical(out$tol, opt.args$tol)
  expect_identical(out$custom, opt.args$custom)

  out2 <- hot_opt_args(opt.args, strategy = "single_iteration")
  expect_identical(out2$nmulti, 1L)
  expect_false(out2$powell.remin)

  out3 <- hot_opt_args(opt.args, strategy = "single_iteration", remin = TRUE)
  expect_identical(out3$nmulti, 1L)
  expect_true(out3$powell.remin)
})

test_that("NOMAD restart summary records remin metadata", {
  attach_summary <- getFromNamespace(".np_attach_nomad_restart_summary", "np")

  bws <- list(fval = 1)
  search_result <- list(
    direction = "min",
    best = list(objective = 0.5),
    best.restart = 2L,
    restart.results = list(
      list(objective = 1.0),
      list(objective = 0.5, remin = TRUE)
    ),
    restart.starts = list(10, 20),
    restart.degree.starts = list(1L, 2L),
    restart.bandwidth.starts = list(0.1, 0.2),
    nomad.remin = TRUE,
    nomad.remin.index = 2L,
    nomad.remin.roundtrip = list(objective = 1.0, degree = 1L)
  )

  out <- attach_summary(bws, search_result)

  expect_true(out$nomad.remin)
  expect_identical(out$nomad.remin.index, 2L)
  expect_identical(out$nomad.remin.roundtrip$degree, 1L)
  expect_identical(out$nomad.best.restart, 2L)
  expect_equal(out$nomad.restart.fval, c(1.0, 0.5))
})

test_that("coordinate search skips incumbent cell revisits within a sweep", {
  degree_search <- getFromNamespace(".np_degree_search", "np")

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
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.degree.search.warn.grid = 3L)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(20),
    x2 = runif(20)
  )
  dat$y <- dat$x1 + dat$x2 + rnorm(nrow(dat), sd = 0.05)

  expect_warning(
    np::npregbw(
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
  old_opts <- options(
    np.messages = FALSE,
    np.tree = FALSE,
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
    np::npregbw(
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
  old_opts <- options(np.messages = TRUE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(20)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  msgs <- with_np_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 20, by = 0.5))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("np"), inherits = FALSE)(
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
  expect_true(any(grepl("best (", msgs, fixed = TRUE)))

  coord_msgs <- with_np_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 40, by = 0.5))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("np"), inherits = FALSE)(
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
