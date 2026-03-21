with_npRmpi_degree_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("npRmpi")
  old <- lapply(names(bindings), function(name) get(name, envir = ns, inherits = FALSE))
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked)
      unlockBinding(name, ns)
    assign(name, bindings[[name]], envir = ns)
    if (was_locked)
      lockBinding(name, ns)
  }

  on.exit({
    for (name in names(old)) {
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked)
        unlockBinding(name, ns)
      assign(name, old[[name]], envir = ns)
      if (was_locked)
        lockBinding(name, ns)
    }
  }, add = TRUE)

  eval(code, envir = parent.frame())
}

capture_degree_messages_only <- function(expr) {
  messages <- character()
  withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  messages
}

degree_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

test_that("npregbw exhaustive degree search matches manual profile minimum", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(28)))
  dat$y <- dat$x^2 + rnorm(nrow(dat), sd = 0.05)

  bw0 <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 0L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  bw1 <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  auto <- npregbw(
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

  expect_s3_class(auto, "rbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))
  expect_identical(nrow(auto$degree.search$trace), auto$degree.search$n.unique)
  expect_identical(auto$degree.search$n.cached, auto$degree.search$n.visits - auto$degree.search$n.unique)

  manual <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  expect_null(manual$degree.search)
})

test_that("npregbw coordinate search can be exhaustively certified on a small grid", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(26),
    x2 = runif(26)
  )
  dat$y <- dat$x1 + dat$x2^2 + rnorm(nrow(dat), sd = 0.05)

  exhaustive <- npregbw(
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
  )
  coordinate <- npregbw(
    y ~ x1 + x2,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = TRUE,
    degree.restarts = 1L,
    degree.max.cycles = 4L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(coordinate$degree.search$mode, "coordinate")
  expect_true(isTRUE(coordinate$degree.search$completed))
  expect_true(isTRUE(coordinate$degree.search$certified))
  expect_equal(as.integer(coordinate$degree), as.integer(exhaustive$degree))
  expect_equal(coordinate$fval, exhaustive$fval, tolerance = 1e-10)
  expect_lte(coordinate$degree.search$best.fval, coordinate$degree.search$baseline.fval + 1e-10)
  expect_identical(nrow(coordinate$degree.search$trace), coordinate$degree.search$n.unique)
  expect_identical(coordinate$degree.search$n.cached, coordinate$degree.search$n.visits - coordinate$degree.search$n.unique)
})

test_that("npregbw automatic degree search enforces pilot guardrails", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(y = rnorm(24), x = runif(24))

  expect_error(
    npregbw(
      y ~ x,
      data = dat,
      regtype = "lc",
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "automatic degree search currently requires regtype='lp'"
  )

  expect_error(
    npregbw(
      y ~ x,
      data = dat,
      regtype = "lp",
      bernstein.basis = FALSE,
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 4L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "degree.max <= 3"
  )

  bw <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    bernstein.basis = FALSE,
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_false(isTRUE(bw$bernstein.basis))
})

test_that("npreg forwards automatic LP degree search through npregbw", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(24))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  fit <- npreg(
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

  expect_s3_class(fit, "npregression")
  expect_s3_class(fit$bws, "rbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})

test_that("npregbw NOMAD degree search backend improves over the baseline", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(24)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 3L
  )

  expect_s3_class(bw, "rbandwidth")
  expect_identical(bw$degree.search$mode, "nomad")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_gte(bw$degree.search$n.unique, 1L)
  expect_lte(bw$degree.search$best.fval, bw$degree.search$baseline.fval + 1e-10)
  expect_identical(length(bw$degree.search$restart.degree.starts), 3L)
  expect_identical(length(bw$degree.search$restart.results), 3L)
  expect_identical(bw$degree.search$restart.results[[1L]]$restart, 1L)
  expect_identical(length(bw$nomad.restart.fval), 3L)
  expect_identical(length(bw$nomad.restart.results), 3L)
  expect_identical(length(bw$nomad.restart.degree.starts), 3L)
  restart.objectives <- unlist(lapply(
    bw$degree.search$restart.results,
    function(x) if (is.null(x$objective)) NA_real_ else as.numeric(x$objective[1L])
  ))
  restart.objectives <- restart.objectives[is.finite(restart.objectives)]
  expect_true(length(restart.objectives) >= 1L)
  expect_lte(bw$degree.search$best.fval, min(restart.objectives) + 1e-10)
  expect_identical(as.numeric(bw$nomad.restart.fval), as.numeric(unlist(lapply(
    bw$nomad.restart.results,
    function(x) if (is.null(x$objective)) NA_real_ else as.numeric(x$objective[1L])
  ))))
  expect_identical(bw$nomad.best.restart, which.min(bw$nomad.restart.fval))
})

test_that("npregbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(24)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))
})

test_that("npregbw NOMAD degree search fails fast when crs is unavailable", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(16), y = rnorm(16))

  expect_error(
    with_npRmpi_degree_bindings(
      list(.np_nomad_require_crs = function() stop("crs missing", call. = FALSE)),
      npregbw(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "nomad",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    ),
    "crs missing"
  )
})

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
  expect_true(any(grepl("Selecting polynomial degree and bw", msgs, fixed = TRUE)))
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

test_that("npregbw NOMAD plus Powell progress keeps stage text on managed lines", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    npRmpi.autodispatch = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260321)
  dat <- data.frame(x = sort(runif(24)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  msgs <- with_npRmpi_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 40, by = 0.25))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("npRmpi"), inherits = FALSE)(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "nomad+powell",
        degree.min = 0L,
        degree.max = 2L,
        degree.verify = FALSE,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L,
        max.bb.eval = 15
      )
    )
  )

  expect_false(any(grepl("^\\[np\\] Starting NOMAD automatic polynomial degree search at degree ", msgs)))
  expect_false(any(grepl("^\\[np\\] Refining NOMAD solution with one Powell hot start at degree ", msgs)))
  expect_true(any(grepl("starting at degree \\(", msgs)))
})
