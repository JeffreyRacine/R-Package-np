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

  user_starts <- build_degree_starts(
    initial = c(1L, 1L, 1L, 1L),
    lower = rep(0L, 4L),
    upper = rep(10L, 4L),
    basis = "tensor",
    nobs = 100L,
    nmulti = 3L,
    random.seed = 42L,
    user_supplied = TRUE
  )

  expect_identical(starts5[seq_len(3L), , drop = FALSE], starts3)
  expect_identical(as.integer(user_starts[1L, ]), c(1L, 1L, 1L, 1L))
  expect_true(all(t(t(starts5) >= rep(0L, 4L))))
  expect_true(all(t(t(starts5) <= rep(10L, 4L))))
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

test_that("NOMAD Powell hot-start option helper controls internal remin explicitly", {
  hot_opt_args <- getFromNamespace(".np_nomad_powell_hotstart_opt_args", "npRmpi")

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
  attach_summary <- getFromNamespace(".np_attach_nomad_restart_summary", "npRmpi")

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

  expect_gt(length(msgs), 0L)
  expect_true(any(grepl("(Selecting degree and bandwidth|NOMAD degree/bw|Exhaustive degree/bw|Auto:NOMAD degree/bw|Auto:exhaustive degree/bw)", msgs)))
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

  expect_gt(length(coord_msgs), 0L)
})

test_that("cell-search baselines belong to the requested degree lattice", {
  degree_search <- getFromNamespace(".np_degree_search", "npRmpi")
  withr::local_options(list(np.messages = FALSE))
  grids <- list(list(1:2), list(0:1, 2:3), list(c(1L, 3L)))
  for (candidates in grids) {
    for (method in c("exhaustive", "coordinate")) {
      for (direction in c("min", "max")) {
        seen <- new.env(parent = emptyenv())
        seen$degrees <- list()
        baseline <- integer(length(candidates))
        result <- degree_search(
          method = method, candidates = candidates,
          baseline_degree = baseline,
          start_degree = vapply(candidates, max, 1L),
          restarts = 1L, verify = TRUE,
          direction = direction,
          eval_fun = function(degree) {
            seen$degrees[[length(seen$degrees) + 1L]] <- degree
            value <- sum(degree)
            list(objective = if (direction == "min") value else -value,
                 payload = list(degree = degree), num.feval = 1L)
          }
        )
        legal <- function(degree) all(vapply(seq_along(candidates), function(j)
          degree[j] %in% candidates[[j]], TRUE))
        expect_true(all(vapply(seen$degrees, legal, TRUE)))
        expect_identical(result$baseline$degree, vapply(candidates, min, 1L))
        expect_identical(result$best$degree, vapply(candidates, min, 1L))
        expect_true(result$completed)
        expect_true(result$certified)
      }
    }
  }
})

test_that("a positive-bound search cannot recover an excluded baseline", {
  degree_search <- getFromNamespace(".np_degree_search", "npRmpi")
  withr::local_options(list(np.messages = FALSE))
  evaluate <- function(degree) {
    if (degree > 0L)
      stop(structure(list(message = "no admissible model", call = NULL),
                     class = c("np_degree_candidate_invalid", "error", "condition")))
    list(objective = 0, payload = list(degree = degree))
  }
  expect_error(degree_search(method = "exhaustive", candidates = list(1:2),
    baseline_degree = 0L, start_degree = 1L, eval_fun = evaluate),
    "no admissible model")

  seen <- new.env(parent = emptyenv())
  seen$count <- 0L
  result <- degree_search(method = "exhaustive", candidates = list(1:2),
    baseline_degree = 0L, start_degree = 1L, eval_fun = function(degree) {
      seen$count <- seen$count + 1L
      if (seen$count > 1L)
        stop(structure(list(message = "interrupt"), class = c("interrupt", "condition")))
      list(objective = degree, payload = list(degree = degree))
    })
  expect_true(result$interrupted)
  expect_identical(result$best$degree, 1L)
  expect_identical(result$baseline$degree, 1L)
})

test_that("degree restarts fill finite grids and refuse impossible counts", {
  f <- getFromNamespace(".np_degree_restart_starts", "npRmpi")
  expect_identical(f(list(0:1, 0:1), 3L, list(c(0L, 0L))),
                   list(c(1L, 1L), c(1L, 0L), c(0L, 1L)))
  # Freeze the deterministic completion order when several points are missing.
  expect_identical(f(list(0:2, 0:2), 9L),
    list(c(0L,0L),c(2L,2L),c(1L,1L),c(1L,2L),c(2L,1L),
         c(0L,1L),c(0L,2L),c(1L,0L),c(2L,0L)))
  # Two exclusions, including a modular seed, on the default degree lattice.
  expect_identical(f(list(0:3,0:3), 6L, list(c(0L,0L),c(1L,1L))),
    list(c(3L,3L),c(1L,2L),c(2L,0L),c(3L,2L),c(0L,1L),c(0L,2L)))
  expect_error(f(list(0:3,0:3), 15L, list(c(0L,0L),c(1L,1L))),
               "exceeds the 14 available distinct")
  expect_error(f(list(0:1), 2L, list(0L)), "exceeds the 1 available distinct")
  expect_error(f(list(0L), 1L, list(0L)), "exceeds the 0 available distinct")
  expect_identical(f(list(0:1), 1L, list(0L, 0L, -1L, c(0L, 1L))), list(1L))
  expect_identical(f(list(0:1), 2L, list(-1L)), list(0L, 1L))
  expect_identical(f(list(4L, 2:3, 9L), 2L), list(c(4L, 2L, 9L), c(4L, 3L, 9L)))
  # Coordinate labels are metadata and must not alter modular-cycle detection.
  for (g in list(list(0:1,0:1), list(0:3,0:3), list(4L,2:3,9L))) {
    named <- setNames(g, paste0("x", seq_along(g)))
    available <- prod(lengths(g))-1L
    exclude <- list(vapply(g, min, 1L))
    expect_identical(f(named, available, exclude), f(g, available, exclude))
    shaped <- array(g, dim = c(length(g), 1L),
                    dimnames = list(names(named), "grid"))
    expect_identical(f(shaped, available, exclude), f(g, available, exclude))
  }
  set.seed(286); before <- .Random.seed
  expect_identical(f(list(0:1), 0L, list(0:1)), list())
  expect_error(f(list(0:1), 3L), "degree.restarts")
  expect_identical(.Random.seed, before)
})

test_that("restart completion preserves legacy prefixes and covers independent grids", {
  f <- getFromNamespace(".np_degree_restart_starts", "npRmpi")
  # An independent finite sequence oracle for these small cardinalities.
  reference <- function(g, exclude) {
    proposed <- list(vapply(g, min, 1L), vapply(g, max, 1L),
                     vapply(g, function(x) x[ceiling(length(x)/2)], 1L))
    for (k in seq_len(prod(lengths(g))))
      proposed[[length(proposed)+1L]] <- vapply(seq_along(g), function(j)
        as.integer(g[[j]][1L + (k*j) %% length(g[[j]])]), 1L)
    out <- list()
    for (z in proposed)
      if (!any(vapply(c(exclude, out), identical, TRUE, y=z))) out[[length(out)+1L]] <- z
    out
  }
  grids <- list(list(0:1, 0:1), list(0:2, 0:2), list(0:1, 1:3),
                list(2L, 0:2, 1:2), list(c(1L,3L), 0:2, 2:4))
  set.seed(738); before <- .Random.seed
  for (g in grids) {
    exclude <- list(vapply(g, min, 1L))
    old <- reference(g, exclude)
    for (k in seq_along(old)) expect_identical(f(g, k, exclude), old[seq_len(k)])
    grid <- as.matrix(expand.grid(g)); available <- nrow(grid)-1L
    starts <- f(g, available, exclude)
    expect_length(starts, available)
    expect_equal(length(unique(vapply(starts, paste, "", collapse=","))), available)
    expect_true(all(vapply(starts, function(z) any(apply(grid,1L,function(v)identical(as.integer(v),z))), TRUE)))
    expect_false(any(vapply(starts, identical, TRUE, y=exclude[[1L]])))
    expect_error(f(g, available+1L, exclude), "available distinct")
  }
  expect_identical(.Random.seed, before)
})

test_that("exhaustive degree searches ignore unused coordinate restarts", {
  f <- getFromNamespace(".np_degree_search", "npRmpi")
  called <- 0L
  objective <- function(degree) {called <<- called+1L; list(objective=(degree-1)^2, payload=degree)}
  z <- f(method="exhaustive", candidates=list(0:1), baseline_degree=0L,
         start_degree=0L, restarts=100L, eval_fun=objective)
  expect_identical(z$best$degree, 1L)
  expect_true(z$completed)
  expect_identical(z$restart.starts, list())
  expect_equal(called, 2L)
  called <- 0L
  expect_error(f(method="coordinate", candidates=list(0:1), baseline_degree=0L,
                 start_degree=0L, restarts=2L, eval_fun=objective), "available distinct")
  expect_identical(called, 0L)
})
