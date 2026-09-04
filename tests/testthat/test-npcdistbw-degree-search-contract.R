library(npRmpi)

skip_slow_npcdistbw_degree_search <- function() {
  skip_on_cran()
}

with_nprmpi_npcdist_degree_bindings <- function(bindings, code) {
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

test_that("npcdistbw exhaustive degree search matches manual profile minimum", {
  skip_slow_npcdistbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(24)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  bw0 <- npcdistbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 0L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  bw1 <- npcdistbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  auto <- npcdistbw(
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

  expect_s3_class(auto, "condbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))

  manual <- npcdistbw(
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

test_that("npcdistbw coordinate search can be exhaustively certified on a small grid", {
  skip_slow_npcdistbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(22),
    x2 = runif(22)
  )
  dat$y <- dat$x1 + dat$x2^2 + rnorm(nrow(dat), sd = 0.08)

  exhaustive <- npcdistbw(
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
  coordinate <- npcdistbw(
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
})

test_that("npcdistbw automatic degree search enforces pilot guardrails", {
  skip_slow_npcdistbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(y = rnorm(20), x = runif(20))

  expect_error(
    npcdistbw(
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

  bw <- npcdistbw(
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
  )

  expect_s3_class(bw, "condbandwidth")
  expect_false(isTRUE(bw$bernstein.basis))
  expect_lte(max(as.integer(bw$degree)), 4L)
})

test_that("npcdist forwards automatic LP degree search through npcdistbw", {
  skip_slow_npcdistbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(20))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  fit <- npcdist(
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

  expect_false(is.null(fit$bws))
  expect_s3_class(fit$bws, "condbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})

test_that("npcdistbw NOMAD degree search backend improves over the baseline", {
  skip_slow_npcdistbw_degree_search()
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(18)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  bw <- npcdistbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L,
    ngrid = 30L
  )

  expect_s3_class(bw, "condbandwidth")
  expect_identical(bw$degree.search$mode, "nomad")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_gte(bw$degree.search$n.unique, 1L)
  expect_lte(bw$degree.search$best.fval, bw$degree.search$baseline.fval + 1e-10)
})

test_that("npcdistbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(18)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  bw <- npcdistbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))
  expect_equal(as.double(bw$total.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
  expect_equal(as.double(bw$degree.search$optim.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
})

test_that("npcdistbw native NOMAD route has an explicit crs availability guard", {
  require_body <- paste(
    deparse(body(get(".npcdistbw_nomad_native_require_crs", envir = asNamespace("npRmpi"), inherits = FALSE))),
    collapse = "\n"
  )

  expect_match(require_body, "requireNamespace\\(\"crs\",\\s*quietly = TRUE\\)")
  expect_true(grepl("packageVersion(\"crs\") < \"0.15.46\"", require_body, fixed = TRUE))
  expect_true(grepl("native npcdist NOMAD route requires crs >= 0.15-46", require_body, fixed = TRUE))
})

test_that("npcdistbw native prepared search result has a unique exact schema", {
  native_result <- list(
    best_point = c(0.4, 1),
    best_degree = 1L,
    first_degree = 0L,
    first_objective = 3,
    message = "ok",
    objective = 2,
    blackbox_evaluations = 4L,
    iterations = 2L,
    solution = c(0.4, 1),
    total_num.feval = 5,
    total_num.feval.fast = 3,
    compiled_callback_calls = 4L,
    best_num.feval = 5,
    best_num.feval.fast = 3,
    official_objective = 2,
    compiled_callback_failures = 0L,
    crs_callback_evaluations = 4L,
    cache_hits = 0L,
    cache_size = 4L,
    total_evaluations = 4L
  )

  with_nprmpi_npcdist_degree_bindings(list(
    .npcdistbw_nomad_degree_native_target = function(...) TRUE,
    .npcdistbw_nomad_native_require_crs = function(...) invisible(TRUE),
    .np_nomad_bw_restart_start_bounds = function(...) {
      list(lower = 0.1, upper = 1)
    },
    .np_nomad_prepare_solver_opts = function(...) list(),
    .npcdistbw_nomad_native_option_vectors = function(...) {
      list(names = character(), values = character())
    },
    .np_nomad_build_starts = function(...) matrix(c(0.4, 1), nrow = 1L),
    .npcdistbw_nomad_native_prepare_args = function(...) list(),
    .np_nomad_native_progress_begin = function(...) new.env(parent = emptyenv()),
    .np_nomad_native_progress_restart = function(...) invisible(NULL),
    .np_nomad_native_progress_end = function(...) invisible(NULL),
    .np_nomad_native_progress_abort = function(...) invisible(NULL),
    npNomadNativeSearchConditionalDistribution = function(...) native_result,
    .npcdistbw_nomad_eval_point = function(...) {
      list(objective = 2, num.feval = 0, num.feval.fast = 0)
    },
    .np_nomad_native_status = function(...) invisible(NULL)
  ), {
    result <- npRmpi:::npRmpiPreparedSearchConditionalDistribution(
      xdat = data.frame(x = seq(0.1, 0.9, length.out = 6L)),
      ydat = data.frame(y = seq(0.2, 0.8, length.out = 6L)),
      template = list(type = "fixed", scaling = TRUE, ybw = 0.4, xbw = numeric()),
      setup = list(
        type = "fixed",
        cont_flat = 1L,
        cont_scale = 1,
        cat_flat = integer(),
        cat_upper = numeric(),
        nobs = 6L
      ),
      reg.args = list(),
      opt.args = list(nomad.opts = list(), nomad.remin = FALSE),
      degree.search = list(
        engine = "nomad",
        start.degree = 0L,
        candidates = list(0:1),
        lower = 0L,
        upper = 1L,
        basis = "glp",
        nobs = 6L,
        start.user = FALSE,
        bernstein.basis = TRUE
      ),
      x0 = c(0.4, 0),
      bbin = c(0L, 1L),
      lb = c(0.1, 0),
      ub = c(1, 1),
      source = "unit",
      reason = "schema"
    )

    expected_names <- c(
      "method", "source", "reason", "direction", "verify", "completed",
      "certified", "interrupted", "baseline", "best", "best_payload",
      "best_point", "n.unique", "n.visits", "n.cached", "nomad.time",
      "powell.time", "optim.time", "grid.size", "best.restart",
      "nomad.remin", "nomad.remin.index", "nomad.remin.roundtrip",
      "restart.starts", "restart.degree.starts", "restart.bandwidth.starts",
      "restart.start.info", "restart.results", "trace",
      "native.diagnostics", "num.feval.total", "num.feval.fast.total"
    )
    expect_identical(anyDuplicated(names(result)), 0L)
    expect_identical(names(result), expected_names)
    expect_identical(result$method, "nomad")
    expect_identical(result$source, "unit")
    expect_identical(result$reason, "schema")
    expect_identical(result$direction, "min")
    expect_identical(result$verify, FALSE)
    expect_identical(result$completed, TRUE)
    expect_identical(result$certified, FALSE)
    expect_identical(result$interrupted, FALSE)
    expect_type(result$baseline, "list")
    expect_type(result$best, "list")
    expect_null(result$best_payload)
    expect_identical(result$best_point, c(0.4, 1))
    expect_identical(result$n.unique, 4L)
    expect_identical(result$n.visits, 4L)
    expect_identical(result$n.cached, 0L)
    expect_type(result$nomad.time, "double")
    expect_identical(result$powell.time, NA_real_)
    expect_type(result$optim.time, "double")
    expect_identical(result$grid.size, NA_integer_)
    expect_identical(result$best.restart, 1L)
    expect_identical(result$nomad.remin, FALSE)
    expect_identical(result$nomad.remin.index, NA_integer_)
    expect_null(result$nomad.remin.roundtrip)
    expect_type(result$restart.starts, "list")
    expect_type(result$restart.degree.starts, "list")
    expect_type(result$restart.bandwidth.starts, "list")
    expect_type(result$restart.start.info, "list")
    expect_type(result$restart.results, "list")
    expect_s3_class(result$trace, "data.frame")
    expect_type(result$native.diagnostics, "list")
    expect_identical(result$num.feval.total, 5)
    expect_identical(result$num.feval.fast.total, 3)

    metadata <- npRmpi:::.np_degree_search_metadata(result, "min")
    expect_identical(metadata$mode, "nomad")
    expect_identical(metadata$source, "unit")
    expect_identical(metadata$reason, "schema")
    expect_identical(metadata$best.degree, 1L)
    expect_identical(metadata$best.fval, 2)
    expect_identical(metadata$best.restart, 1L)

    attached <- npRmpi:::.npcdistbw_attach_degree_search(
      bws = structure(list(), class = "condbandwidth"),
      search_result = result
    )
    expect_s3_class(attached, "condbandwidth")
    expect_identical(attached$degree.search, metadata)
    expect_identical(unserialize(serialize(result, NULL)), result)
    expect_identical(unserialize(serialize(attached, NULL)), attached)
  })
})
