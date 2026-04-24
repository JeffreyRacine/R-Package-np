npreg_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

npreg_fit_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

npreg_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

normalize_npreg_progress_output <- function(output) {
  text <- paste(output, collapse = "\n")
  text <- gsub("\r", "\n", text, fixed = TRUE)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  lines <- sub("[[:space:]]+$", "", lines)
  lines[nzchar(lines)]
}

extract_npreg_case_block <- function(lines, case = "npreg") {
  start <- grep(sprintf("^CASE_START %s$", case), lines)
  done <- grep(sprintf("^CASE_DONE %s$", case), lines)

  if (!length(start) || !length(done)) {
    return(character())
  }

  lines[seq.int(start[[1L]], done[[1L]])]
}

npreg_case_progress_lines <- function(lines, case = "npreg") {
  block <- extract_npreg_case_block(lines, case)
  block[grepl("^\\[npRmpi\\] ", block)]
}

expect_npreg_powell_progress_surface <- function(lines) {
  powell.lines <- lines[grepl("^\\[npRmpi\\] Refining bandwidth \\(", lines)]
  info <- paste(lines, collapse = "\n")

  expect_true(length(powell.lines) > 0L, info = info)
  expect_false(any(grepl(
    "^\\[npRmpi\\] Bandwidth selection \\(Refining NOMAD solution",
    lines
  )), info = info)
  expect_false(any(grepl("best (", powell.lines, fixed = TRUE)), info = info)
  expect_true(any(grepl(
    "^\\[npRmpi\\] Refining bandwidth \\(elapsed [0-9]+\\.[0-9]s, degree \\([^)]*\\), iter [0-9]+\\)$",
    powell.lines
  )), info = info)
}

ensure_npreg_subprocess_nprmpi_lib <- local({
  lib.path.cache <- NULL

  function() {
    if (!is.null(lib.path.cache) && dir.exists(lib.path.cache)) {
      return(lib.path.cache)
    }

    pkg.root <- tryCatch(
      normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
      error = function(e) ""
    )
    if (!nzchar(pkg.root)) {
      return(NULL)
    }

    lib.path.cache <<- tempfile("npRmpi-npreg-subprocess-lib-")
    dir.create(lib.path.cache, recursive = TRUE, showWarnings = FALSE)

    cmd <- file.path(R.home("bin"), "R")
    out <- suppressWarnings(system2(
      cmd,
      c("CMD", "INSTALL", "--no-test-load", "-l", lib.path.cache, pkg.root),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(out, "status")
    if (is.null(status)) {
      status <- 0L
    }

    if (status != 0L) {
      warning(paste(out, collapse = "\n"))
      unlink(lib.path.cache, recursive = TRUE, force = TRUE)
      lib.path.cache <<- NULL
      return(NULL)
    }

    lib.path.cache
  }
})

npreg_local_env <- function() {
  skip_on_cran()
  skip_if_not(
    tolower(Sys.getenv("NOT_CRAN", "")) %in% c("true", "1", "yes"),
    "extended local npreg progress contract"
  )
  lib.path <- ensure_npreg_subprocess_nprmpi_lib()
  skip_if(is.null(lib.path), "local npRmpi install unavailable for subprocess smoke")

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    "NP_RMPI_NO_REUSE_SLAVES=1"
  )
}

npreg_subprocess_lib <- function(env) {
  r_libs <- env[grepl("^R_LIBS=", env)]
  if (!length(r_libs)) {
    return(NULL)
  }

  lib_paths <- strsplit(sub("^R_LIBS=", "", r_libs[[1L]]), .Platform$path.sep, fixed = TRUE)[[1L]]
  lib_paths[[1L]]
}

run_npreg_cmd_subprocess <- function(cmd, args = character(), timeout = 60L, env = character()) {
  out <- suppressWarnings(system2(
    cmd,
    args,
    stdout = TRUE,
    stderr = TRUE,
    timeout = timeout,
    env = env
  ))
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list(status = as.integer(status), output = out)
}

npreg_is_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

run_npreg_session_progress_contract <- function(timeout = 120L) {
  env <- npreg_local_env()

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(",
      "  np.messages = TRUE,",
      "  np.tree = FALSE,",
      "  np.progress.bandwidth.enhanced = TRUE,",
      "  np.progress.start.grace.known.sec = 0,",
      "  np.progress.start.grace.unknown.sec = 0,",
      "  np.progress.interval.known.sec = 0,",
      "  np.progress.interval.unknown.sec = 0,",
      "  width = 80",
      ")",
      "Sys.setenv(COLUMNS = '80')",
      "assignInNamespace('.np_progress_is_interactive', function() TRUE, ns = 'npRmpi')",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(20260404)",
      "n <- 24L",
      "dat <- data.frame(x = seq(-0.9, 0.9, length.out = n))",
      "dat$y <- sin(2 * pi * dat$x) + 0.35 * dat$x",
      "cat('CASE_START npreg\\n')",
      "fit <- npreg(y ~ x, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'npregression'))",
      "cat('CASE_DONE npreg\\n')"
    ),
    timeout = timeout,
    env = env
  )

  list(
    status = res$status,
    lines = normalize_npreg_progress_output(res$output),
    raw = res$output
  )
}

run_npreg_attach_progress_contract <- function(timeout = 120L) {
  env_common <- npreg_local_env()
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-npreg-attach-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(",
    "  np.messages = TRUE,",
    "  np.tree = FALSE,",
    "  np.progress.bandwidth.enhanced = TRUE,",
    "  np.progress.start.grace.known.sec = 0,",
    "  np.progress.start.grace.unknown.sec = 0,",
    "  np.progress.interval.known.sec = 0,",
    "  np.progress.interval.unknown.sec = 0,",
    "  width = 140",
    ")",
    "Sys.setenv(COLUMNS = '140')",
    "assignInNamespace('.np_progress_is_interactive', function() TRUE, ns = 'npRmpi')",
    "is.master <- isTRUE(npRmpi.init(mode = 'attach', quiet = TRUE, autodispatch = TRUE))",
    "if (is.master) {",
    "  on.exit({",
    "    try(npRmpi.quit(mode = 'attach'), silent = TRUE)",
    "    try(mpi.quit(), silent = TRUE)",
    "  }, add = TRUE)",
    "  set.seed(20260404)",
    "  n <- 24L",
    "  dat <- data.frame(x = seq(-0.9, 0.9, length.out = n))",
    "  dat$y <- sin(2 * pi * dat$x) + 0.35 * dat$x",
    "  cat('CASE_START npreg\\n')",
    "  fit <- npreg(y ~ x, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L)",
    "  stopifnot(inherits(fit, 'npregression'))",
    "  cat('CASE_DONE npreg\\n')",
    "}"
  ), script, useBytes = TRUE)

  res <- run_npreg_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = timeout,
    env = c(
      env_common,
      "R_PROFILE_USER=",
      "R_PROFILE=",
      "FI_TCP_IFACE=en0",
      "FI_PROVIDER=tcp",
      "FI_SOCKETS_IFACE=en0"
    )
  )
  if (res$status != 0L) {
    res <- run_npreg_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = timeout,
      env = c(
        env_common,
        "R_PROFILE_USER=",
        "R_PROFILE=",
        "FI_TCP_IFACE=lo0",
        "FI_PROVIDER=tcp",
        "FI_SOCKETS_IFACE=lo0"
      )
    )
  }

  if (res$status != 0L && npreg_is_mpi_init_env_failure(res$output)) {
    skip("MPI runtime interface unavailable in this environment for attach-mode smoke")
  }

  list(
    status = res$status,
    lines = normalize_npreg_progress_output(res$output),
    raw = res$output
  )
}

run_npreg_profile_progress_contract <- function(timeout = 120L) {
  env_common <- npreg_local_env()
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  lib.path <- npreg_subprocess_lib(env_common)
  skip_if(is.null(lib.path), "local npRmpi install unavailable for subprocess smoke")
  profile.path <- file.path(lib.path, "npRmpi", "Rprofile")
  skip_if(!file.exists(profile.path), "npRmpi profile template unavailable in subprocess lib")

  script <- tempfile("npRmpi-npreg-profile-", fileext = ".R")
  batch_file <- tempfile("npRmpi-npreg-profile-", fileext = ".Rout")
  on.exit(unlink(c(script, batch_file)), add = TRUE)
  writeLines(c(
    "mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE)",
    "mpi.bcast.cmd(options(",
    "  np.messages = TRUE,",
    "  np.tree = FALSE,",
    "  np.progress.bandwidth.enhanced = TRUE,",
    "  np.progress.start.grace.known.sec = 0,",
    "  np.progress.start.grace.unknown.sec = 0,",
    "  np.progress.interval.known.sec = 0,",
    "  np.progress.interval.unknown.sec = 0,",
    "  width = 140",
    "), caller.execute = TRUE)",
    "mpi.bcast.cmd(Sys.setenv(COLUMNS = '140'), caller.execute = TRUE)",
    "mpi.bcast.cmd(assignInNamespace('.np_progress_is_interactive', function() TRUE, ns = 'npRmpi'), caller.execute = TRUE)",
    "set.seed(20260404)",
    "n <- 24L",
    "dat <- data.frame(x = seq(-0.9, 0.9, length.out = n))",
    "dat$y <- sin(2 * pi * dat$x) + 0.35 * dat$x",
    "mpi.bcast.Robj2slave(dat)",
    "cat('CASE_START npreg\\n')",
    "fit <- mpi.bcast.cmd(npreg(y ~ x, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L), caller.execute = TRUE)",
    "stopifnot(inherits(fit, 'npregression'))",
    "cat('CASE_DONE npreg\\n')",
    "mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)"
  ), script, useBytes = TRUE)

  res <- run_npreg_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "R"), "CMD", "BATCH", "--no-save", script, batch_file),
    timeout = timeout,
    env = c(
      env_common,
      sprintf("R_PROFILE_USER=%s", profile.path),
      "R_PROFILE=",
      "NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=120",
      "FI_TCP_IFACE=en0",
      "FI_PROVIDER=tcp",
      "FI_SOCKETS_IFACE=en0"
    )
  )
  if (res$status != 0L) {
    res <- run_npreg_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "R"), "CMD", "BATCH", "--no-save", script, batch_file),
      timeout = timeout,
      env = c(
        env_common,
        sprintf("R_PROFILE_USER=%s", profile.path),
        "R_PROFILE=",
        "NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=120",
        "FI_TCP_IFACE=lo0",
        "FI_PROVIDER=tcp",
        "FI_SOCKETS_IFACE=lo0"
      )
    )
  }

  if (res$status != 0L && npreg_is_mpi_init_env_failure(res$output)) {
    skip("MPI runtime interface unavailable in this environment for profile-mode smoke")
  }

  batch_lines <- if (file.exists(batch_file)) readLines(batch_file, warn = FALSE) else character()

  list(
    status = res$status,
    lines = normalize_npreg_progress_output(batch_lines),
    raw = c(batch_lines, res$output)
  )
}

capture_npreg_fit_progress_trace <- function(expr,
                                             force_renderer = "single_line",
                                             now = function() 0,
                                             interactive = TRUE,
                                             master = TRUE) {
  expr_env <- parent.frame()
  expr <- substitute(expr)
  trace <- list()

  recorder <- function(snapshot, event = c("render", "finish", "abort")) {
    event <- match.arg(event)
    if (isTRUE(getFromNamespace(".np_progress_is_message_muffled", "npRmpi")())) {
      return(invisible(snapshot))
    }
    trace[[length(trace) + 1L]] <<- list(
      event = event,
      renderer = force_renderer,
      id = snapshot$id,
      kind = snapshot$kind,
      current = snapshot$current,
      total = snapshot$total,
      detail = snapshot$detail,
      line = snapshot$line,
      started_at = snapshot$started_at,
      now = snapshot$now,
      last_width = snapshot$last_width
    )
    invisible(snapshot)
  }

  value <- with_nprmpi_progress_bindings(
    c(
      if (!is.null(force_renderer)) {
        list(.np_progress_renderer_for_surface = function(surface, capability) force_renderer)
      } else {
        list()
      },
      list(
        .np_progress_render_legacy = recorder,
        .np_progress_render_single_line = recorder,
        .np_progress_is_interactive = function() interactive,
        .np_progress_is_master = function() master,
        .np_progress_now = now,
        .np_progress_output_width = function() 500L
      )
    ),
    {
      reset <- getFromNamespace(".np_progress_reset_registry", "npRmpi")
      reset()
      on.exit(reset(), add = TRUE)
      eval(expr, envir = expr_env)
    }
  )

  list(
    value = value,
    trace = trace,
    final_line = if (length(trace)) trace[[length(trace)]]$line else NULL
  )
}

make_npreg_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 24L
  dat <- data.frame(x = seq(-0.9, 0.9, length.out = n))
  dat$y <- sin(2 * pi * dat$x) + 0.35 * dat$x

  list(
    dat = dat,
    tx = dat["x"],
    y = dat$y,
    bw = npregbw(
      xdat = dat["x"],
      ydat = dat$y,
      bws = 0.22,
      bandwidth.compute = FALSE
    )
  )
}

test_that("npreg direct bws fit emits known-total fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npreg_fit_progress_trace(
    npreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)

  expect_s3_class(actual$value, "npregression")
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/24 \\([0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})

test_that("npreg direct bws fit stays silent below start grace without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0.75,
    np.progress.interval.known.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npreg_fit_progress_trace(
    npreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    now = npreg_fit_progress_time_values(c(0, 0.2, 0.4, 0.6))
  )

  expect_length(actual$trace, 0L)
})

test_that("npreg bw to fit route hands off immediately into single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npreg_fit_progress_trace(
    npreg(
      y ~ x,
      data = fixture$dat,
      nmulti = 1L
    ),
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )
  fit.finish.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )

  expect_s3_class(actual$value, "npregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npreg nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npreg_fit_progress_trace(
    npreg(
      y ~ x,
      data = fixture$dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Selecting degree and bandwidth", lines)
  powell.pos <- grep("^\\[npRmpi\\] Refining bandwidth \\(", lines)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )
  fit.finish.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )

  expect_s3_class(actual$value, "npregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_npreg_powell_progress_surface(lines)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("predict.npregression re-entry emits fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npreg_fit_progress_fixture()

  fit <- npreg(
    bws = fixture$bw,
    txdat = fixture$tx,
    tydat = fixture$y
  )

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npreg_fit_progress_trace(
    predict(fit, newdata = fixture$dat[c(2L, 7L), "x", drop = FALSE]),
    now = npreg_fit_progress_time_counter()
  )

  lines <- npreg_fit_progress_lines(actual)

  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})

test_that("session npreg nomad route keeps visible powell handoff in subprocess", {
  skip_if_not_installed("crs")

  actual <- run_npreg_session_progress_contract()
  block <- npreg_case_progress_lines(actual$lines)
  bandwidth.pos <- grep("^\\[npRmpi\\] Selecting degree and bandwidth", block)
  powell.pos <- grep("^\\[npRmpi\\] Refining bandwidth \\(", block)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    block
  )

  expect_identical(actual$status, 0L, info = paste(actual$raw, collapse = "\n"))
  expect_true(length(bandwidth.pos) > 0L, info = paste(block, collapse = "\n"))
  expect_true(length(powell.pos) > 0L, info = paste(block, collapse = "\n"))
  expect_identical(length(fit.start.pos), 1L, info = paste(block, collapse = "\n"))
  expect_npreg_powell_progress_surface(block)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
})

test_that("attach npreg nomad route does not duplicate fit-start lines in subprocess", {
  skip_if_not_installed("crs")
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")

  actual <- run_npreg_attach_progress_contract()
  block <- npreg_case_progress_lines(actual$lines)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    block
  )

  expect_identical(actual$status, 0L, info = paste(actual$raw, collapse = "\n"))
  expect_npreg_powell_progress_surface(block)
  expect_identical(length(fit.start.pos), 1L, info = paste(block, collapse = "\n"))
})

test_that("profile npreg nomad route keeps visible powell handoff in subprocess", {
  skip_if_not_installed("crs")
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_PROFILE_TEST"), "1"),
              "set NP_RMPI_ENABLE_PROFILE_TEST=1 to run profile-mode smoke")

  actual <- run_npreg_profile_progress_contract()
  block <- npreg_case_progress_lines(actual$lines)
  bandwidth.pos <- grep("^\\[npRmpi\\] Selecting degree and bandwidth", block)
  powell.pos <- grep("^\\[npRmpi\\] Refining bandwidth \\(", block)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    block
  )

  expect_identical(actual$status, 0L, info = paste(actual$raw, collapse = "\n"))
  expect_true(length(bandwidth.pos) > 0L, info = paste(block, collapse = "\n"))
  expect_true(length(powell.pos) > 0L, info = paste(block, collapse = "\n"))
  expect_identical(length(fit.start.pos), 1L, info = paste(block, collapse = "\n"))
  expect_npreg_powell_progress_surface(block)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
})
