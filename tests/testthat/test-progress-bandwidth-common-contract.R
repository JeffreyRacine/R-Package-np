progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

shadow_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

normalize_bandwidth_progress_output <- function(output) {
  text <- paste(output, collapse = "\n")
  text <- gsub("\r", "\n", text, fixed = TRUE)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  lines <- sub("[[:space:]]+$", "", lines)
  lines[nzchar(lines)]
}

extract_bandwidth_case_block <- function(lines, case) {
  start <- grep(sprintf("^CASE_START %s$", case), lines)
  done <- grep(sprintf("^CASE_DONE %s$", case), lines)

  if (!length(start) || !length(done)) {
    return(character())
  }

  lines[seq.int(start[[1L]], done[[1L]])]
}

bandwidth_case_progress_lines <- function(lines, case) {
  block <- extract_bandwidth_case_block(lines, case)
  block[grepl("^\\[npRmpi\\] ", block)]
}

bandwidth_adjacent_duplicates <- function(lines) {
  if (length(lines) < 2L) {
    return(character())
  }

  duplicated <- lines[-1L] == lines[-length(lines)]
  lines[-1L][duplicated]
}

ensure_bandwidth_subprocess_nprmpi_lib <- local({
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

    lib.path.cache <<- tempfile("npRmpi-bandwidth-subprocess-lib-")
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

bandwidth_local_env <- function() {
  skip_on_cran()
  skip_if_not(
    tolower(Sys.getenv("NOT_CRAN", "")) %in% c("true", "1", "yes"),
    "extended local bandwidth progress contract"
  )
  lib.path <- ensure_bandwidth_subprocess_nprmpi_lib()
  skip_if(is.null(lib.path), "local npRmpi install unavailable for subprocess smoke")

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    "NP_RMPI_NO_REUSE_SLAVES=1"
  )
}

bandwidth_subprocess_lib <- function(env) {
  r_libs <- env[grepl("^R_LIBS=", env)]
  if (!length(r_libs)) {
    return(NULL)
  }

  lib_paths <- strsplit(sub("^R_LIBS=", "", r_libs[[1L]]), .Platform$path.sep, fixed = TRUE)[[1L]]
  lib_paths[[1L]]
}

run_bandwidth_cmd_subprocess <- function(cmd, args = character(), timeout = 60L, env = character()) {
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

bandwidth_is_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

run_bandwidth_session_contract <- function(timeout = 120L) {
  env <- bandwidth_local_env()

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
      "n <- 12L",
      "x <- sort(runif(n))",
      "y <- x + 0.35 * x^2 + rnorm(n, sd = 0.12)",
      "dat <- data.frame(y = y, x = x)",
      "cat('CASE_START npregbw\\n')",
      "bw.reg <- npregbw(y ~ x, data = dat, regtype = 'lc', bwmethod = 'cv.ls', nmulti = 2L)",
      "cat('CASE_DONE npregbw\\n')",
      "cat('CASE_START npcdensbw\\n')",
      "bw.cdens <- npcdensbw(y ~ x, data = dat, bwmethod = 'cv.ls', nmulti = 2L)",
      "cat('CASE_DONE npcdensbw\\n')"
    ),
    timeout = timeout,
    env = env
  )

  list(
    status = res$status,
    lines = normalize_bandwidth_progress_output(res$output),
    raw = res$output
  )
}

run_bandwidth_attach_contract <- function(timeout = 120L) {
  env_common <- bandwidth_local_env()
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-bandwidth-attach-", fileext = ".R")
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
    "  n <- 12L",
    "  x <- sort(runif(n))",
    "  y <- x + 0.35 * x^2 + rnorm(n, sd = 0.12)",
    "  dat <- data.frame(y = y, x = x)",
    "  cat('CASE_START npregbw\\n')",
    "  bw.reg <- npregbw(y ~ x, data = dat, regtype = 'lc', bwmethod = 'cv.ls', nmulti = 2L)",
    "  cat('CASE_DONE npregbw\\n')",
    "  cat('CASE_START npcdensbw\\n')",
    "  bw.cdens <- npcdensbw(y ~ x, data = dat, bwmethod = 'cv.ls', nmulti = 2L)",
    "  cat('CASE_DONE npcdensbw\\n')",
    "}"
  ), script, useBytes = TRUE)

  res <- run_bandwidth_cmd_subprocess(
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
    res <- run_bandwidth_cmd_subprocess(
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

  if (res$status != 0L && bandwidth_is_mpi_init_env_failure(res$output)) {
    skip("MPI runtime interface unavailable in this environment for attach-mode smoke")
  }

  list(
    status = res$status,
    lines = normalize_bandwidth_progress_output(res$output),
    raw = res$output
  )
}

run_bandwidth_profile_contract <- function(timeout = 120L) {
  env_common <- bandwidth_local_env()
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  lib.path <- bandwidth_subprocess_lib(env_common)
  skip_if(is.null(lib.path), "local npRmpi install unavailable for subprocess smoke")
  profile.path <- file.path(lib.path, "npRmpi", "Rprofile")
  skip_if(!file.exists(profile.path), "npRmpi profile template unavailable in subprocess lib")

  script <- tempfile("npRmpi-bandwidth-profile-", fileext = ".R")
  batch_file <- tempfile("npRmpi-bandwidth-profile-", fileext = ".Rout")
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
    "n <- 12L",
    "x <- sort(runif(n))",
    "y <- x + 0.35 * x^2 + rnorm(n, sd = 0.12)",
    "dat <- data.frame(y = y, x = x)",
    "mpi.bcast.Robj2slave(dat)",
    "cat('CASE_START npregbw\\n')",
    "bw.reg <- mpi.bcast.cmd(npregbw(y ~ x, data = dat, regtype = 'lc', bwmethod = 'cv.ls', nmulti = 2L), caller.execute = TRUE)",
    "cat('CASE_DONE npregbw\\n')",
    "cat('CASE_START npcdensbw\\n')",
    "bw.cdens <- mpi.bcast.cmd(npcdensbw(y ~ x, data = dat, bwmethod = 'cv.ls', nmulti = 2L), caller.execute = TRUE)",
    "cat('CASE_DONE npcdensbw\\n')",
    "mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)"
  ), script, useBytes = TRUE)

  res <- run_bandwidth_cmd_subprocess(
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
    res <- run_bandwidth_cmd_subprocess(
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

  if (res$status != 0L && bandwidth_is_mpi_init_env_failure(res$output)) {
    skip("MPI runtime interface unavailable in this environment for profile-mode smoke")
  }

  batch_lines <- if (file.exists(batch_file)) readLines(batch_file, warn = FALSE) else character()

  list(
    status = res$status,
    lines = normalize_bandwidth_progress_output(batch_lines),
    raw = c(batch_lines, res$output)
  )
}

test_that("bandwidth compaction never abbreviates the Bandwidth selection prefix", {
  compact_bandwidth <- getFromNamespace(".np_progress_compact_bandwidth_line", "npRmpi")

  compacted <- compact_bandwidth(
    "[npRmpi] Bandwidth selection (multistart 2/2, iteration 84, elapsed 10.0s, 99.9%, eta 0.0s)",
    max_width = 74L
  )

  expect_true(grepl("Bandwidth selection \\(", compacted))
  expect_false(grepl("Bandwidth sel \\(", compacted))
})

test_that("unknown-bound NOMAD restart detail drops synthetic percent and eta", {
  nomad_detail <- getFromNamespace(".np_nomad_progress_detail", "npRmpi")

  line <- nomad_detail(
    current_degree = 5L,
    best_record = list(degree = 1L),
    iteration = 77L,
    cumulative_iteration = 1234L,
    restart_index = 2L,
    nmulti = 2L,
    restart_durations = c(18.8),
    elapsed = 18.8
  )

  expect_identical(
    line,
    "multistart 2/2, iteration 77 (1234), elapsed 18.8s, deg (5), best (1)"
  )
  expect_false(grepl("%|eta ", line))
})

test_that("dark-launched bandwidth engine preserves nmulti=1 iteration heartbeats", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "npRmpi")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting density bandwidth", {
        activity_bw(28L)
        activity_bw(64L)
        7
      })
      expect_identical(value, 7)
    },
    now = progress_time_values(c(0, 1.2, 3.4, 5.6))
  )

  lines <- shadow_lines(actual)

  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection\\.\\.\\.$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(iteration 28, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(iteration 64, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_false(any(grepl("multistart", lines, fixed = TRUE)))
})

test_that("dark-launched bandwidth engine switches from iteration to estimate mode after first completion", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "npRmpi")
  set_total <- getFromNamespace(".np_progress_bandwidth_set_total", "npRmpi")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "npRmpi")
  step_bw <- getFromNamespace(".np_progress_bandwidth_multistart_step", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting density bandwidth", {
        set_total(2L)
        activity_bw(28L)
        step_bw(1L, 2L)
        activity_bw(56L)
        activity_bw(84L)
        11
      })
      expect_identical(value, 11)
    },
    now = progress_time_values(c(0, 1.0, 4.0, 7.0, 10.0, 12.0))
  )

  lines <- shadow_lines(actual)

  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 1/2, iteration 28, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, elapsed 4\\.0s, 50\\.0%, eta 4\\.0s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, iteration 56, elapsed 7\\.0s, 87\\.5%, eta 1\\.0s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, iteration 84, elapsed 10\\.0s, 99\\.9%, eta 0\\.0s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, elapsed 12\\.0s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("dark-launched completion-estimate heartbeats keep the unknown-total throttle cadence", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "npRmpi")
  set_total <- getFromNamespace(".np_progress_bandwidth_set_total", "npRmpi")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "npRmpi")
  step_bw <- getFromNamespace(".np_progress_bandwidth_multistart_step", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting regression bandwidth", {
        set_total(2L)
        activity_bw(20L)
        step_bw(1L, 2L)
        activity_bw(21L)
        activity_bw(22L)
        13
      })
      expect_identical(value, 13)
    },
    now = progress_time_values(c(0, 1.0, 4.0, 5.0, 6.2, 8.4))
  )

  complete_trace <- actual$trace[grepl("eta ", shadow_lines(actual), fixed = TRUE)]
  complete_lines <- vapply(complete_trace, `[[`, character(1L), "line")
  complete_times <- vapply(complete_trace, `[[`, numeric(1L), "now")

  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, elapsed 4\\.0s, 50\\.0%, eta 4\\.0s\\)$", complete_lines)))
  expect_false(any(grepl("elapsed 5\\.0s", complete_lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, iteration 22, elapsed 6\\.2s, 77\\.5%, eta 1\\.8s\\)$", complete_lines)))
  expect_true(all(diff(complete_times[seq_len(min(2L, length(complete_times)))]) >= 2.0))
})

test_that("bandwidth progress can carry a coordinator context label such as degree", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "npRmpi")
  set_total <- getFromNamespace(".np_progress_bandwidth_set_total", "npRmpi")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "npRmpi")
  step_bw <- getFromNamespace(".np_progress_bandwidth_multistart_step", "npRmpi")
  set_context <- getFromNamespace(".np_progress_bandwidth_set_context", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      set_context("deg (1,0)")
      on.exit(set_context(NULL), add = TRUE)
      value <- select_bw("Selecting regression bandwidth", {
        set_total(2L)
        activity_bw(20L)
        step_bw(1L, 2L)
        activity_bw(21L)
        13
      })
      expect_identical(value, 13)
    },
    now = progress_time_values(c(0, 1.0, 4.0, 6.2, 8.4))
  )

  lines <- shadow_lines(actual)

  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(deg \\(1,0\\)\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(deg \\(1,0\\), multistart 1/2, iteration 20, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(deg \\(1,0\\), multistart 2/2, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(deg \\(1,0\\), multistart 2/2, iteration 21, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(deg \\(1,0\\), multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("session selector routes avoid adjacent duplicate bandwidth lines", {
  actual <- run_bandwidth_session_contract()

  expect_identical(actual$status, 0L, info = paste(actual$raw, collapse = "\n"))
  expect_length(
    bandwidth_adjacent_duplicates(bandwidth_case_progress_lines(actual$lines, "npregbw")),
    0L
  )
  expect_length(
    bandwidth_adjacent_duplicates(bandwidth_case_progress_lines(actual$lines, "npcdensbw")),
    0L
  )
})

test_that("attach selector routes avoid adjacent duplicate bandwidth lines", {
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")

  actual <- run_bandwidth_attach_contract()

  expect_identical(actual$status, 0L, info = paste(actual$raw, collapse = "\n"))
  expect_length(
    bandwidth_adjacent_duplicates(bandwidth_case_progress_lines(actual$lines, "npregbw")),
    0L
  )
  expect_length(
    bandwidth_adjacent_duplicates(bandwidth_case_progress_lines(actual$lines, "npcdensbw")),
    0L
  )
})

test_that("profile selector routes avoid adjacent duplicate bandwidth lines", {
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_PROFILE_TEST"), "1"),
              "set NP_RMPI_ENABLE_PROFILE_TEST=1 to run profile-mode smoke")

  actual <- run_bandwidth_profile_contract()

  expect_identical(actual$status, 0L, info = paste(actual$raw, collapse = "\n"))
  expect_length(
    bandwidth_adjacent_duplicates(bandwidth_case_progress_lines(actual$lines, "npregbw")),
    0L
  )
  expect_length(
    bandwidth_adjacent_duplicates(bandwidth_case_progress_lines(actual$lines, "npcdensbw")),
    0L
  )
})
