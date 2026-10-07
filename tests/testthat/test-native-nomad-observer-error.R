test_that("local native progress errors preserve their conditions after cleanup", {
  skip_if_not_installed("crs")
  pkg <- getNamespaceName(environment(npregbw))
  work <- tempfile("native-observer-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fixture <- normalizePath(test_path("fixtures", "native-nomad-observer.R"))
  input <- file.path(work, "input.rds")
  spec <- list(pkg = pkg, libpath = .libPaths())
  # Exercise R's real message-less interrupt through its public API on both
  # Windows and Unix. This tests delivery through the native observer boundary,
  # not OS-specific keyboard or GUI event handling.
  # Calibrate R's complete toolchain, including configured wrappers/flags.
  # An unavailable toolchain must not suppress compiler-independent checks.
  writeLines("void np_observer_toolchain_probe(void) {}",
             file.path(work, "toolchain-probe.c"))
  old <- setwd(work)
  probe <- tryCatch(suppressWarnings(system2(file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "toolchain-probe.c"), stdout = TRUE, stderr = TRUE)),
    error = function(e) structure(conditionMessage(e), status = 1L))
  setwd(old)
  compiler.available <- is.null(attr(probe, "status")) &&
    file.exists(file.path(work, paste0("toolchain-probe", .Platform$dynlib.ext)))
  if (compiler.available) {
    source <- normalizePath(test_path("fixtures", "native-nomad-interrupt.c"))
    file.copy(source, file.path(work, basename(source)))
    old <- setwd(work)
    build <- suppressWarnings(system2(file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", "native-nomad-interrupt.c"), stdout = TRUE, stderr = TRUE))
    setwd(old)
    status <- attr(build, "status")
    expect_true(is.null(status) || identical(status, 0L),
                info = paste(build, collapse = "\n"))
    spec$dll <- file.path(work, paste0("native-nomad-interrupt", .Platform$dynlib.ext))
    if (!file.exists(spec$dll)) spec$dll <- NULL
  }
  saveRDS(spec, input)
  lines <- c(paste0("spec <- readRDS(", deparse(input), ")"),
    ".libPaths(spec$libpath)", paste0("source(", deparse(fixture), ")"))
  if (pkg == "npRmpi") {
    env <- npRmpi_subprocess_env()
    skip_if(is.null(env))
    child <- npRmpi_run_rscript_subprocess(lines, timeout = 45L,
                                         env = env, cleanup = FALSE)
  } else {
    script <- file.path(work, "probe.R")
    writeLines(lines, script)
    log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE, timeout = 45L))
    status <- attr(log, "status")
    child <- list(status = if (is.null(status)) 0L else status, output = log)
  }
  expect_identical(child$status, 0L, info = paste(child$output, collapse = "\n"))
  expect_true(any(grepl("NATIVE_OBSERVER_PASS", child$output, fixed = TRUE)),
              info = paste(child$output, collapse = "\n"))
  if (!is.null(spec$dll)) {
    # An ordinary error handler must not swallow Stop. These children abort at
    # top level; the explicit-interrupt-handler/recovery cases ran above.
    fixture <- normalizePath(test_path("fixtures", "native-nomad-real-interrupt.R"))
    for (wrapper in c("top", "try", "error", "error-hook", "interrupt-hook", "traceback", "traceback-empty", "resume")) {
      spec$wrap <- if (wrapper == "traceback-empty") "traceback" else wrapper
      saveRDS(spec, input)
      lines <- c(paste0("spec <- readRDS(", deparse(input), ")"),
        ".libPaths(spec$libpath)", paste0("source(", deparse(fixture), ")"))
      # A top-level interrupt unwinds source() itself. Inspect its traceback in
      # the next top-level expression, as in an ordinary R script.
      if (wrapper == "traceback-empty") lines <- c(lines,
        "if (exists('.Traceback', envir = baseenv(), inherits = FALSE)) assign('.Traceback', NULL, envir = baseenv())")
      if (wrapper %in% c("traceback", "traceback-empty")) lines <- c(lines,
        "if (!length(.traceback())) q(save = 'no', status = 4L)",
        "cat('REAL_INTERRUPT_TRACEBACK_PASS\n')",
        "if (spec$pkg == 'npRmpi') get('mpi.finalize', asNamespace(spec$pkg))()")
      if (pkg == "npRmpi") {
        child <- npRmpi_run_rscript_subprocess(lines, timeout = 30L,
                                             env = env, cleanup = FALSE)
      } else {
        writeLines(lines, script)
        log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
          c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE, timeout = 30L))
        status <- attr(log, "status")
        child <- list(status = if (is.null(status)) 0L else status, output = log)
      }
      info <- paste(wrapper, paste(child$output, collapse = "\n"))
      expected <- switch(wrapper, `error-hook` = 3L, traceback = 0L,
                         `traceback-empty` = 4L, 1L)
      expect_identical(child$status, expected, info = info)
      if (wrapper %in% c("error-hook", "traceback", "traceback-empty"))
        expect_true(any(grepl("REAL_INTERRUPT_ERROR_HOOK", child$output, fixed = TRUE)), info = info)
      if (wrapper == "resume")
        expect_true(any(grepl("cannot be resumed", child$output, fixed = TRUE)), info = info)
      if (wrapper == "interrupt-hook")
        expect_true(any(grepl("REAL_INTERRUPT_HOOK", child$output, fixed = TRUE)), info = info)
      if (wrapper == "traceback-empty")
        expect_false(any(grepl("REAL_INTERRUPT_TRACEBACK_PASS", child$output, fixed = TRUE)), info = info)
      if (wrapper == "traceback")
        expect_true(any(grepl("REAL_INTERRUPT_TRACEBACK_PASS", child$output, fixed = TRUE)), info = info)
      expect_true(any(grepl("REAL_OBSERVER_INTERRUPT", child$output, fixed = TRUE)), info = info)
      # q() in the user's error hook exits before ordinary R on.exit handlers,
      # just as it does for an interrupt outside a search.
      if (wrapper != "error-hook")
        expect_true(any(grepl("REAL_OBSERVER_CALLER_UNWOUND", child$output, fixed = TRUE)), info = info)
      expect_false(any(grepl("REAL_OBSERVER_INCORRECTLY_CONTINUED", child$output, fixed = TRUE)), info = info)
      expect_false(any(grepl("bad error message", child$output, fixed = TRUE)), info = info)
    }
  }
  if (!compiler.available)
    skip("Compiler-independent observer checks passed; native interrupt fixture requires R's configured C compiler")
})
