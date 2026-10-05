test_that("local native progress errors preserve their conditions after cleanup", {
  skip_if_not_installed("crs")
  pkg <- getNamespaceName(environment(npregbw))
  work <- tempfile("native-observer-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fixture <- normalizePath(test_path("fixtures", "native-nomad-observer.R"))
  input <- file.path(work, "input.rds")
  spec <- list(pkg = pkg, libpath = .libPaths())
  # The existing small child also protects real GUI-style message-less
  # interrupts. Keep platform-specific event-loop injection out of Windows.
  if (.Platform$OS.type == "unix") {
    source <- normalizePath(test_path("fixtures", "native-nomad-unwind.c"))
    file.copy(source, file.path(work, basename(source)))
    old <- setwd(work)
    build <- suppressWarnings(system2(file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", "native-nomad-unwind.c"), stdout = TRUE, stderr = TRUE))
    setwd(old)
    status <- attr(build, "status")
    expect_true(is.null(status) || identical(status, 0L),
                info = paste(build, collapse = "\n"))
    spec$dll <- file.path(work, paste0("native-nomad-unwind", .Platform$dynlib.ext))
    if (!file.exists(spec$dll)) return(invisible(NULL))
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
    for (wrapper in c("top", "try", "error")) {
      spec$wrap <- wrapper
      saveRDS(spec, input)
      lines <- c(paste0("spec <- readRDS(", deparse(input), ")"),
        ".libPaths(spec$libpath)", paste0("source(", deparse(fixture), ")"))
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
      expect_identical(child$status, 1L, info = info)
      expect_true(any(grepl("REAL_OBSERVER_INTERRUPT", child$output, fixed = TRUE)), info = info)
      expect_true(any(grepl("REAL_OBSERVER_CALLER_UNWOUND", child$output, fixed = TRUE)), info = info)
      expect_false(any(grepl("REAL_OBSERVER_INCORRECTLY_CONTINUED", child$output, fixed = TRUE)), info = info)
      expect_false(any(grepl("bad error message", child$output, fixed = TRUE)), info = info)
    }
  }
})
