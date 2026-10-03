test_that("witnessed native errors and flag interrupts cannot poison NOMAD", {
  skip_if_not_installed("crs")
  # R_PolledEvents is a Unix API. Windows GUI interruption needs its platform
  # release lane; do not pretend that a pre-entry timer tests a native callback.
  skip_on_os("windows")
  pkg <- getNamespaceName(environment(npregbw))
  work <- tempfile("native-unwind-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  source <- normalizePath(test_path("fixtures", "native-nomad-unwind.c"))
  fixture <- normalizePath(test_path("fixtures", "native-nomad-unwind.R"))
  file.copy(source, file.path(work, basename(source)))
  old <- setwd(work)
  on.exit(setwd(old), add = TRUE)
  build <- suppressWarnings(system2(file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "native-nomad-unwind.c"), stdout = TRUE, stderr = TRUE))
  status <- attr(build, "status")
  expect_true(is.null(status) || identical(status, 0L),
              info = paste(build, collapse = "\n"))
  dll <- file.path(work, paste0("native-nomad-unwind", .Platform$dynlib.ext))
  if (!file.exists(dll)) return(invisible(NULL))
  input <- file.path(work, "input.rds")
  saveRDS(list(pkg = pkg, dll = dll, libpath = .libPaths()), input)
  lines <- c(paste0("spec <- readRDS(", deparse(input), ")"),
    ".libPaths(spec$libpath)", paste0("source(", deparse(fixture), ")"))
  if (pkg == "npRmpi") {
    env <- npRmpi_subprocess_env()
    skip_if(is.null(env))
    child <- npRmpi_run_rscript_subprocess(lines, timeout = 30L,
                                         env = env, cleanup = FALSE)
  } else {
    script <- file.path(work, "probe.R")
    writeLines(lines, script)
    log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE, timeout = 30L))
    status <- attr(log, "status")
    child <- list(status = if (is.null(status)) 0L else status, output = log)
  }
  expect_identical(child$status, 0L, info = paste(child$output, collapse = "\n"))
  expect_true(any(grepl("NATIVE_UNWIND_PASS", child$output, fixed = TRUE)),
              info = paste(child$output, collapse = "\n"))
})
