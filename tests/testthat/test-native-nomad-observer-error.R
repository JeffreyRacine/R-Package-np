test_that("local native progress errors preserve their conditions after cleanup", {
  skip_if_not_installed("crs")
  pkg <- getNamespaceName(environment(npregbw))
  work <- tempfile("native-observer-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fixture <- normalizePath(test_path("fixtures", "native-nomad-observer.R"))
  input <- file.path(work, "input.rds")
  saveRDS(list(pkg = pkg, libpath = .libPaths()), input)
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
})
