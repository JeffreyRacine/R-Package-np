test_that("an R time-limit error cannot poison a subsequent native NOMAD solve", {
  skip_if_not_installed("crs")
  pkg <- getNamespaceName(environment(npregbw))
  runner <- function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    options(np.messages = FALSE, np.tree = FALSE)
    run <- function(fun) {
      if (pkg == "npRmpi")
        get(".npRmpi_with_local_regression", asNamespace(pkg))(fun())
      else fun()
    }
    set.seed(62031)
    small <- data.frame(x = stats::runif(60), z = stats::runif(60))
    small$y <- sin(4 * small$x) + small$z + stats::rnorm(60, sd = .2)
    big <- small[rep(seq_len(nrow(small)), length.out = 2400), ]
    big$x <- big$x + seq_len(nrow(big)) * 1e-6
    signature <- function(b) list(b$bw, b$xbw, b$ybw, b$fval, b$num.feval)
    results <- list()
    for (family in c("npregbw", "npcdensbw")) {
      fun <- get(family, asNamespace(pkg))
      search <- function(dat, budget) run(function() {
        fun(y ~ x + z, data = dat, regtype = "lp", degree = c(1L, 1L),
            bwsolver = "mads", nmulti = 1L,
            nomad.opts = list(MAX_BB_EVAL = budget))
      })
      set.seed(42)
      before <- signature(search(small, 12L))
      caught <- tryCatch({
        setTimeLimit(elapsed = .01, transient = TRUE)
        search(big, 4000L)
      }, error = identity, finally = setTimeLimit(cpu = Inf, elapsed = Inf))
      set.seed(42)
      after <- signature(search(small, 12L))
      results[[family]] <- list(is_error = inherits(caught, "error"),
        message = if (inherits(caught, "error")) conditionMessage(caught) else "",
        before = before, after = after)
    }
    if (pkg == "npRmpi") get("mpi.finalize", asNamespace(pkg))()
    results
  }
  environment(runner) <- baseenv()
  input <- tempfile("native-error-input-", fileext = ".rds")
  output <- tempfile("native-error-output-", fileext = ".rds")
  script <- tempfile("native-error-script-", fileext = ".R")
  on.exit(unlink(c(input, output, script)), add = TRUE)
  saveRDS(list(run = runner, pkg = pkg, libpath = .libPaths()), input)
  lines <- c(paste0("spec <- readRDS(", deparse(input), ")"),
             ".libPaths(spec$libpath)",
             paste0("saveRDS(spec$run(spec$pkg), ", deparse(output), ")"))
  if (pkg == "npRmpi") {
    env <- npRmpi_subprocess_env()
    skip_if(is.null(env))
    child <- npRmpi_run_rscript_subprocess(lines, timeout = 30L,
                                         env = env, cleanup = FALSE)
  } else {
    writeLines(lines, script)
    log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE, timeout = 30L))
    status <- attr(log, "status")
    child <- list(status = if (is.null(status)) 0L else status, output = log)
  }
  expect_equal(child$status, 0L, info = paste(child$output, collapse = "\n"))
  expect_true(file.exists(output), info = paste(child$output, collapse = "\n"))
  if (!file.exists(output)) return(invisible(NULL))
  result <- readRDS(output)
  for (family in names(result)) {
    expect_true(result[[family]]$is_error, info = family)
    expect_match(result[[family]]$message, "time limit", fixed = TRUE)
    expect_identical(result[[family]]$after, result[[family]]$before, info = family)
  }
})
