local_npRmpi_fresh_subprocess_env <- function(extra = character()) {
  pkg.root <- tryCatch(
    normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
    error = function(e) ""
  )
  if (!nzchar(pkg.root))
    return(NULL)

  lib.path <- tempfile("npRmpi-northstar-lib-")
  dir.create(lib.path, recursive = TRUE, showWarnings = FALSE)

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-test-load", "-l", lib.path, pkg.root),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L
  if (status != 0L) {
    unlink(lib.path, recursive = TRUE, force = TRUE)
    return(NULL)
  }

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    "NP_RMPI_NO_REUSE_SLAVES=1",
    extra
  )
}

test_that("npRmpi rotate bootstrap npudens payload matches np north star", {
  skip_on_cran()
  skip_if_not_installed("np")

  env <- local_npRmpi_fresh_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  ok_tag <- "NPUDENS_ROTATE_NORTHSTAR_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(np))",
    "run_case <- function() {",
    "  options(np.messages = FALSE, plot.par.mfrow = FALSE)",
    "  set.seed(42)",
    "  n <- 200L",
    "  x <- runif(n, -1, 1)",
    "  y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))",
    "  fit_np <- np::npudens(~ y + x)",
    "  set.seed(4242)",
    "  out_np <- suppressWarnings(np:::plot.npdensity(",
    "    fit_np,",
    "    plot.behavior = 'data',",
    "    plot.errors.method = 'bootstrap',",
    "    view = 'rotate',",
    "    plot.errors.boot.method = 'inid',",
    "    plot.errors.boot.num = 39L,",
    "    plot.errors.type = 'all'",
    "  ))",
    "  suppressPackageStartupMessages(library(npRmpi))",
    "  options(npRmpi.autodispatch = TRUE)",
    "  npRmpi.init(nslaves = 1L, quiet = TRUE)",
    "  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "  fit_mpi <- npRmpi::npudens(~ y + x)",
    "  set.seed(4242)",
    "  out_mpi <- suppressWarnings(npRmpi:::plot.npdensity(",
    "    fit_mpi,",
    "    plot.behavior = 'data',",
    "    plot.errors.method = 'bootstrap',",
    "    view = 'rotate',",
    "    plot.errors.boot.method = 'inid',",
    "    plot.errors.boot.num = 39L,",
    "    plot.errors.type = 'all'",
    "  ))",
    "  stopifnot(length(out_np) == length(out_mpi))",
    "  for (nm in names(out_np)) {",
    "    stopifnot(isTRUE(all.equal(out_mpi[[nm]]$dens, out_np[[nm]]$dens, tolerance = 1e-8)))",
    "    rng_np <- range(out_np[[nm]]$derr, finite = TRUE)",
    "    rng_mpi <- range(out_mpi[[nm]]$derr, finite = TRUE)",
    "    stopifnot(length(rng_np) == 2L, length(rng_mpi) == 2L)",
    "    stopifnot(isTRUE(all.equal(rng_mpi, rng_np, tolerance = 0.05, scale = pmax(abs(rng_np), 1e-12))))",
    "  }",
    "}",
    "run_case()",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
