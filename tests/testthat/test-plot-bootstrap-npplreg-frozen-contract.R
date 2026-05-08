npplreg_semiparam_frozen_contract_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  ok_tag <- "NPPLREG_SEMIPARAM_FROZEN_CONTRACT_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(np))",
    "suppressPackageStartupMessages(library(npRmpi))",
    "run_case <- function() {",
    "  npRmpi.init(nslaves = 1, quiet = TRUE)",
    "  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
    "  txdat <- data.frame(x = c(0.05, 0.05, 0.25, 0.25, 0.60, 0.60, 0.90, 0.90))",
    "  tzdat <- data.frame(z = c(0.10, 0.10, 0.35, 0.35, 0.65, 0.65, 0.95, 0.95))",
    "  y <- with(txdat, 0.5 + 0.8 * x + sin(2 * pi * tzdat$z))",
    "  bw <- npplregbw(xdat = txdat, zdat = tzdat, ydat = y, bws = matrix(c(3L, 3L), nrow = 2L, ncol = 1L), bwtype = 'adaptive_nn', bandwidth.compute = FALSE, regtype = 'll')",
    "  ns <- asNamespace('npRmpi')",
    "  orig <- getFromNamespace('.np_inid_boot_from_plreg', 'npRmpi')",
    "  local.fit <- getFromNamespace('.np_plot_plreg_local_fit', 'npRmpi')",
    "  boot.method <- getS3method('compute.bootstrap.errors', 'plbandwidth', envir = ns)",
    "  modes <- character()",
    "  seen.ex <- list()",
    "  seen.ez <- list()",
    "  seen.helper.t0 <- list()",
    "  assignInNamespace('.np_inid_boot_from_plreg', function(..., mode = c('exact', 'frozen')) {",
    "    mode <- match.arg(mode)",
    "    dots <- list(...)",
    "    modes <<- c(modes, mode)",
    "    seen.ex[[length(seen.ex) + 1L]] <<- dots$exdat",
    "    seen.ez[[length(seen.ez) + 1L]] <<- dots$ezdat",
    "    seen.helper.t0[length(seen.helper.t0) + 1L] <<- list(dots$t0.override)",
    "    t0.local <- if (is.null(dots$t0.override)) rep(0, nrow(dots$exdat)) else as.double(dots$t0.override)",
    "    B.local <- if (is.null(dots$B)) 1L else as.integer(dots$B)",
    "    list(t = matrix(rep(t0.local, each = B.local), nrow = B.local), t0 = t0.local)",
    "  }, ns = 'npRmpi')",
    "  on.exit(assignInNamespace('.np_inid_boot_from_plreg', orig, ns = 'npRmpi'), add = TRUE)",
    "  bootstrap.plot <- NULL",
    "  invisible(capture.output(bootstrap.plot <- plot(bw, xdat = txdat, ydat = y, zdat = tzdat, neval = 3L, output = 'data', errors = 'bootstrap', bootstrap = 'inid', boot_control = np_boot_control(nonfixed = 'frozen'), B = 41L, band = 'pointwise')))",
    "  stopifnot(length(modes) >= 1L)",
    "  stopifnot(all(modes == 'frozen'))",
    "  expected.slice <- local.fit(bws = bw, xdat = txdat, ydat = y, zdat = tzdat, exdat = seen.ex[[1L]], ezdat = seen.ez[[1L]])",
    "  center.check <- boot.method(xdat = txdat, ydat = y, zdat = tzdat, exdat = seen.ex[[1L]], ezdat = seen.ez[[1L]], gradients = FALSE, slice.index = 1L, plot.errors.boot.method = 'inid', t0.override = as.vector(expected.slice$mean), plot.errors.boot.nonfixed = 'frozen', plot.errors.boot.blocklen = 1L, plot.errors.boot.num = 41L, plot.errors.center = 'estimate', plot.errors.type = 'pointwise', plot.errors.alpha = 0.05, bws = bw)",
    "  expected.err <- cbind(as.vector(expected.slice$mean), -as.vector(expected.slice$mean), rep(NA_real_, length(expected.slice$mean)))",
    "  stopifnot(isTRUE(all.equal(center.check$boot.err, expected.err, tolerance = 1e-10)))",
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
}

test_that("npRmpi npplreg frozen plot route forwards the frozen mode and center", {
  npplreg_semiparam_frozen_contract_case()
})
