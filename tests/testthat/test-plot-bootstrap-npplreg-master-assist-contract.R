test_that("npplreg LL/LP inid plot bootstrap uses fused master-assist fanout", {
  skip_on_cran()

  trace.file <- tempfile("npplreg-plot-master-assist-", fileext = ".tsv")
  env <- npRmpi_subprocess_env(c(
    "NP_RMPI_NO_REUSE_SLAVES=1",
    paste0("NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE=", trace.file)
  ))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  nslaves <- 3L
  ok_tag <- "NPPLREG_PLOT_MASTER_ASSIST_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "run_case <- function() {",
    sprintf("  npRmpi.init(nslaves = %d, quiet = TRUE)", nslaves),
    "  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
    "  set.seed(7301)",
    "  n <- 24L",
    "  xdat <- data.frame(x = rnorm(n))",
    "  zdat <- data.frame(z = runif(n))",
    "  y <- 1 + 0.5 * xdat$x + sin(zdat$z) + rnorm(n, sd = 0.1)",
    "  bw <- npplregbw(xdat = xdat, zdat = zdat, ydat = y,",
    "                  bws = matrix(c(0.5, 0.5), nrow = 2L, ncol = 1L),",
    "                  bwtype = 'fixed', bandwidth.compute = FALSE,",
    "                  regtype = 'll')",
    "  boot.plot <- NULL",
    "  invisible(capture.output(boot.plot <- plot(",
    "    bw, xdat = xdat, ydat = y, zdat = zdat, neval = 4L,",
    "    output = 'data', errors = 'bootstrap', bootstrap = 'inid',",
    "    B = 41L, band = 'pointwise'",
    "  )))",
    "}",
    "run_case()",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(lines = lines, timeout = 180L, env = env)

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(file.exists(trace.file))

  trace <- readLines(trace.file, warn = FALSE)
  starts <- trace[
    grepl("what=inid-plreg-fixed-fused", trace, fixed = TRUE) &
      grepl("event=fanout.start", trace, fixed = TRUE)
  ]
  expect_true(length(starts) > 0L, info = paste(trace, collapse = "\n"))
  expect_true(all(grepl("master_local_chunk=TRUE", starts, fixed = TRUE)),
              info = paste(starts, collapse = "\n"))
  expect_true(all(grepl(sprintf("workers=%d", nslaves), starts, fixed = TRUE)),
              info = paste(starts, collapse = "\n"))
  assists <- trace[
    grepl("what=inid-plreg-fixed-fused", trace, fixed = TRUE) &
      grepl("event=fanout.master_assist.start", trace, fixed = TRUE)
  ]
  expect_true(length(assists) > 0L, info = paste(trace, collapse = "\n"))
  expect_true(all(grepl("scheduler=static_bundle", assists, fixed = TRUE)),
              info = paste(assists, collapse = "\n"))
  for (dest in seq_len(nslaves)) {
    expect_true(any(grepl("what=inid-plreg-fixed-fused", trace, fixed = TRUE) &
                      grepl("event=fanout.send.initial", trace, fixed = TRUE) &
                      grepl(sprintf("dest=%d", dest), trace, fixed = TRUE)),
                info = paste(trace, collapse = "\n"))
  }
  expect_true(any(grepl("what=inid-plreg-fixed-fused", trace, fixed = TRUE) &
                    grepl("event=fanout.master_local_chunk.done", trace, fixed = TRUE)),
              info = paste(trace, collapse = "\n"))
})

test_that("npplreg LC inid plot bootstrap uses fused master-assist fanout", {
  skip_on_cran()

  trace.file <- tempfile("npplreg-lc-plot-master-assist-", fileext = ".tsv")
  env <- npRmpi_subprocess_env(c(
    "NP_RMPI_NO_REUSE_SLAVES=1",
    paste0("NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE=", trace.file)
  ))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  nslaves <- 3L
  ok_tag <- "NPPLREG_LC_PLOT_MASTER_ASSIST_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "run_case <- function() {",
    sprintf("  npRmpi.init(nslaves = %d, quiet = TRUE)", nslaves),
    "  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
    "  set.seed(7302)",
    "  n <- 24L",
    "  xdat <- data.frame(x = rnorm(n))",
    "  zdat <- data.frame(z = runif(n))",
    "  y <- 1 + 0.5 * xdat$x + sin(zdat$z) + rnorm(n, sd = 0.1)",
    "  bw <- npplregbw(xdat = xdat, zdat = zdat, ydat = y,",
    "                  bws = matrix(c(0.5, 0.5), nrow = 2L, ncol = 1L),",
    "                  bwtype = 'fixed', bandwidth.compute = FALSE,",
    "                  regtype = 'lc')",
    "  boot.plot <- NULL",
    "  invisible(capture.output(boot.plot <- plot(",
    "    bw, xdat = xdat, ydat = y, zdat = zdat, neval = 4L,",
    "    output = 'data', errors = 'bootstrap', bootstrap = 'inid',",
    "    B = 41L, band = 'pointwise'",
    "  )))",
    "}",
    "run_case()",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(lines = lines, timeout = 180L, env = env)

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(file.exists(trace.file))

  trace <- readLines(trace.file, warn = FALSE)
  starts <- trace[
    grepl("what=inid-plreg-lc-fixed-fused", trace, fixed = TRUE) &
      grepl("event=fanout.start", trace, fixed = TRUE)
  ]
  expect_true(length(starts) > 0L, info = paste(trace, collapse = "\n"))
  expect_true(all(grepl("master_local_chunk=TRUE", starts, fixed = TRUE)),
              info = paste(starts, collapse = "\n"))
  expect_true(all(grepl(sprintf("workers=%d", nslaves), starts, fixed = TRUE)),
              info = paste(starts, collapse = "\n"))
  assists <- trace[
    grepl("what=inid-plreg-lc-fixed-fused", trace, fixed = TRUE) &
      grepl("event=fanout.master_assist.start", trace, fixed = TRUE)
  ]
  expect_true(length(assists) > 0L, info = paste(trace, collapse = "\n"))
  expect_true(all(grepl("scheduler=static_bundle", assists, fixed = TRUE)),
              info = paste(assists, collapse = "\n"))
  for (dest in seq_len(nslaves)) {
    expect_true(any(grepl("what=inid-plreg-lc-fixed-fused", trace, fixed = TRUE) &
                      grepl("event=fanout.send.initial", trace, fixed = TRUE) &
                      grepl(sprintf("dest=%d", dest), trace, fixed = TRUE)),
                info = paste(trace, collapse = "\n"))
  }
  expect_true(any(grepl("what=inid-plreg-lc-fixed-fused", trace, fixed = TRUE) &
                    grepl("event=fanout.master_local_chunk.done", trace, fixed = TRUE)),
              info = paste(trace, collapse = "\n"))
})
