test_that("npscoef inid plot bootstrap uses all slaves plus master assist", {
  skip_on_cran()

  trace.file <- tempfile("npscoef-plot-master-assist-", fileext = ".tsv")
  env <- npRmpi_subprocess_env(c(
    "NP_RMPI_NO_REUSE_SLAVES=1",
    paste0("NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE=", trace.file)
  ))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  ok_tag <- "NPSCOEF_PLOT_MASTER_ASSIST_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "run_case <- function() {",
    "  npRmpi.init(nslaves = 2, quiet = TRUE)",
    "  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
    "  options(np.plot.inid.chunk.size = 10L)",
    "  set.seed(7302)",
    "  n <- 28L",
    "  xdat <- data.frame(x = rnorm(n))",
    "  zdat <- data.frame(z = runif(n))",
    "  y <- (1 + 0.5 * xdat$x) * sin(zdat$z) + rnorm(n, sd = 0.1)",
    "  bw <- npscoefbw(xdat = xdat, zdat = zdat, ydat = y,",
    "                  bws = 0.5, bwtype = 'fixed',",
    "                  bandwidth.compute = FALSE, regtype = 'll')",
    "  boot.plot <- NULL",
    "  invisible(capture.output(boot.plot <- plot(",
    "    bw, xdat = xdat, ydat = y, zdat = zdat, neval = 4L,",
    "    coef = FALSE, behavior = 'data', errors = 'bootstrap',",
    "    bootstrap = 'inid', B = 31L, band = 'pointwise'",
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
    grepl("what=inid-scoef-localpoly", trace, fixed = TRUE) &
      grepl("event=fanout.start", trace, fixed = TRUE)
  ]
  expect_true(length(starts) > 0L, info = paste(trace, collapse = "\n"))
  expect_true(all(grepl("workers=2", starts, fixed = TRUE)),
              info = paste(starts, collapse = "\n"))
  expect_true(all(grepl("master_local_chunk=TRUE", starts, fixed = TRUE)),
              info = paste(starts, collapse = "\n"))
  expect_true(any(grepl("what=inid-scoef-localpoly", trace, fixed = TRUE) &
                    grepl("event=fanout.send.initial", trace, fixed = TRUE) &
                    grepl("dest=1", trace, fixed = TRUE)),
              info = paste(trace, collapse = "\n"))
  expect_true(any(grepl("what=inid-scoef-localpoly", trace, fixed = TRUE) &
                    grepl("event=fanout.send.initial", trace, fixed = TRUE) &
                    grepl("dest=2", trace, fixed = TRUE)),
              info = paste(trace, collapse = "\n"))
  expect_true(any(grepl("what=inid-scoef-localpoly", trace, fixed = TRUE) &
                    grepl("event=fanout.master_local_chunk.done", trace, fixed = TRUE)),
              info = paste(trace, collapse = "\n"))
})
