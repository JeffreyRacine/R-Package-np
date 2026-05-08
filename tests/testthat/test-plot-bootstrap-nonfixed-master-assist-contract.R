test_that("nonfixed exact and frozen plot bootstraps use all slaves plus master assist", {
  skip_on_cran()

  trace.dir <- tempfile("plot-nonfixed-master-assist-")
  dir.create(trace.dir, recursive = TRUE)
  env <- npRmpi_subprocess_env("NP_RMPI_NO_REUSE_SLAVES=1")
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  nslaves <- 3L
  ok_tag <- "PLOT_NONFIXED_MASTER_ASSIST_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    sprintf("trace.dir <- %s", deparse(trace.dir)),
    "run_case <- function(name, expr) {",
    "  trace <- file.path(trace.dir, paste0(name, '.tsv'))",
    "  options(npRmpi.bootstrap.transport.trace.file = trace)",
    "  on.exit(options(npRmpi.bootstrap.transport.trace.file = ''), add = TRUE)",
    "  invisible(capture.output(invisible(force(expr))))",
    "}",
    sprintf("npRmpi.init(nslaves = %d, quiet = TRUE)", nslaves),
    "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
    "options(np.plot.inid.chunk.size = 10L)",
    "set.seed(7303)",
    "n <- 40L",
    "x <- runif(n, -1, 1)",
    "z <- runif(n, -1, 1)",
    "y <- sin(x) + 0.5 * z + rnorm(n, sd = 0.2)",
    "run_case('npudens_exact', {",
    "  bw <- npudensbw(dat = data.frame(x = x), bws = 5L,",
    "                  bwtype = 'generalized_nn', bandwidth.compute = FALSE)",
    "  plot(bw, xdat = data.frame(x = x), neval = 8L, output = 'data',",
    "       errors = 'bootstrap', bootstrap = 'inid', B = 43L, band = 'pointwise')",
    "})",
    "run_case('npudist_exact', {",
    "  bw <- npudistbw(dat = data.frame(x = x), bws = 5L,",
    "                  bwtype = 'generalized_nn', bandwidth.compute = FALSE)",
    "  plot(bw, xdat = data.frame(x = x), neval = 8L, output = 'data',",
    "       errors = 'bootstrap', bootstrap = 'inid', B = 43L, band = 'pointwise')",
    "})",
    "run_case('npindex_exact', {",
    "  bw <- npindexbw(xdat = data.frame(x = x, z = z), ydat = y,",
    "                  bws = c(1, 0.5, 5L), bwtype = 'generalized_nn',",
    "                  bandwidth.compute = FALSE, regtype = 'll')",
    "  plot(bw, xdat = data.frame(x = x, z = z), ydat = y, neval = 8L,",
    "       output = 'data', gradients = FALSE, errors = 'bootstrap',",
    "       bootstrap = 'inid', B = 43L, band = 'pointwise')",
    "})",
    "run_case('npcdens_frozen', {",
    "  bw <- npcdensbw(xdat = data.frame(x = x), ydat = data.frame(y = y),",
    "                  bws = c(5L, 5L), bwtype = 'generalized_nn',",
    "                  bandwidth.compute = FALSE, regtype = 'lc')",
    "  plot(bw, xdat = data.frame(x = x), ydat = data.frame(y = y),",
    "       neval = 8L, output = 'data', errors = 'bootstrap',",
    "       bootstrap = 'inid', boot_control = np_boot_control(nonfixed = 'frozen'),",
    "       B = 43L, band = 'pointwise')",
    "})",
    "run_case('npcdist_frozen', {",
    "  bw <- npcdistbw(xdat = data.frame(x = x), ydat = data.frame(y = y),",
    "                  bws = c(5L, 5L), bwtype = 'generalized_nn',",
    "                  bandwidth.compute = FALSE, regtype = 'lc')",
    "  plot(bw, xdat = data.frame(x = x), ydat = data.frame(y = y),",
    "       neval = 8L, output = 'data', errors = 'bootstrap',",
    "       bootstrap = 'inid', boot_control = np_boot_control(nonfixed = 'frozen'),",
    "       B = 43L, band = 'pointwise')",
    "})",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(lines = lines, timeout = 240L, env = env)

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))

  expected <- c(
    npudens_exact = "inid-ksum-unconditional-exact",
    npudist_exact = "inid-ksum-unconditional-exact",
    npindex_exact = "inid-index-exact",
    npcdens_frozen = "inid-hat-frozen-conditional",
    npcdist_frozen = "inid-hat-frozen-conditional"
  )
  for (case in names(expected)) {
    trace.file <- file.path(trace.dir, paste0(case, ".tsv"))
    expect_true(file.exists(trace.file), info = case)
    trace <- readLines(trace.file, warn = FALSE)
    starts <- trace[
      grepl(paste0("what=", expected[[case]]), trace, fixed = TRUE) &
        grepl("event=fanout.start", trace, fixed = TRUE)
    ]
    expect_true(length(starts) > 0L, info = paste(trace, collapse = "\n"))
    expect_true(all(grepl(sprintf("workers=%d", nslaves), starts, fixed = TRUE)),
                info = paste(starts, collapse = "\n"))
    expect_true(all(grepl("master_local_chunk=TRUE", starts, fixed = TRUE)),
                info = paste(starts, collapse = "\n"))
    assists <- trace[
      grepl(paste0("what=", expected[[case]]), trace, fixed = TRUE) &
        grepl("event=fanout.master_assist.start", trace, fixed = TRUE)
    ]
    expect_true(length(assists) > 0L, info = paste(trace, collapse = "\n"))
    expect_true(all(grepl("scheduler=static_bundle", assists, fixed = TRUE)),
                info = paste(assists, collapse = "\n"))
    for (dest in seq_len(nslaves)) {
      expect_true(any(grepl(paste0("what=", expected[[case]]), trace, fixed = TRUE) &
                        grepl("event=fanout.send.initial", trace, fixed = TRUE) &
                        grepl(sprintf("dest=%d", dest), trace, fixed = TRUE)),
                  info = paste(trace, collapse = "\n"))
    }
    expect_true(any(grepl(paste0("what=", expected[[case]]), trace, fixed = TRUE) &
                      grepl("event=fanout.master_local_chunk.done", trace, fixed = TRUE)),
                info = paste(trace, collapse = "\n"))
  }
})
