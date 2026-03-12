test_that("fixed conditional localpoly helper routes all counts sources through grouped core", {
  helper <- getFromNamespace(".np_inid_boot_from_conditional_localpoly_fixed", "npRmpi")
  body.txt <- paste(deparse(body(helper), width.cutoff = 500L), collapse = " ")

  expect_match(body.txt, "\\.np_inid_boot_from_conditional_localpoly_fixed_precompute\\(")
  expect_match(body.txt, "\\.np_inid_boot_from_conditional_localpoly_fixed_core\\(")
  expect_false(grepl("\\.np_inid_boot_from_conditional_localpoly_fixed_rowwise\\(", body.txt))
})

test_that("fixed conditional localpoly precompute uses direct fixed x-kernel weights in npRmpi", {
  precompute <- getFromNamespace(".np_inid_boot_from_conditional_localpoly_fixed_precompute", "npRmpi")
  body.txt <- paste(deparse(body(precompute), width.cutoff = 500L), collapse = " ")

  expect_match(body.txt, "\\.np_kernel_weights_direct\\(")
})

test_that("npRmpi fixed conditional grouped precompute/core smoke completes in subprocess", {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=FALSE, np.plot.inid.chunk.size=2L)",
      "np.ns <- asNamespace('npRmpi')",
      "pre <- get('.np_inid_boot_from_conditional_localpoly_fixed_precompute', envir=np.ns, inherits=FALSE)",
      "core <- get('.np_inid_boot_from_conditional_localpoly_fixed_core', envir=np.ns, inherits=FALSE)",
      "set.seed(6033222)",
      "n <- 68L",
      "tx <- data.frame(x = rnorm(n))",
      "ty <- data.frame(y = rnorm(n))",
      "x.grid <- seq(min(tx$x), max(tx$x), length.out = 7L)",
      "y.grid <- seq(min(ty$y), max(ty$y), length.out = 6L)",
      "grid <- expand.grid(y = y.grid, x = x.grid)",
      "ex <- data.frame(x = grid$x)",
      "ey <- data.frame(y = grid$y)",
      "B <- 7L",
      "counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))",
      "bw <- npcdensbw(xdat=tx, ydat=ty, bws=c(0.38, 0.38), bandwidth.compute=FALSE, bwtype='fixed', regtype='ll')",
      "state <- pre(xdat=tx, ydat=ty, exdat=ex, eydat=ey, bws=bw, cdf=FALSE)",
      "out <- core(state=state, B=B, counts=counts)",
      "stopifnot(is.list(out), is.matrix(out$t), length(out$t0) == nrow(ex), nrow(out$t) == B)",
      "cat('GROUPED_PRECORE_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("GROUPED_PRECORE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("npRmpi fixed conditional plot bootstrap covers fixed and geom in subprocess", {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260312)",
      "n <- 60L",
      "x <- rnorm(n)",
      "y <- sin(x) + rnorm(n, sd=0.25)",
      "bw.dens <- npcdensbw(xdat=data.frame(x=x), ydat=data.frame(y=y), regtype='lp', degree=2, bwtype='fixed', bandwidth.compute=FALSE, bws=c(0.45, 0.45), basis='glp', bernstein.basis=FALSE)",
      "fit.dens <- npcdens(bws=bw.dens)",
      "bw.dist <- npcdistbw(xdat=data.frame(x=x), ydat=data.frame(y=y), regtype='lp', degree=2, bwtype='fixed', bandwidth.compute=FALSE, bws=c(0.45, 0.45), basis='glp', bernstein.basis=FALSE)",
      "fit.dist <- npcdist(bws=bw.dist)",
      "for (boot.method in c('fixed', 'geom')) {",
      "  out.dens <- plot(fit.dens, plot.errors.method='bootstrap', plot.errors.boot.method=boot.method, plot.errors.boot.num=9L, plot.errors.boot.blocklength=3L, plot.behavior='data')",
      "  out.dist <- plot(fit.dist, plot.errors.method='bootstrap', plot.errors.boot.method=boot.method, plot.errors.boot.num=9L, plot.errors.boot.blocklength=3L, plot.behavior='data')",
      "  stopifnot(is.list(out.dens), length(out.dens) > 0L, is.list(out.dist), length(out.dist) > 0L)",
      "  cat('PLOT_BOOT_OK', boot.method, '\\n')",
      "}"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PLOT_BOOT_OK fixed", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PLOT_BOOT_OK geom", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
