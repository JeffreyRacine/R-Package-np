test_that("extended generalized-NN regression fits preserve training identity", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess proof")

  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE, np.tree=FALSE, np.extendednn=TRUE)",
    "npRmpi.init(nslaves=1L, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "x <- data.frame(x=c(-1.8,-0.9,-0.25,0.1,0.65,1.4,2.3))",
    "y <- sin(1.3*x$x)+seq_len(nrow(x))/50",
    "k <- nrow(x)+2L",
    "bw <- npregbw(xdat=x, ydat=y, bws=k, bwmethod='cv.ls', bwtype='generalized_nn', bwscaling=FALSE, regtype='lc', ckertype='gaussian', bandwidth.compute=FALSE)",
    "expected <- vapply(seq_len(nrow(x)), function(index) {",
    "  distance <- sort(abs(x$x[-index]-x$x[[index]]), method='radix')",
    "  radius <- max(distance)*k/length(distance)",
    "  weight <- dnorm((x$x[[index]]-x$x)/radius)/radius",
    "  sum(weight*y)/sum(weight)",
    "}, numeric(1L))",
    "for (tree in c(FALSE,TRUE)) {",
    "  options(np.tree=tree)",
    "  stopifnot(isTRUE(all.equal(fitted(npreg(bws=bw, txdat=x, tydat=y)), expected, tolerance=2e-12)))",
    "}",
    "cat('NPRMPI_REGRESSION_EXTENDED_GNN_OK\\n')"
  )

  result <- npRmpi_run_rscript_subprocess(
    lines=lines, timeout=90L, env=env)
  info <- paste(result$output, collapse="\n")
  expect_equal(result$status, 0L, info=info)
  expect_true(any(grepl("NPRMPI_REGRESSION_EXTENDED_GNN_OK",
                        result$output, fixed=TRUE)), info=info)
})
