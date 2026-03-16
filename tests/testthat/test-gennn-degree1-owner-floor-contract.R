gennn_degree1_owner_floor_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "GENNN_DEGREE1_OWNER_FLOOR_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "set.seed(20260315)",
    "n <- 60L",
    "x <- rnorm(n)",
    "y <- x^2 + rnorm(n, sd = 0.2)",
    "tx <- data.frame(x = x)",
    "ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))",
    "make_reg_bw <- function(regtype) {",
    "  args <- list(xdat = tx, ydat = y, regtype = regtype, bwtype = 'generalized_nn', bws = 2L, bandwidth.compute = FALSE)",
    "  if (identical(regtype, 'lp')) {",
    "    args$basis <- 'glp'",
    "    args$degree <- 1L",
    "    args$bernstein.basis <- FALSE",
    "  }",
    "  do.call(npregbw, args)",
    "}",
    "bw.ll <- make_reg_bw('ll')",
    "bw.lp <- make_reg_bw('lp')",
    "fit.ll <- npreg(bws = bw.ll, txdat = tx, tydat = y, exdat = ex)",
    "fit.lp <- npreg(bws = bw.lp, txdat = tx, tydat = y, exdat = ex)",
    "H.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, output = 'matrix')",
    "H.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex, output = 'matrix')",
    "hat.apply.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, y = y, output = 'apply')",
    "hat.apply.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex, y = y, output = 'apply')",
    "stopifnot(isTRUE(all.equal(as.vector(fit.ll$mean), as.vector(fit.lp$mean), tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(as.vector(hat.apply.ll), as.vector(hat.apply.lp), tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.ll %*% y), as.vector(hat.apply.ll), tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.lp %*% y), as.vector(hat.apply.lp), tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.ll %*% y), as.vector(fit.ll$mean), tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.lp %*% y), as.vector(fit.lp$mean), tolerance = 1e-10)))",
    "n <- 80L",
    "z <- runif(n)",
    "x <- 0.7 * z + rnorm(n, sd = 0.2)",
    "y <- 1 + 1.1 * x + sin(2 * pi * z) + rnorm(n, sd = 0.05)",
    "tx <- data.frame(x = x)",
    "tz <- data.frame(z = z)",
    "ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))",
    "ez <- data.frame(z = seq(min(z), max(z), length.out = 31L))",
    "bws.fixed <- matrix(c(2L, 9L), nrow = 2L)",
    "make_pl_bw <- function(regtype) {",
    "  args <- list(xdat = tx, zdat = tz, ydat = y, regtype = regtype, bwtype = 'generalized_nn', bws = bws.fixed, bandwidth.compute = FALSE)",
    "  if (identical(regtype, 'lp')) {",
    "    args$basis <- 'glp'",
    "    args$degree <- 1L",
    "    args$bernstein.basis <- FALSE",
    "  }",
    "  do.call(npplregbw, args)",
    "}",
    "bw.ll <- make_pl_bw('ll')",
    "bw.lp <- make_pl_bw('lp')",
    "fit.ll <- npplreg(bws = bw.ll, txdat = tx, tzdat = tz, tydat = y, exdat = ex, ezdat = ez)",
    "fit.lp <- npplreg(bws = bw.lp, txdat = tx, tzdat = tz, tydat = y, exdat = ex, ezdat = ez)",
    "H.ll <- npplreghat(bws = bw.ll, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, output = 'matrix')",
    "H.lp <- npplreghat(bws = bw.lp, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, output = 'matrix')",
    "hat.apply.ll <- npplreghat(bws = bw.ll, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, y = y, output = 'apply')",
    "hat.apply.lp <- npplreghat(bws = bw.lp, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, y = y, output = 'apply')",
    "stopifnot(isTRUE(all.equal(as.vector(fit.ll$mean), as.vector(fit.lp$mean), tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(as.vector(hat.apply.ll), as.vector(hat.apply.lp), tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.ll %*% y), as.vector(hat.apply.ll), tolerance = 1e-7)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.lp %*% y), as.vector(hat.apply.lp), tolerance = 1e-7)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.ll %*% y), as.vector(fit.ll$mean), tolerance = 1e-7)))",
    "stopifnot(isTRUE(all.equal(as.vector(H.lp %*% y), as.vector(fit.lp$mean), tolerance = 1e-7)))",
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

test_that("npRmpi generalized-nn degree-1 owner routes stay aligned at NN floor", {
  gennn_degree1_owner_floor_case()
})
