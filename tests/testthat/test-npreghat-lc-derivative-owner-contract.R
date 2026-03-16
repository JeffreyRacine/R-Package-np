npreghat_lc_derivative_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_LC_DERIV_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "set.seed(410)",
    "n <- 60L",
    "x <- rnorm(n)",
    "y <- x^2 + rnorm(n, sd = 0.2)",
    "train <- data.frame(x = x)",
    "eval <- data.frame(x = seq(min(x), max(x), length.out = 31L))",
    "run_case <- function(bt, bv, tol = 1e-10) {",
    "  bw <- npregbw(y ~ x, data = data.frame(train, y = y), regtype = 'lc', bwtype = bt, bandwidth.compute = FALSE, bws = bv)",
    "  Ht <- npreghat(bws = bw, txdat = train, s = 1L)",
    "  He <- npreghat(bws = bw, txdat = train, exdat = eval, s = 1L)",
    "  at <- npreghat(bws = bw, txdat = train, y = y, s = 1L, output = 'apply')",
    "  ae <- npreghat(bws = bw, txdat = train, exdat = eval, y = y, s = 1L, output = 'apply')",
    "  gt <- npreg(bws = bw, gradients = TRUE)",
    "  ge <- npreg(bws = bw, exdat = eval, gradients = TRUE)",
    "  stopifnot(isTRUE(all.equal(drop(Ht %*% y), at, tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(Ht %*% y), gt$grad[, 1L], tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(He %*% y), ae, tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(He %*% y), ge$grad[, 1L], tolerance = tol)))",
    "}",
    "run_case('fixed', 0.45)",
    "run_case('generalized_nn', 9)",
    "run_case('adaptive_nn', 9)",
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

test_that("npRmpi lc derivative owner matches npreg and apply on common cells", {
  npreghat_lc_derivative_owner_case()
})
