test_that("beta nearest-neighbor kernel sums honor local execution", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for beta NN subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(7321)",
      "train <- data.frame(x = runif(31), z = runif(31))",
      "eval <- data.frame(x = seq(0.08, 0.92, length.out = 7),",
      "                   z = seq(0.91, 0.09, length.out = 7))",
      "rhs <- matrix(rnorm(nrow(train)), ncol = 1L)",
      "for (bwtype in c('generalized_nn', 'adaptive_nn')) {",
      "  bw <- npudistbw(",
      "    dat = train, bws = c(5, 6), bwtype = bwtype,",
      "    ckertype = 'beta', ckerbound = 'range',",
      "    bandwidth.compute = FALSE",
      "  )",
      "  for (operator.name in c('normal', 'integral')) {",
      "    operator <- rep.int(operator.name, ncol(train))",
      "    got <- npRmpi:::.np_exact_operator_apply(",
      "      bw, train, eval, operator, rhs, 'beta NN local execution test'",
      "    )",
      "    direct <- npRmpi:::.np_direct_operator_apply(",
      "      bw, train, eval, operator, rhs, 'beta NN direct oracle'",
      "    )",
      "    stopifnot(isTRUE(all.equal(",
      "      unname(got), unname(direct), tolerance = 2e-13",
      "    )))",
      "  }",
      "}",
      "cat('BETA_NN_LOCAL_EXECUTION_OK\\n')"
    ),
    timeout = 30L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_true(res$status %in% c(0L, 137L), info = info)
  expect_true(any(grepl("BETA_NN_LOCAL_EXECUTION_OK", res$output,
                        fixed = TRUE)),
              info = info)
})
