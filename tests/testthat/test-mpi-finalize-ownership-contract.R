test_that("compute sources no longer own MPI finalization", {
  files <- c(
    testthat::test_path("..", "..", "src", "jksum.c"),
    testthat::test_path("..", "..", "src", "kernele.c"),
    testthat::test_path("..", "..", "src", "statmods.c")
  )
  skip_if_not(all(file.exists(files)), "source C files unavailable in installed test context")

  for (path in files) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_false(grepl("MPI_Finalize\\s*\\(", text),
                 info = paste("unexpected compute-layer MPI_Finalize in", path))
  }
})

test_that("session route remains reusable after constant-variable bandwidth failure", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(123)",
      "xbad <- data.frame(x = rep(1, 10))",
      "ybad <- rnorm(10)",
      "msg <- tryCatch({ npregbw(xdat = xbad, ydat = ybad, nmulti = 1L); '' }, error = conditionMessage)",
      "stopifnot(nzchar(msg))",
      "x <- data.frame(x = seq(0, 1, length.out = 20))",
      "y <- sin(2*pi*x$x)",
      "bw <- npregbw(xdat = x, ydat = y, bws = 0.25, bandwidth.compute = FALSE)",
      "fit <- npreg(bws = bw, txdat = x, tydat = y)",
      "stopifnot(inherits(fit, 'npregression'))",
      "cat('MPI_FINALIZE_OWNERSHIP_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("MPI_FINALIZE_OWNERSHIP_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
