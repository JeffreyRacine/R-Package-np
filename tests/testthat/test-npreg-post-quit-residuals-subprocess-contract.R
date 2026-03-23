test_that("npregression residuals work after npRmpi.quit(force=TRUE)", {
  skip_if(.mpi_check_context())

  env <- npRmpi_subprocess_env()
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "set.seed(20260323)",
      "x <- data.frame(x = seq(0, 1, length.out = 12))",
      "y <- sin(2 * pi * x$x)",
      "fit <- npreg(txdat = x, tydat = y, residuals = FALSE)",
      "res.live <- residuals(npreg(bws = fit$bws, residuals = TRUE))",
      "npRmpi.quit(force = TRUE)",
      "res.post <- residuals(fit)",
      "stopifnot(identical(as.numeric(res.live), as.numeric(res.post)))",
      "cat('post_quit_residuals_ok\\n')"
    ),
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
  expect_true(any(grepl("post_quit_residuals_ok", out$output, fixed = TRUE)))
})
