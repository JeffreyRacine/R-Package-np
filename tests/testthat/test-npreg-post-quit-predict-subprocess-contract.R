test_that("npregression predict works after npRmpi.quit(force=TRUE)", {
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
      "nd <- data.frame(x = seq(0.1, 0.9, length.out = 5))",
      "fit <- npreg(txdat = x, tydat = y)",
      "pred.live <- predict(fit, newdata = nd)",
      "npRmpi.quit(force = TRUE)",
      "pred.post <- predict(fit, newdata = nd)",
      "stopifnot(identical(as.numeric(pred.live), as.numeric(pred.post)))",
      "cat('post_quit_predict_ok\\n')"
    ),
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
  expect_true(any(grepl("post_quit_predict_ok", out$output, fixed = TRUE)))
})
