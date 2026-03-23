test_that("npregression object-fed data plot works after npRmpi.quit(force=TRUE)", {
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
      "fit <- npreg(txdat = x, tydat = y)",
      "plot.live <- plot(fit, plot.behavior = \"data\", perspective = FALSE, view = \"fixed\")",
      "npRmpi.quit(force = TRUE)",
      "plot.post <- plot(fit, plot.behavior = \"data\", perspective = FALSE, view = \"fixed\")",
      "stopifnot(identical(as.numeric(plot.live$r1$mean), as.numeric(plot.post$r1$mean)))",
      "stopifnot(identical(as.numeric(plot.live$r1$eval[[1L]]), as.numeric(plot.post$r1$eval[[1L]])))",
      "stopifnot(identical(as.numeric(plot.live$r1$bws$bw), as.numeric(plot.post$r1$bws$bw)))",
      "cat('post_quit_plot_ok\\n')"
    ),
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
  expect_true(any(grepl("post_quit_plot_ok", out$output, fixed = TRUE)))
})
