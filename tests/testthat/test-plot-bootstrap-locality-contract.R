library(npRmpi)

test_that("session-route regression bootstrap plot stays off local bootstrap wrapper", {
  npRmpi.init(nslaves = 1, quiet = TRUE)
  on.exit(npRmpi.quit(), add = TRUE)

  set.seed(9420)
  xdat <- data.frame(x = rnorm(20))
  ydat <- xdat$x + rnorm(20)

  bw <- npregbw(xdat = xdat, ydat = ydat, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  fit <- npreg(txdat = xdat, tydat = ydat, bws = bw)

  ctr <- new.env(parent = emptyenv())
  ctr$n <- 0L
  trace(
    ".npRmpi_with_local_bootstrap",
    where = asNamespace("npRmpi"),
    tracer = bquote(assign("n", get("n", envir = .(ctr)) + 1L, envir = .(ctr))),
    print = FALSE
  )
  on.exit(try(untrace(".npRmpi_with_local_bootstrap", where = asNamespace("npRmpi")), silent = TRUE), add = TRUE)

  out <- plot(
    fit,
    xdat = xdat,
    ydat = ydat,
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.num = 5
  )

  expect_type(out, "list")
  expect_identical(ctr$n, 0L)
})

test_that("regression plot direct eval stays on local compiled owner", {
  plot_eval <- getFromNamespace(".np_plot_regression_eval", "npRmpi")
  fn.body <- paste(deparse(body(plot_eval), width.cutoff = 500L), collapse = " ")

  expect_match(fn.body, "\\.np_plot_with_local_compiled_eval\\(")
})

test_that("wild regression bootstrap stays on local regression wrappers", {
  boot_fun <- getFromNamespace("compute.bootstrap.errors.rbandwidth", "npRmpi")
  fn.body <- paste(deparse(body(boot_fun), width.cutoff = 500L), collapse = " ")

  expect_match(fn.body, "fit\\.mean\\.train")
  expect_match(fn.body, "\\.npRmpi_with_local_regression\\(suppressWarnings\\(npreg\\.rbandwidth\\(")
  expect_match(fn.body, "\\.npRmpi_with_local_regression\\(suppressWarnings\\(npreghat\\.rbandwidth\\(")
})
