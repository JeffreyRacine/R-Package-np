test_that("plot works with autodispatch for non-bootstrap paths", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(21)
  n <- 60
  d <- data.frame(x = rnorm(n), y = rnorm(n))

  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  out <- suppressWarnings(plot(bw,
                               persp = FALSE,
                               view = "fixed",
                               plot.behavior = "data",
                               plot.errors.method = "none"))

  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("plot bootstrap path works under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(22)
  n <- 60
  d <- data.frame(x = rnorm(n), y = rnorm(n))

  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)

  out <- suppressWarnings(
    plot(bw,
         persp = FALSE,
         view = "fixed",
         plot.behavior = "data",
         plot.errors.method = "bootstrap",
         plot.errors.boot.num = 9)
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("autodispatch keeps formula bws usable for condensity plot()", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(42)
  n <- 120
  x <- rnorm(n)
  y <- rnorm(n)

  bw <- npcdensbw(y ~ x, nmulti = 1)
  fit <- npcdens(bws = bw)

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(fit, perspective = FALSE, plot.errors.method = "none"), NA)
  expect_false(grepl("\\.__npRmpi_autod_", paste(deparse(fit$bws$call), collapse = " ")))
})

test_that("autodispatch keeps npreg formula fits plotable", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(44)
  n <- 80
  d <- data.frame(
    x1 = runif(n),
    x2 = runif(n)
  )
  d$y <- sin(pi * (d$x1 + d$x2))^4 * sin(pi * d$x1)^2 + rnorm(n, sd = 0.1)

  fit <- npreg(y ~ x1 + x2, data = d, nomad = TRUE, nmulti = 1)

  expect_false(is.null(fit$bws$formula))
  expect_identical(as.character(fit$bws$call[[1L]]), "npregbw")

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(fit), NA)
})

test_that("npRmpi S3 methods stay registered after np namespace loads", {
  skip_if_not_installed("np")
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(45)
  n <- 80
  x <- runif(n, -1, 1)
  y <- rbeta(n, 1, 1)

  fit <- npreg(y ~ x, nomad = TRUE, nmulti = 1)

  expect_true("np" %in% loadedNamespaces())
  for (spec in list(
    c("plot", "npregression"),
    c("summary", "npregression"),
    c("summary", "rbandwidth"),
    c("print", "npregression"),
    c("predict", "npregression")
  )) {
    meth <- getS3method(spec[1L], spec[2L])
    expect_identical(environmentName(environment(meth)), "npRmpi")
  }

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(fit), NA)
})

test_that("autodispatch default condensity plot options stay scalar-safe", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(43)
  n <- 100
  x <- rnorm(n)
  y <- rnorm(n)

  bw <- npcdensbw(y ~ x, bws = c(1.0, 1.0), bandwidth.compute = FALSE)
  fit <- npcdens(bws = bw)

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(fit, perspective = FALSE), NA)
})

test_that("plot engine defaults: regression-class uses wild, unsupervised uses inid", {
  reg.engines <- c(
    ".np_plot_rbandwidth_engine",
    ".np_plot_scbandwidth_engine",
    ".np_plot_plbandwidth_engine",
    ".np_plot_sibandwidth_engine"
  )
  unsup.engines <- c(
    ".np_plot_bandwidth_engine",
    ".np_plot_dbandwidth_engine",
    ".np_plot_conbandwidth_engine",
    ".np_plot_condbandwidth_engine"
  )

  for (nm in reg.engines) {
    fn <- getFromNamespace(nm, "npRmpi")
    defaults <- eval(formals(fn)$plot.errors.boot.method)
    expect_identical(defaults[1L], "wild")
  }

  for (nm in unsup.engines) {
    fn <- getFromNamespace(nm, "npRmpi")
    defaults <- eval(formals(fn)$plot.errors.boot.method)
    expect_identical(defaults[1L], "inid")
  }
})
