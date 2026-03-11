test_that("plot contract: asymptotic consumer payloads are supported in npRmpi session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(922)
  n <- 55
  x1 <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y.cont <- sin(2 * pi * x1) + 0.35 * x2 + rnorm(n, sd = 0.08)
  y.bin <- as.numeric(x1 + x2 + rnorm(n, sd = 0.2) > 1)

  sbw <- npindexbw(
    xdat = data.frame(x1 = x1, x2 = x2),
    ydat = y.cont,
    bws = c(1, 6, 6),
    bandwidth.compute = FALSE,
    method = "ichimura",
    bwtype = "adaptive_nn",
    regtype = "ll"
  )
  sbw.bin <- npindexbw(
    xdat = data.frame(x1 = x1, x2 = x2),
    ydat = y.bin,
    bws = c(1, 6, 6),
    bandwidth.compute = FALSE,
    method = "kleinspady",
    bwtype = "adaptive_nn",
    regtype = "ll"
  )
  pbw <- npplregbw(
    xdat = data.frame(x = x1),
    zdat = data.frame(z = z),
    ydat = y.cont,
    bws = matrix(c(6, 6), nrow = 2L),
    bandwidth.compute = FALSE,
    bwmethod = "cv.aic",
    bwtype = "adaptive_nn",
    regtype = "ll"
  )

  out.si <- suppressWarnings(plot(
    sbw,
    xdat = data.frame(x1 = x1, x2 = x2),
    ydat = y.cont,
    plot.behavior = "data",
    plot.errors.method = "asymptotic",
    plot.errors.type = "pointwise",
    perspective = FALSE
  ))[[1]]
  out.si.bin <- suppressWarnings(plot(
    sbw.bin,
    xdat = data.frame(x1 = x1, x2 = x2),
    ydat = y.bin,
    plot.behavior = "data",
    plot.errors.method = "asymptotic",
    plot.errors.type = "pointwise",
    perspective = FALSE
  ))[[1]]
  out.pl <- suppressWarnings(plot(
    pbw,
    xdat = data.frame(x = x1),
    ydat = y.cont,
    zdat = data.frame(z = z),
    plot.behavior = "data",
    plot.errors.method = "asymptotic",
    plot.errors.type = "pointwise",
    perspective = FALSE
  ))[[1]]

  expect_s3_class(out.si, "singleindex")
  expect_s3_class(out.si.bin, "singleindex")
  expect_s3_class(out.pl, "plregression")
  expect_true(all(is.finite(out.si$mean)))
  expect_true(all(is.finite(out.si$merr)))
  expect_true(all(is.finite(out.si.bin$mean)))
  expect_true(all(is.finite(out.si.bin$merr)))
  expect_true(all(is.finite(out.pl$mean)))
  expect_true(all(is.finite(out.pl$merr)))
})
