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
    output = "data",
    errors = "asymptotic",
    band = "pointwise",
    perspective = FALSE
  ))
  out.si.bin <- suppressWarnings(plot(
    sbw.bin,
    xdat = data.frame(x1 = x1, x2 = x2),
    ydat = y.bin,
    output = "data",
    errors = "asymptotic",
    band = "pointwise",
    perspective = FALSE
  ))
  out.pl <- suppressWarnings(plot(
    pbw,
    xdat = data.frame(x = x1),
    ydat = y.cont,
    zdat = data.frame(z = z),
    output = "data",
    errors = "asymptotic",
    band = "pointwise",
    perspective = FALSE
  ))

  expect_s3_class(out.si, "np_plot_data")
  expect_s3_class(out.si.bin, "np_plot_data")
  expect_s3_class(out.pl, "np_plot_data")
  expect_true(all(is.finite(plot_target(out.si)$estimate)))
  expect_true(all(is.finite(plot_band(out.si, "pointwise"))))
  expect_true(all(is.finite(plot_target(out.si.bin)$estimate)))
  expect_true(all(is.finite(plot_band(out.si.bin, "pointwise"))))
  expect_true(all(is.finite(plot_target(out.pl)$estimate)))
  expect_true(all(is.finite(plot_band(out.pl, "pointwise"))))
})
