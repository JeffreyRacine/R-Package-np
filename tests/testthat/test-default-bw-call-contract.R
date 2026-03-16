test_that("default estimators build bandwidth calls in caller frame", {
  skip_if_not(spawn_mpi_slaves(1))
  on.exit(close_mpi_slaves(), add = TRUE)

  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  run_dens <- function() {
    set.seed(101)
    x <- runif(30)
    npudens(~x, data = data.frame(x = x), bwmethod = "normal-reference")
  }

  run_dist <- function() {
    set.seed(202)
    x <- runif(35)
    npudist(~x, data = data.frame(x = x), bwmethod = "normal-reference")
  }

  run_reg <- function() {
    set.seed(303)
    x <- runif(45)
    y <- sin(2 * pi * x) + rnorm(45, sd = 0.2)
    npreg(txdat = x, tydat = y, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  }

  fit_dens <- run_dens()
  fit_dist <- run_dist()
  fit_reg <- run_reg()

  expect_s3_class(fit_dens, "npdensity")
  expect_s3_class(fit_dist, "npdistribution")
  expect_s3_class(fit_reg, "npregression")

  expect_true(any(grepl("bandwidth$", class(fit_dens$bws))))
  expect_true(any(grepl("bandwidth$", class(fit_dist$bws))))
  expect_true(any(grepl("bandwidth$", class(fit_reg$bws))))
})

test_that("npreg default bandwidth selection remains finite after bounded conditional selectors", {
  skip_if_not(spawn_mpi_slaves(1))
  on.exit(close_mpi_slaves(), add = TRUE)

  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260316)
  n <- 48
  x <- runif(n)
  y <- runif(n)

  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerorder = 2L,
    cykerorder = 2L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )

  expect_true(all(is.finite(as.numeric(bw.cd$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd$ybw))))

  fit.reg <- npreg(
    txdat = x,
    tydat = sin(2 * pi * x) + rnorm(n, sd = 0.15),
    regtype = "lc",
    bwmethod = "cv.ls",
    nmulti = 1
  )

  expect_s3_class(fit.reg, "npregression")
  expect_true(any(grepl("bandwidth$", class(fit.reg$bws))))
  expect_equal(length(fit.reg$mean), NROW(fit.reg$eval))
})
