test_that("default estimators build bandwidth calls in caller frame", {
  skip_if_not(spawn_mpi_slaves(1))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old.opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old.opts), add = TRUE)

  run_dens <- function() {
    set.seed(101)
    x <- runif(30)
    npudens(tdat = x, bwmethod = "normal-reference")
  }

  run_dist <- function() {
    set.seed(202)
    x <- runif(35)
    npudist(tdat = x, bwmethod = "normal-reference")
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
