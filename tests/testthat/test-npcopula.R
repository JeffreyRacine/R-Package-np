test_that("npcopula basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  # npcopula can take a dbandwidth (from npudistbw) for copula 
  # or bandwidth (from npudensbw) for copula density
  
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Copula CDF
  cop <- npcopula(data=faithful, bws=bw)
  expect_s3_class(cop, "data.frame")
  expect_true("copula" %in% names(cop))
  expect_true(all(cop$copula >= 0 & cop$copula <= 1))
  
  # Copula density
  bw_dens <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  cop_dens <- npcopula(data=faithful, bws=bw_dens)
  expect_s3_class(cop_dens, "data.frame")
  expect_true("copula" %in% names(cop_dens))
  expect_true(all(cop_dens$copula >= 0))
})

test_that("npcopula u-grid path completes with active slaves", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(1)
  n <- 120L
  mydat <- data.frame(
    x = rnorm(n),
    y = 0.8 * rnorm(n) + rnorm(n, sd = 0.2)
  )
  u <- data.frame(x = c(0.25, 0.5, 0.75), y = c(0.25, 0.5, 0.75))

  bw <- npudistbw(dat = mydat, bws = c(0.8, 0.8), bandwidth.compute = FALSE)
  cop <- npcopula(bws = bw, data = mydat, u = u, n.quasi.inv = 60)

  expect_s3_class(cop, "data.frame")
  expect_true(all(c("copula", "u1", "u2", "x", "y") %in% names(cop)))
  expect_equal(nrow(cop), 9L)
  expect_true(all(is.finite(cop$copula)))
})
