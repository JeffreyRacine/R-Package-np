test_that("bw generics route named data args without bws to NULL methods", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260222)
  x <- rnorm(30)
  y <- x + rnorm(30)
  dat <- data.frame(x = x, y = y)

  bw_reg <- npregbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.aic", nmulti = 1)
  bw_dens <- npudensbw(dat = dat, bwmethod = "normal-reference")
  bw_dist <- npudistbw(dat = dat, bwmethod = "normal-reference")

  expect_s3_class(bw_reg, "rbandwidth")
  expect_s3_class(bw_dens, "bandwidth")
  expect_s3_class(bw_dist, "dbandwidth")
})

test_that("bw object dispatch remains intact", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260222)
  x <- rnorm(25)
  y <- x + rnorm(25)
  dat <- data.frame(x = x, y = y)

  bw_reg0 <- npregbw(xdat = x, ydat = y, bws = 0.5, bandwidth.compute = FALSE, regtype = "lc")
  bw_reg1 <- npregbw(bws = bw_reg0, xdat = x, ydat = y, bandwidth.compute = FALSE, regtype = "lc")
  expect_equal(as.numeric(bw_reg1$bw), as.numeric(bw_reg0$bw))

  bw_dens0 <- npudensbw(dat = dat, bws = c(0.4, 0.6), bandwidth.compute = FALSE)
  bw_dens1 <- npudensbw(bws = bw_dens0, dat = dat, bandwidth.compute = FALSE)
  expect_equal(as.numeric(bw_dens1$bw), as.numeric(bw_dens0$bw))

  bw_dist0 <- npudistbw(dat = dat, bws = c(0.4, 0.6), bandwidth.compute = FALSE)
  bw_dist1 <- npudistbw(bws = bw_dist0, dat = dat, bandwidth.compute = FALSE)
  expect_equal(as.numeric(bw_dist1$bw), as.numeric(bw_dist0$bw))
})
