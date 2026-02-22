test_that("extended bw generics route named data args without bws to NULL methods", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260222)
  n <- 35
  x <- rnorm(n)
  z <- rnorm(n)
  y <- x + z + rnorm(n)

  x1 <- data.frame(x = x)
  x2 <- data.frame(x = x, z = z)
  y1 <- data.frame(y = y)
  z1 <- data.frame(z = z)

  bw_cdens <- npcdensbw(xdat = x1, ydat = y1, bwmethod = "normal-reference")
  bw_cdist <- npcdistbw(xdat = x1, ydat = y1, bwmethod = "normal-reference")
  bw_pl <- npplregbw(xdat = x1, ydat = y, zdat = z1, nmulti = 1)
  bw_si <- npindexbw(xdat = x2, ydat = y, nmulti = 1)
  bw_sc <- npscoefbw(xdat = x1, ydat = y, zdat = z1, nmulti = 1)

  expect_s3_class(bw_cdens, "conbandwidth")
  expect_s3_class(bw_cdist, "condbandwidth")
  expect_s3_class(bw_pl, "plbandwidth")
  expect_s3_class(bw_si, "sibandwidth")
  expect_s3_class(bw_sc, "scbandwidth")
})
