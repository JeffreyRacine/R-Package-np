test_that("generalized NN regression bandwidth is stable on duplicated continuous data", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  x <- data.frame(x = rep(1:5, each = 4))
  y <- c(
    1, 1.1, 0.9, 1.05,
    2, 2.1, 1.9, 2.05,
    3, 3.1, 2.9, 3.05,
    4, 4.2, 3.8, 4.1,
    5, 5.2, 4.8, 5.1
  )

  bw <- npregbw(
    xdat = x,
    ydat = y,
    bwtype = "generalized_nn",
    nmulti = 1
  )

  expect_equal(as.numeric(bw$bw), 2)
  expect_equal(as.numeric(bw$fval), 0.1717472951192434, tolerance = 1e-14)
})

test_that("adaptive NN regression bandwidth is stable on duplicated continuous data", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  x <- data.frame(x = rep(1:5, each = 4))
  y <- c(
    1, 1.1, 0.9, 1.05,
    2, 2.1, 1.9, 2.05,
    3, 3.1, 2.9, 3.05,
    4, 4.2, 3.8, 4.1,
    5, 5.2, 4.8, 5.1
  )

  bw <- npregbw(
    xdat = x,
    ydat = y,
    bwtype = "adaptive_nn",
    nmulti = 1
  )

  expect_equal(as.numeric(bw$bw), 3)
  expect_equal(as.numeric(bw$fval), 0.1717472951192434, tolerance = 1e-14)
})
