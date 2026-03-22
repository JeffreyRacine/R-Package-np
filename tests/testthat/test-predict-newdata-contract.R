test_that("predict aliases newdata to native eval args for default npreg/npudens/npudist/npindex", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(70)
  y <- rnorm(70)
  nd <- data.frame(x = c(0.1, 0.3, 0.7))

  bw.reg <- npregbw(
    xdat = data.frame(x = x),
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.reg <- npreg(bws = bw.reg, txdat = data.frame(x = x), tydat = y)
  expect_equal(
    as.numeric(predict(fit.reg, newdata = nd)),
    as.numeric(predict(fit.reg, exdat = nd)),
    tolerance = 1e-12
  )

  bw.den <- npudensbw(
    dat = data.frame(x = x),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.den <- npudens(bws = bw.den, tdat = data.frame(x = x))
  expect_equal(
    as.numeric(predict(fit.den, newdata = nd)),
    as.numeric(predict(fit.den, edat = nd)),
    tolerance = 1e-12
  )

  bw.dist <- npudistbw(
    dat = data.frame(x = x),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.dist <- npudist(bws = bw.dist, tdat = data.frame(x = x))
  expect_equal(
    as.numeric(predict(fit.dist, newdata = nd)),
    as.numeric(predict(fit.dist, edat = nd)),
    tolerance = 1e-12
  )

  x2 <- runif(70)
  nd.si <- data.frame(x = c(0.15, 0.35), x2 = c(0.4, 0.8))
  bw.si <- npindexbw(
    xdat = data.frame(x = x, x2 = x2),
    ydat = y,
    bws = c(0.25, 0.25, 1),
    bandwidth.compute = FALSE
  )
  fit.si <- npindex(
    bws = bw.si,
    txdat = data.frame(x = x, x2 = x2),
    tydat = y
  )
  expect_equal(
    as.numeric(predict(fit.si, newdata = nd.si)),
    as.numeric(predict(fit.si, exdat = nd.si)),
    tolerance = 1e-12
  )
})

test_that("predict aliases newdata to exdat/eydat for default npcdens/npcdist", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(60)
  y <- runif(60)
  nd <- data.frame(y = c(0.2, 0.5), x = c(0.1, 0.8))

  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.2),
    bandwidth.compute = FALSE
  )
  fit.cd <- npcdens(
    bws = bw.cd,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y)
  )
  expect_equal(
    as.numeric(predict(fit.cd, newdata = nd)),
    as.numeric(predict(fit.cd, exdat = nd["x"], eydat = nd["y"])),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit.cd, newdata = data.frame(x = c(0.1, 0.2))),
    "must include columns"
  )

  bw.cdist <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.2),
    bandwidth.compute = FALSE
  )
  fit.cdist <- npcdist(
    bws = bw.cdist,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y)
  )
  expect_equal(
    as.numeric(predict(fit.cdist, newdata = nd)),
    as.numeric(predict(fit.cdist, exdat = nd["x"], eydat = nd["y"])),
    tolerance = 1e-12
  )
})

test_that("predict aliases newdata to the explicit-evaluation slice route for npcdens/npcdist", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260322)
  x <- runif(70, -1, 1)
  y <- sin(2 * pi * x) + rnorm(70, sd = 0.2)
  nd <- rbind(
    data.frame(y = c(-0.7, -0.15, 0.45), x = rep(-0.35, 3L)),
    data.frame(y = c(-0.25, 0.2, 0.75, 1.0), x = rep(0.4, 4L))
  )
  ctrl <- list(mode = "slice", slice.grid.size = 21L, slice.extend.factor = 0)

  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.27, 0.21),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit.cd <- npcdens(
    bws = bw.cd,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )
  expect_equal(
    as.numeric(predict(fit.cd, newdata = nd, proper.control = ctrl)),
    as.numeric(predict(
      fit.cd,
      exdat = nd["x"],
      eydat = nd["y"],
      proper.control = ctrl
    )),
    tolerance = 1e-10
  )

  bw.cdist <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.27, 0.21),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit.cdist <- npcdist(
    bws = bw.cdist,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )
  expect_equal(
    as.numeric(predict(fit.cdist, newdata = nd, proper.control = ctrl)),
    as.numeric(predict(
      fit.cdist,
      exdat = nd["x"],
      eydat = nd["y"],
      proper.control = ctrl
    )),
    tolerance = 1e-10
  )
})

test_that("predict aliases newdata to exdat/ezdat for default npscoef/npplreg", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260227)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 1.5 * z + rnorm(n, sd = 0.05)
  nd <- data.frame(x = c(0.2, 0.5, 0.8), z = c(0.1, 0.4, 0.9))

  bw.sc <- npscoefbw(
    xdat = data.frame(x = x),
    ydat = y,
    zdat = data.frame(z = z),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.sc <- npscoef(
    bws = bw.sc,
    txdat = data.frame(x = x),
    tydat = y,
    tzdat = data.frame(z = z)
  )
  expect_equal(
    as.numeric(predict(fit.sc, newdata = nd)),
    as.numeric(predict(fit.sc, exdat = nd["x"], ezdat = nd["z"])),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit.sc, newdata = data.frame(x = nd$x)),
    "must include columns"
  )

  bw.pl <- npplregbw(
    xdat = data.frame(z = z),
    ydat = y,
    zdat = data.frame(x = x),
    bws = matrix(c(0.25, 0.25), nrow = 2),
    bandwidth.compute = FALSE
  )
  fit.pl <- npplreg(
    bws = bw.pl,
    txdat = data.frame(z = z),
    tydat = y,
    tzdat = data.frame(x = x)
  )
  expect_equal(
    as.numeric(predict(fit.pl, newdata = nd)),
    as.numeric(predict(fit.pl, exdat = nd["z"], ezdat = nd["x"])),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit.pl, newdata = data.frame(x = nd$x)),
    "must include columns"
  )
})
