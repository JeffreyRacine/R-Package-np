test_that("named bws formula dispatch matches positional density/distribution routes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260323)
  d <- data.frame(
    x = seq(0.1, 0.9, length.out = 9),
    y = seq(0.2, 1.0, length.out = 9)
  )
  nd.u <- data.frame(x = c(0.2, 0.5, 0.8))
  nd.c <- data.frame(x = c(0.2, 0.5, 0.8), y = c(0.25, 0.55, 0.85))

  bw.ud <- npudensbw(~ x, data = d, newdata = nd.u)
  ud.pos <- npudens(bws = bw.ud, newdata = nd.u)
  ud.named <- npudens(bws = ~ x, data = d, newdata = nd.u)
  expect_equal(as.numeric(ud.named$dens), as.numeric(ud.pos$dens), tolerance = 0)
  expect_equal(as.numeric(ud.named$derr), as.numeric(ud.pos$derr), tolerance = 0)

  bw.udist <- npudistbw(~ x, data = d, newdata = nd.u)
  udist.pos <- npudist(bws = bw.udist, newdata = nd.u)
  udist.named <- npudist(bws = ~ x, data = d, newdata = nd.u)
  expect_equal(as.numeric(udist.named$dist), as.numeric(udist.pos$dist), tolerance = 0)
  expect_equal(as.numeric(udist.named$derr), as.numeric(udist.pos$derr), tolerance = 0)

  bw.cd <- npcdensbw(y ~ x, data = d, newdata = nd.c)
  cd.pos <- npcdens(bws = bw.cd, data = d, newdata = nd.c)
  cd.named <- npcdens(bws = y ~ x, data = d, newdata = nd.c)
  expect_lt(max(abs(as.numeric(cd.named$condens) - as.numeric(cd.pos$condens))), 1e-2)
  expect_lt(max(abs(as.numeric(cd.named$conderr) - as.numeric(cd.pos$conderr))), 1e-2)

  bw.cdist <- npcdistbw(y ~ x, data = d, newdata = nd.c)
  cdist.pos <- npcdist(bws = bw.cdist, data = d, newdata = nd.c)
  cdist.named <- npcdist(bws = y ~ x, data = d, newdata = nd.c)
  expect_lt(max(abs(as.numeric(cdist.named$condist) - as.numeric(cdist.pos$condist))), 1e-2)
  expect_lt(max(abs(as.numeric(cdist.named$conderr) - as.numeric(cdist.pos$conderr))), 1e-2)
})
