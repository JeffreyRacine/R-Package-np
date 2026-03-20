test_that("npcdens stores continuous scale factors distinctly from bandwidths", {
  skip_if_not(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 80
  x <- runif(n)
  y <- rbeta(n, 1, 1)

  fit <- npcdens(
    y ~ x,
    cxkerbound = "range",
    cykerbound = "range",
    regtype = "lp",
    degree = 3,
    bwmethod = "cv.ml",
    nmulti = 1
  )

  bw <- fit$bws
  expect_false(isTRUE(all.equal(unname(bw$bandwidth$x), unname(bw$sfactor$x))))
  expect_false(isTRUE(all.equal(unname(bw$bandwidth$y), unname(bw$sfactor$y))))

  expected.x <- unname(bw$xbw / (as.numeric(bw$sdev[seq_len(bw$xncon)]) * bw$nconfac))
  expected.y <- unname(bw$ybw / (as.numeric(bw$sdev[bw$xncon + seq_len(bw$yncon)]) * bw$nconfac))

  expect_equal(unname(bw$sfactor$x), expected.x, tolerance = 1e-10)
  expect_equal(unname(bw$sfactor$y), expected.y, tolerance = 1e-10)
  expect_equal(unname(bw$sumNum$x), expected.x, tolerance = 1e-10)
  expect_equal(unname(bw$sumNum$y), expected.y, tolerance = 1e-10)
})
