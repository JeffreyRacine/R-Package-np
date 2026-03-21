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

test_that("npcdens cell cv.ml exhaustive search maximizes over degrees", {
  skip_if_not(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 50
  x <- runif(n)
  y <- rbeta(n, 1, 1)

  manual <- lapply(0:3, function(d) {
    npcdens(
      txdat = data.frame(x = x),
      tydat = data.frame(y = y),
      cxkerbound = "range",
      cykerbound = "range",
      regtype = "lp",
      degree = d,
      bwtype = "fixed",
      bwmethod = "cv.ml",
      nmulti = 1
    )$bws
  })
  manual.fval <- vapply(manual, function(bw) as.numeric(bw$fval[1L]), numeric(1L))
  expected.degree <- which.max(manual.fval) - 1L

  auto <- npcdens(
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    cxkerbound = "range",
    cykerbound = "range",
    regtype = "lp",
    degree.min = 0,
    degree.max = 3,
    degree.select = "exhaustive",
    search.engine = "cell",
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 1
  )

  expect_equal(auto$bws$degree, expected.degree)
  expect_equal(as.numeric(auto$bws$fval[1L]), max(manual.fval), tolerance = 1e-8)
})

test_that("npcdens nomad cv.ml retains the best objective found", {
  skip_if_not(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 60
  x <- runif(n)
  y <- rbeta(n, 1, 1)

  auto <- npcdens(
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    cxkerbound = "range",
    cykerbound = "range",
    regtype = "lp",
    degree.min = 0,
    degree.max = 3,
    degree.select = "coordinate",
    search.engine = "nomad+powell",
    degree.verify = FALSE,
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 2
  )

  expect_true(is.finite(auto$bws$degree.search$baseline.fval))
  expect_true(is.finite(auto$bws$degree.search$best.fval))
  expect_true(auto$bws$degree.search$best.fval >= auto$bws$degree.search$baseline.fval - 1e-8)
  expect_equal(as.numeric(auto$bws$fval[1L]), as.numeric(auto$bws$degree.search$best.fval), tolerance = 1e-8)
  expect_identical(auto$bws$degree.search$direction, "max")
  expect_identical(
    as.numeric(auto$bws$nomad.restart.fval),
    as.numeric(vapply(auto$bws$nomad.restart.results, `[[`, numeric(1L), "objective"))
  )
  expect_equal(as.integer(auto$bws$nomad.best.restart), as.integer(which.max(auto$bws$nomad.restart.fval)))
  expect_lte(
    auto$bws$nomad.restart.fval[auto$bws$nomad.best.restart],
    as.numeric(auto$bws$degree.search$best.fval) + 1e-6
  )
})
