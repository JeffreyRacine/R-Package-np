test_that("npconmode proper helper enforces binary complement probabilities", {
  helper <- getFromNamespace(".npConmodeProperProbabilities", "npRmpi")
  raw <- matrix(c(-0.2, 1.2,
                  0.4, 0.8,
                  0.3, 0.3), ncol = 2L, byrow = TRUE)

  out <- helper(raw, levels = c("0", "1"), proper = TRUE)

  expect_true(isTRUE(out$proper.requested))
  expect_true(isTRUE(out$proper.applied))
  expect_equal(rowSums(out$probabilities), rep(1, nrow(raw)), tolerance = 1e-12)
  expect_true(all(out$probabilities >= -1e-12))
  expect_equal(out$probabilities[, 2L], 1 - out$probabilities[, 1L], tolerance = 1e-12)
  expect_identical(out$proper.info$reason, "projected")
})

test_that("npconmode proper defaults follow the canonical regression type", {
  effective <- getFromNamespace(".npConmodeEffectiveProper", "npRmpi")

  expect_false(effective(list(regtype = "lc"), NULL))
  expect_true(effective(list(regtype = "ll"), NULL))
  expect_true(effective(list(regtype = "lp"), NULL))
  expect_false(effective(list(regtype = "lp", regtype.engine = "lc"), NULL))
  expect_true(effective(list(regtype = "lc", regtype.engine = "lp"), NULL))
  expect_false(effective(list(regtype = "lp"), FALSE))
  expect_true(effective(list(regtype = "lc"), TRUE))
})

test_that("npconmode formula route forwards NOMAD shortcut controls", {
  skip_if_not_installed("crs")
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260511)
  n <- 30L
  d <- data.frame(
    x = seq(-1, 1, length.out = n)
  )
  d$y <- factor(rbinom(n, 1L, plogis(1.5 * d$x)))

  capture.output(
    fit <- npconmode(
      y ~ x,
      data = d,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L
    )
  )

  expect_s3_class(fit, "conmode")
  expect_true(isTRUE(fit$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit$bws$regtype.engine), "lp")
})
