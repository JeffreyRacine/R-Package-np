test_that("calibration names map to retained family owners", {
  for (family in c("density", "distribution", "regression")) {
    expect_identical(.np_plot_calibration_method("symmetric-sd", family), "absolute")
    expect_identical(.np_plot_calibration_method("asymmetric-sd", family),
      if (family == "density") "percentile-tails" else "separate-tails")
    expect_identical(.np_plot_calibration_method("tail-envelope", family), "tail-envelope")
  }
  expect_identical(np_boot_control()$calibration, "default")
  expect_identical(np_boot_control(calibration = "smoothed-percentile")$calibration,
                   "smoothed-percentile")
})

test_that("calibrated stream allocation is deterministic and preserves caller RNG", {
  set.seed(771932)
  entry <- .Random.seed
  a <- .np_udens_calibrated_streams(entry)
  expect_identical(a, .np_udens_calibrated_streams(entry))
  expect_identical(entry, .Random.seed)
  expect_length(a$outer, 39L)
  starts <- c(a$original, unlist(a$outer, recursive = FALSE))
  expect_identical(anyDuplicated(vapply(starts, paste, "", collapse = ",")), 0L)
})

# Full public 39-outer chains, MPI ownership and saved-reference equivalence
# are qualified in the scoped delivery harness, not repeated in every check.
