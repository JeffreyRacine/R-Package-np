test_that("saved calibration controls name a family-specific migration", {
  advise <- function(family, value) {
    old <- np_boot_control()
    old$calibration <- value
    .np_plot_family_interval_request(list(boot.control = old), family)
  }
  expect_error(advise("npudens", "separate-tails"), "use interval='separate-tails'", fixed = TRUE)
  expect_error(advise("npudist", "separate-tails"), "use interval='asymmetric-sd'", fixed = TRUE)
  for (family in c("npplreg", "npindex", "npscoef", "npqreg", "nplsqreg", "npconmode", "npcopula"))
    expect_error(advise(family, "symmetric-sd"), "band='simultaneous', interval='symmetric-sd'", fixed = TRUE)
  expect_error(advise("npudist", "legacy"), "np_udist_bias_control", fixed = TRUE)
  expect_error(advise("npudens", "basic-envelope"), "unavailable for npudens", fixed = TRUE)
  expect_error(advise("npudist", "percentile-tails"), "unavailable for npudist", fixed = TRUE)
})
