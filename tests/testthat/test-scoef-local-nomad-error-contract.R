test_that("local fixed smooth-coefficient search preserves unexpected preparation errors", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  withr::local_options(list(np.messages = FALSE))
  i <- seq_len(48L)
  d <- data.frame(x = sin(i*sqrt(2)), z = cos(i*sqrt(3)))
  d$y <- d$x*(1+d$z) + sin(i*sqrt(5))/5
  original <- structure(list(message = "unexpected smooth-coefficient preparation error",
    call = quote(moment_preparation()), token = 643L),
    class = c("scoef_preparation_test_error", "error", "condition"))
  real.moment <- .npscoefbw_nomad_moment_state
  for (visit in c(1L, 2L)) {
    calls <- 0L
    got <- testthat::with_mocked_bindings({
      tryCatch(.npRmpi_with_local_regression(
        npscoefbw(y ~ x | z, data = d, nomad = TRUE,
          degree.min = 0L, degree.max = 1L, nmulti = 1L, random.seed = 42L)),
        error = identity)
    }, .npscoefbw_nomad_moment_state = function(...) {
      calls <<- calls + 1L
      if (calls == visit) stop(original)
      real.moment(...)
    }, .package = "npRmpi")
    expect_identical(got, original)
    expect_identical(calls, visit)
  }
})
