test_that("response reconstruction preserves numeric values and rejects categorical provenance", {
  value <- list(mean = c(1, 2), resid = c(0.2, -0.1))
  expected <- value$mean + value$resid
  for (type in list(NULL, list(icon = TRUE, iuno = FALSE, iord = FALSE)))
    expect_identical(.np_plot_reconstruct_response(value, list(ydati = type), "plot"), expected)
  for (type in list(list(iuno = TRUE, iord = FALSE), list(iuno = FALSE, iord = TRUE)))
    expect_error(.np_plot_reconstruct_response(value, list(ydati = type), "plot"),
                 "must be a numeric ", fixed = TRUE)
  expect_identical(value, list(mean = c(1, 2), resid = c(0.2, -0.1)))
})

test_that("every fitted regression plot checks provenance before rebuilding outcomes", {
  # Exercise the actual plot wrappers with an otherwise valid numeric fit and
  # retained categorical metadata, the shape of the pre-Numeric-Y failure.
  d <- data.frame(x = seq(-1, 1, length.out = 28L))
  d$y <- sin(d$x) + cos(7 * d$x) * 0.1
  d$z <- cos(3 * d$x)
  objects <- list(npreg(y ~ x, data = d, bws = .3, bandwidth.compute = FALSE, residuals = TRUE),
    npscoef(y ~ z | x, data = d, bws = .3, bandwidth.compute = FALSE, residuals = TRUE),
    npplreg(y ~ z | x, data = d, bws = matrix(c(.3, .3), 2L, 1L), bandwidth.compute = FALSE, residuals = TRUE))
  for (object in objects) {
    object$bws$ydati$iuno[] <- TRUE
    object$bws$ydati$iord[] <- FALSE
    object$bws$ydati$icon[] <- FALSE
    expect_error(plot(object, output = "data", neval = 5), "must be a numeric ", fixed = TRUE)
  }
})
