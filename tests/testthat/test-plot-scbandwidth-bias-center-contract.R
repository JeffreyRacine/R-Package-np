test_that("smooth coefficient bootstrap bias centers are component-specific", {
  set.seed(42)
  n <- 80L
  x <- runif(n)
  z <- runif(n, min = -2, max = 2)
  y <- x * exp(z) * (1 + rnorm(n, sd = 0.2))

  model <- npscoef(
    y ~ x | z,
    regtype = "lp",
    degree = 0L,
    bws = 2
  )

  set.seed(123)
  out <- suppressWarnings(plot(
    model,
    center = "bias-corrected",
    errors = "bootstrap",
    bootstrap = "wild",
    band = "pointwise",
    B = 19L,
    data_overlay = FALSE,
    perspective = FALSE,
    output = "data"
  ))

  expect_length(out, 2L)
  expect_equal(
    vapply(out, function(x) length(x$bias.corrected), integer(1L)),
    vapply(out, function(x) length(x$mean), integer(1L))
  )
  expect_true(all(vapply(out, function(x) all(is.finite(x$bias.corrected)), logical(1L))))
  expect_false(isTRUE(all.equal(out[[1L]]$bias.corrected, out[[2L]]$bias.corrected)))
})

test_that("smooth coefficient bias-corrected plot renders center legend", {
  set.seed(42)
  n <- 80L
  x <- runif(n)
  z <- runif(n, min = -2, max = 2)
  y <- x * exp(z) * (1 + rnorm(n, sd = 0.2))

  model <- npscoef(
    y ~ x | z,
    regtype = "lp",
    degree = 0L,
    bws = 2
  )

  pdf.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf.file)
  on.exit(grDevices::dev.off(), add = TRUE)

  options(plot.par.mfrow = FALSE)
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par), add = TRUE)
  graphics::par(mfrow = c(1L, 2L))

  expect_error(
    suppressWarnings(plot(
      model,
      center = "bias-corrected",
      errors = "bootstrap",
      bootstrap = "wild",
      band = "simultaneous",
      B = 29L,
      data_overlay = FALSE,
      perspective = FALSE
    )),
    NA
  )
})
