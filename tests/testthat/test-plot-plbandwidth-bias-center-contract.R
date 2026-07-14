test_that("partially linear bootstrap bias centers are component-specific", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 80L
  x1 <- rnorm(n)
  x2 <- factor(rbinom(n, 1L, 0.5))
  z1 <- factor(rbinom(n, 1L, 0.5))
  z2 <- rnorm(n)
  y <- 1 + x1 + as.integer(x2) + as.integer(z1) + sin(z2) + rnorm(n)

  model <- npplreg(
    y ~ x1 + x2 | z1 + z2,
    regtype = "lp",
    degree = 0L,
    nmulti = 1L
  )

  for (boot in c("wild", "inid", "fixed", "geom")) {
    set.seed(123)
    out <- suppressWarnings(plot(
      model,
      center = "bias-corrected",
      errors = "bootstrap",
      bootstrap = boot,
      band = "pointwise",
      B = 19L,
      data_overlay = FALSE,
      output = "data"
    ))

    expect_length(out, 4L)
    bc.lengths <- vapply(out, function(x) length(x$bias.corrected), integer(1L))
    mean.lengths <- vapply(out, function(x) length(x$mean), integer(1L))
    expect_equal(bc.lengths, mean.lengths)
    expect_true(all(vapply(out, function(x) all(is.finite(x$bias.corrected)), logical(1L))))

    continuous.panels <- vapply(out, function(x) length(x$bias.corrected) > 2L, logical(1L))
    expect_equal(sum(continuous.panels), 2L)
    continuous.centers <- lapply(out[continuous.panels], `[[`, "bias.corrected")
    expect_false(isTRUE(all.equal(continuous.centers[[1L]], continuous.centers[[2L]])))
  }
})

test_that("partially linear bias-corrected plot renders continuous and factor panels", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 80L
  x1 <- rnorm(n)
  x2 <- factor(rbinom(n, 1L, 0.5))
  z1 <- factor(rbinom(n, 1L, 0.5))
  z2 <- rnorm(n)
  y <- 1 + x1 + as.integer(x2) + as.integer(z1) + sin(z2) + rnorm(n)

  model <- npplreg(
    y ~ x1 + x2 | z1 + z2,
    regtype = "lp",
    degree = 0L,
    nmulti = 1L
  )

  withr::local_pdf(NULL)
  withr::local_options(list(plot.par.mfrow = FALSE))
  withr::local_par(list(mfrow = c(2L, 2L)))

  expect_error(
    suppressWarnings(plot(
      model,
      center = "bias-corrected",
      errors = "bootstrap",
      bootstrap = "wild",
      band = "simultaneous",
      B = 29L,
      data_overlay = FALSE
    )),
    NA
  )
})
