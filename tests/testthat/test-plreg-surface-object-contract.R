test_that("partially linear surface payloads retain named dimensions", {
  set.seed(752)
  n <- 61L
  d <- data.frame(x = rnorm(n), z = runif(n))
  d$y <- 2 * d$x + cos(d$z) + rnorm(n, sd = .2)
  bw <- npplregbw(xdat = d["x"], zdat = d["z"], ydat = d$y,
                  bws = matrix(.55, 2, 1), bandwidth.compute = FALSE)
  fit <- npplreg(bw, se = TRUE)
  for (common in c(FALSE, TRUE)) {
    for (surface in c(FALSE, TRUE)) {
      payload <- plot(fit, output = "data", neval = 5L,
                      errors = "none", perspective = surface,
                      common.scale = common)
      expect_s3_class(payload, "np_plot_data")
      expect_false(any(c("R2", "MSE", "fit") %in% names(payload)))
      expect_output(print(payload), "Plot data: npplreg")
      for (j in seq_along(payload$panels)) {
        one <- as.data.frame(payload, panel = j)
        ex <- plot_coordinates(payload, j, "X")
        ez <- plot_coordinates(payload, j, "Z")
        expect_identical(nrow(one), if (surface) 25L else 5L)
        expect_identical(names(ex), "x")
        expect_identical(names(ez), "z")
        ref <- npplreg(bw, exdat = ex, ezdat = ez)
        expect_equal(one$estimate, fitted(ref), tolerance = 0)
      }
    }
  }
})

test_that("native PL plot reconstruction preserves non-syntactic names", {
  set.seed(892)
  x <- data.frame(value = rnorm(63)); z <- data.frame(value = runif(63))
  names(x) <- "linear x"; names(z) <- "smooth z"
  y <- 2*x[[1]] + cos(z[[1]]) + rnorm(63, sd = .2)
  b <- npplregbw(xdat = x, zdat = z, ydat = y,
                 bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  fit <- npplreg(b, se = TRUE)
  for (obj in list(b, fit)) for (surface in c(FALSE, TRUE)) {
    implicit <- plot(obj, output = "data", perspective = surface,
                     neval = 5L, errors = "none")
    explicit <- plot(b, xdat = x, ydat = y, zdat = z,
                     output = "data", perspective = surface,
                     neval = 5L, errors = "none")
    for (j in seq_along(implicit$panels)) {
      expect_identical(names(plot_coordinates(implicit, j, "X")), names(x))
      expect_identical(names(plot_coordinates(implicit, j, "Z")), names(z))
      expect_equal(implicit$panels[[j]]$grid, explicit$panels[[j]]$grid, tolerance = 0)
      expect_equal(implicit$panels[[j]]$values, explicit$panels[[j]]$values, tolerance = 0)
    }
    # Coefficients belong to the fitted model, never to a plot-data facsimile.
    expect_identical(names(coef(fit)), names(x))
  }
  names(x) <- "replacement x"; names(z) <- "replacement z"
  replacement <- plot(b, xdat = x, ydat = y, zdat = z,
    output = "data", perspective = TRUE, neval = 5L, errors = "none")
  expect_identical(names(plot_coordinates(replacement, 1L, "X")), names(x))
  expect_identical(names(plot_coordinates(replacement, 1L, "Z")), names(z))
})
