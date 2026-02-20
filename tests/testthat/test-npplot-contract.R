test_that("npplot contract: regression and conditional estimators return data payloads", {
  skip_if_not_installed("np")

  set.seed(101)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- rnorm(n)

  rbw <- npregbw(y ~ x + z, nmulti = 1)
  rout <- suppressWarnings(
    plot(
      rbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.num = 9,
      plot.errors.type = "all",
      view = "fixed"
    )
  )

  expect_type(rout, "list")
  expect_true(length(rout) > 0)
  expect_true(all(vapply(rout, inherits, logical(1), "npregression")))
  expect_true(all(vapply(rout, function(xi) !is.null(xi$mean), logical(1))))
  expect_true(all(vapply(rout, function(xi) !is.null(xi$merr), logical(1))))
})

test_that("npplot contract: npcdens/npcdist support bootstrap all in data mode", {
  skip_if_not_installed("np")

  set.seed(102)
  n <- 80
  x <- runif(n)
  y <- rnorm(n)

  cbw <- npcdensbw(y ~ x, nmulti = 1, bwmethod = "cv.ls")
  cout <- suppressWarnings(
    plot(
      cbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.num = 9,
      plot.errors.type = "all",
      view = "fixed"
    )
  )

  expect_type(cout, "list")
  expect_true(length(cout) > 0)
  expect_true(all(vapply(cout, inherits, logical(1), "condensity")))
  expect_true(all(vapply(cout, function(xi) !is.null(xi$condens), logical(1))))
  expect_true(all(vapply(cout, function(xi) !is.null(xi$conderr), logical(1))))

  dbw <- npcdistbw(y ~ x, nmulti = 1, bwmethod = "cv.ls")
  dout <- suppressWarnings(
    plot(
      dbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.num = 9,
      plot.errors.type = "all",
      view = "fixed"
    )
  )

  expect_type(dout, "list")
  expect_true(length(dout) > 0)
  expect_true(all(vapply(dout, inherits, logical(1), "condistribution")))
  expect_true(all(vapply(dout, function(xi) !is.null(xi$condist), logical(1))))
  expect_true(all(vapply(dout, function(xi) !is.null(xi$conderr), logical(1))))
})

test_that("npplot contract: plot.errors.alpha is enforced in (0,0.5)", {
  skip_if_not_installed("np")

  set.seed(103)
  n <- 60
  x <- runif(n)
  y <- rnorm(n)
  bw <- npregbw(y ~ x, nmulti = 1)

  expect_error(
    suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.alpha = 0.5
      )
    ),
    "must lie in \\(0,0\\.5\\)"
  )

  expect_error(
    suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.alpha = 0
      )
    ),
    "must lie in \\(0,0\\.5\\)"
  )
})
