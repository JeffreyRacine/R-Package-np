test_that("integral kernel sums expose ordinary versus bandwidth-scaled CDFs", {
  z <- c(-1, 1)
  rhs <- c(2, 5)
  h <- 0.5
  expected <- drop(stats::pnorm(outer(0, z, "-") / h) %*% rhs)

  ordinary <- npksum(
    txdat = z,
    exdat = 0,
    tydat = rhs,
    operator = "integral",
    bws = h,
    bandwidth.divide = TRUE
  )$ksum
  scaled <- npksum(
    txdat = z,
    exdat = 0,
    tydat = rhs,
    operator = "integral",
    bws = h,
    bandwidth.divide = FALSE
  )$ksum

  ## The compiled Gaussian CDF approximation agrees with pnorm() to roughly
  ## 1e-10 on this platform; the factor-of-h defect is orders of magnitude
  ## larger than that numerical envelope.
  expect_equal(ordinary, expected, tolerance = 1e-10)
  expect_equal(scaled, h * expected, tolerance = 1e-10)
})

test_that("npregivderiv owns ordinary-CDF adjoint normalization", {
  real_npksum <- getFromNamespace("npksum", "np")
  integral_bandwidth_divide <- logical()

  local_mocked_bindings(
    npksum = function(...) {
      args <- list(...)
      if (identical(args$operator, "integral")) {
        integral_bandwidth_divide <<-
          c(integral_bandwidth_divide, args$bandwidth.divide)
      }
      do.call(real_npksum, args)
    },
    .package = "np"
  )

  set.seed(431)
  n <- 60L
  z <- runif(n, -1, 1)
  w <- z + rnorm(n, sd = 0.15)
  y <- z^2 + rnorm(n, sd = 0.08)

  fit <- suppressWarnings(npregivderiv(
    y, z, w,
    iterate.max = 2L,
    bandwidth.divide = FALSE
  ))

  expect_s3_class(fit, "npregivderiv")
  expect_length(integral_bandwidth_divide, 2L)
  expect_true(all(integral_bandwidth_divide))
})

test_that("npregivderiv monotonicity guard uses only computed norms", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_match(
    src,
    "is\\.monotone\\.increasing\\(norm\\.stop\\[seq_len\\(N\\)\\]\\)",
    perl = TRUE
  )
  expect_false(grepl(
    "!is\\.monotone\\.increasing\\(norm\\.stop\\)\\s*\\{",
    src,
    perl = TRUE
  ))
})

test_that("npregivderiv centers both adjoint terms on the fitted residual", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  fitted_second_term <- gregexpr(
    "S\\.z\\*mean\\.predicted\\.E\\.mu\\.w",
    src,
    perl = TRUE
  )[[1L]]

  expect_length(fitted_second_term[fitted_second_term > 0L], 2L)
  expect_false(grepl("S\\.z\\*mean\\.mu", src, perl = TRUE))
  expect_false(grepl("mean\\.mu <- mean\\(mu\\)", src, perl = TRUE))
})
