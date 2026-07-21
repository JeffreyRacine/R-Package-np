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

test_that("npregivderiv adjoint dots remove only operator-owned names", {
  filter.dots <- getFromNamespace(".np_iv_deriv_adjoint_dots", "np")
  dots <- structure(
    list(1, "user-u", 2, "user-o", 3, FALSE, 4, 5),
    names = c("", "ukertype", "other", "okertype", NA_character_,
              "bandwidth.divide", "other", "")
  )

  observed <- filter.dots(dots)
  expected <- dots[c(1L, 3L, 5L, 7L, 8L)]

  expect_identical(observed, expected)
  expect_identical(filter.dots(unname(dots)), unname(dots))
})

test_that("categorical regression kernels do not collide with the adjoint", {
  real.npreg <- getFromNamespace("npreg", "np")
  real.npksum <- getFromNamespace("npksum", "np")
  regression.kernels <- list()
  adjoint.names <- list()
  adjoint.values <- list()

  local_mocked_bindings(
    npreg = function(...) {
      args <- list(...)
      regression.kernels[[length(regression.kernels) + 1L]] <<-
        args[c("ukertype", "okertype")]
      do.call(real.npreg, args)
    },
    npksum = function(...) {
      args <- list(...)
      if(identical(args$operator, "integral")) {
        adjoint.names[[length(adjoint.names) + 1L]] <<- names(args)
        adjoint.values[[length(adjoint.values) + 1L]] <<-
          args[c("bandwidth.divide", "ukertype", "okertype")]
      }
      do.call(real.npksum, args)
    },
    .package = "np"
  )

  set.seed(432)
  n <- 48L
  latent <- rnorm(n)
  z <- 0.4 * latent + rnorm(n, sd = 0.25)
  y <- z^2 + rnorm(n, sd = 0.06)
  w <- data.frame(
    wc = latent,
    wu = factor(ifelse(latent > 0, "high", "low")),
    wo = ordered(cut(latent, c(-Inf, -0.3, 0.3, Inf),
                     labels = c("low", "mid", "high")))
  )

  fit <- suppressWarnings(suppressMessages(npregivderiv(
    y, z, w, iterate.max = 2L, nmulti = 1L,
    ukertype = "aitchisonaitken", okertype = "wangvanryzin"
  )))

  expect_s3_class(fit, "npregivderiv")
  expect_true(length(regression.kernels) >= 4L)
  expect_true(all(vapply(
    regression.kernels,
    function(x) identical(x$ukertype, "aitchisonaitken") &&
      identical(x$okertype, "wangvanryzin"),
    logical(1L)
  )))
  expect_length(adjoint.names, 2L)
  expect_true(all(vapply(adjoint.names, function(x) {
    sum(x == "bandwidth.divide") == 1L &&
      sum(x == "ukertype") == 1L && sum(x == "okertype") == 1L
  }, logical(1L))))
  expect_true(all(vapply(adjoint.values, function(x) {
    identical(x$bandwidth.divide, TRUE) &&
      identical(x$ukertype, "liracine") &&
      identical(x$okertype, "liracine")
  }, logical(1L))))
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
