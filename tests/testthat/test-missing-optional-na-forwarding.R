# Broad family/transport matrices are retained in the R84 campaign evidence.
test_that("missing optional NA arguments preserve omission without swallowing policies", {
  if (environmentName(environment(npreg)) == "npRmpi") {
    skip_on_cran()
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  withr::local_options(np.messages = FALSE)
  d <- data.frame(x = seq(-1, 1, length.out = 24L))
  d$y <- sin(3 * d$x)
  d$y[4L] <- NA_real_
  calls <- new.env(parent = emptyenv())
  calls$data <- calls$factory <- calls$policy <- 0L
  factory <- function() {
    calls$factory <- calls$factory + 1L
    function(x) {
      calls$policy <- calls$policy + 1L
      stats::na.omit(x)
    }
  }
  fresh <- function(na.action)
    npreg(y ~ x, data = { calls$data <- calls$data + 1L; d },
          na.action = na.action, bws = .4, bandwidth.compute = FALSE)
  forwarded <- function(na.action) fresh(na.action = na.action)
  expected <- npreg(y ~ x, data = d, bws = .4, bandwidth.compute = FALSE)
  set.seed(8404)
  seed <- .Random.seed
  expect_equal(fitted(forwarded()), fitted(expected), tolerance = 0)
  expect_identical(.Random.seed, seed)
  calls$data <- 0L
  expect_equal(fitted(fresh(factory())), fitted(expected), tolerance = 0)
  expect_identical(calls$data, 1L)
  expect_identical(calls$factory, 1L)
  expect_identical(calls$policy, 1L)

  b <- expected$bws
  retained <- function(na.action) npreg(b, na.action = na.action)
  replacement <- function(na.action) npreg(b, data = d, na.action = na.action)
  expect_equal(fitted(retained()), fitted(npreg(b)), tolerance = 0)
  expect_equal(fitted(replacement()), fitted(npreg(b, data = d)), tolerance = 0)
  calls$factory <- 0L
  expect_error(retained(factory()), "na.action override is not supported")
  expect_identical(calls$factory, 0L)
  expect_error(retained(NULL), "na.action override is not supported")

  native <- function(na.action)
    npregbw(xdat = d$x, ydat = replace(d$y, 4L, 0),
            bws = .4, bandwidth.compute = FALSE, na.action = na.action)
  expect_identical(native()$bw,
    npregbw(xdat = d$x, ydat = replace(d$y, 4L, 0),
            bws = .4, bandwidth.compute = FALSE)$bw)
  expect_error(native(factory()), "supported only with the formula interface")
  expect_identical(calls$factory, 0L)

  test <- function(na.action)
    npsigtest(b, data = d, na.action = na.action, B = 9L, boot.method = "iid")
  set.seed(8405)
  got <- test()
  seed <- .Random.seed
  set.seed(8405)
  ref <- npsigtest(b, data = d, B = 9L, boot.method = "iid")
  fields <- c("In", "In.bootstrap", "P", "bootstrap.executed")
  expect_equal(got[fields], ref[fields], tolerance = 0)
  expect_identical(.Random.seed, seed)
})


test_that("native refusal checks use original dot positions", {
  pkg <- if ("package:npRmpi" %in% search()) "npRmpi" else "np"
  if (pkg == "npRmpi") {
    skip_on_cran()
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  auxiliary <- getFromNamespace(".npregbw_auxiliary", pkg)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  policy <- function() { counter$n <- counter$n + 1L; stats::na.omit }
  expect_error(auxiliary(xdat = 1:8, ydat = 1:8, bws = .4,
    bandwidth.compute = FALSE, na.action = policy()),
    "supported only with the formula interface")
  expect_identical(counter$n, 0L)
  d <- data.frame(x = seq(0, 1, length.out = 24L))
  f <- function(na.action) npudensbw(dat = d, bws = .4,
    bandwidth.compute = FALSE, na.action = na.action)
  expect_identical(f()$bw,
    npudensbw(dat = d, bws = .4, bandwidth.compute = FALSE)$bw)
})


test_that("default regression fit omits a forwarded missing policy", {
  if (identical(environmentName(environment(npreg)), "npRmpi")) {
    skip_on_cran()
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  i <- seq_len(24L)
  d <- data.frame(x = sin(i / 5), y = cos(i / 3))
  f <- function(na.action) npreg(txdat = d["x"], tydat = d$y,
    bandwidth.compute = FALSE, na.action = na.action)
  expected <- npreg(txdat = d["x"], tydat = d$y, bandwidth.compute = FALSE)
  expect_identical(fitted(f()), fitted(expected))
})
