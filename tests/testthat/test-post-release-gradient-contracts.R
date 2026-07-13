test_that("npindex accepts only valid first-order gradient requests", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old.opts <- options(np.messages = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260713)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- x$x1 + x$x2 + rnorm(n, sd = 0.05)
  bw <- npindexbw(
    xdat = x,
    ydat = y,
    bws = c(1, 1, 0.25),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L
  )
  fit <- npindex(bws = bw, txdat = x, tydat = y, gradients = TRUE)

  expect_identical(gradients(fit, gradient.order = 1L), gradients(fit))
  expect_identical(gradients(fit, gradient.order = c(1L, 1L)), gradients(fit))

  invalid <- list(0L, -1L, 1.5, NA_real_, Inf, integer(0L), "2", TRUE,
                  c(1L, 2L))
  for (value in invalid) {
    expect_error(
      gradients(fit, gradient.order = value),
      "finite positive integers|only first-order"
    )
    expect_error(
      npindex(bws = bw, txdat = x, tydat = y, gradients = TRUE,
              gradient.order = value),
      "finite positive integers|only first-order"
    )
  }

  expect_error(
    plot(bw, xdat = x, ydat = y, output = "data", gradient_order = "2"),
    "positive numeric values|finite positive integers"
  )
})

test_that("nplsqreg gradient accessors enforce the stored derivative order", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old.opts <- options(np.messages = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260713)
  n <- 48L
  x <- data.frame(x1 = runif(n, -0.9, 0.9),
                  x2 = runif(n, -0.8, 0.8))
  y <- 0.5 + x$x1^3 - 0.7 * x$x2^2 + rnorm(n, sd = 0.03)
  bw <- nplsqregbw(
    xdat = x,
    ydat = y,
    bws = c(0.32, 0.32),
    bandwidth.compute = FALSE,
    scale = rep(0.35, n),
    delta = 0.5,
    tau = c(0.25, 0.75),
    regtype = "lp",
    degree = c(3L, 2L)
  )
  fit <- nplsqreg(bw, gradients = TRUE, gradient.order = c(2L, 1L))

  expect_identical(
    gradients(fit, gradient.order = c(2L, 1L)),
    fit$quantgrad
  )
  expect_identical(
    gradients(fit, errors = TRUE, gradient.order = c(2L, 1L)),
    fit$quantgerr
  )
  expect_error(
    gradients(fit, gradient.order = c(1L, 1L)),
    "differs from the derivative order stored"
  )
  expect_error(
    gradients(fit, gradient.order = c(4L, 1L)),
    "differs from the derivative order stored"
  )
  expect_error(gradients(fit, gradient.order = 0L), "positive integer")
  expect_error(gradients(fit, gradient.order = "2"), "finite positive integers")

  broken <- fit
  broken$fit[[1L]]$gradient.order <- c(1L, 1L)
  expect_error(
    gradients(broken, gradient.order = c(2L, 1L)),
    "stored lsqregression and child orders are inconsistent"
  )

  partial <- suppressWarnings(nplsqreg(
    bw,
    gradients = TRUE,
    gradient.order = c(4L, 1L)
  ))
  expect_true(all(is.na(partial$quantgrad[, 1L, ])))
  expect_true(any(is.finite(partial$quantgrad[, 2L, ])))
  expect_identical(
    suppressWarnings(gradients(partial, gradient.order = c(4L, 1L))),
    partial$quantgrad
  )
  expect_error(
    gradients(partial, gradient.order = c(2L, 1L)),
    "differs from the derivative order stored"
  )
})

test_that("nplsqreg plot treats gradient_order as an order, not a flag", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old.opts <- options(np.messages = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260714)
  n <- 40L
  x <- data.frame(x1 = runif(n, -0.9, 0.9),
                  x2 = runif(n, -0.8, 0.8))
  y <- x$x1^3 - x$x2^2 + rnorm(n, sd = 0.03)
  bw <- nplsqregbw(
    xdat = x,
    ydat = y,
    bws = c(0.32, 0.32),
    bandwidth.compute = FALSE,
    scale = rep(0.35, n),
    delta = 0.5,
    tau = c(0.25, 0.75),
    regtype = "lp",
    degree = c(3L, 2L)
  )
  fit <- nplsqreg(bw, gradients = TRUE, gradient.order = c(2L, 1L))

  out <- plot(
    fit,
    output = "data",
    neval = 7L,
    gradients = TRUE,
    gradient_order = c(2L, 1L),
    errors = "none",
    perspective = FALSE
  )
  expect_identical(names(out), c("cd1", "cd2"))
  expect_s3_class(out$cd1, "lsqregression")
  expect_identical(out$cd1$gradient.order, c(2L, 1L))
})
