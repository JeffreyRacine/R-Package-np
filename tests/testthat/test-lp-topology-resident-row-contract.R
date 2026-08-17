test_that("fixed LP topology objectives agree with independent hat oracles", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  old <- options(np.messages = FALSE, np.tree = FALSE, np.acceleration = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260811L)
  n <- 97L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(4 * pi * x$x1) + 0.5 * x$x2 + rnorm(n, sd = 0.15)
  evaluate <- getFromNamespace(".npregbw_eval_only", "npRmpi")
  cases <- list(
    list(basis = "additive", bernstein = FALSE,
         degree = c(1L, 1L), kernel = "gaussian"),
    list(basis = "additive", bernstein = TRUE,
         degree = c(1L, 1L), kernel = "epanechnikov"),
    list(basis = "tensor", bernstein = FALSE,
         degree = c(1L, 1L), kernel = "gaussian"),
    list(basis = "tensor", bernstein = TRUE,
         degree = c(1L, 1L), kernel = "epanechnikov"),
    list(basis = "additive", bernstein = FALSE,
         degree = c(1L, 1L), kernel = "uniform"),
    list(basis = "additive", bernstein = FALSE,
         degree = c(2L, 2L), kernel = "gaussian"),
    list(basis = "additive", bernstein = TRUE,
         degree = c(2L, 2L), kernel = "epanechnikov"),
    list(basis = "tensor", bernstein = TRUE,
         degree = c(4L, 0L), kernel = "gaussian"),
    list(basis = "additive", bernstein = TRUE,
         degree = c(2L, 2L), kernel = "beta", bound = "range")
  )

  for (case in cases) {
    bw <- npregbw(
      xdat = x, ydat = y, regtype = "lp", bwmethod = "cv.ls",
      bwtype = "fixed", ckertype = case$kernel,
      ckerbound = if (is.null(case$bound)) "none" else case$bound,
      bws = c(0.30, 0.30), bandwidth.compute = FALSE,
      degree = case$degree, degree.select = "manual",
      basis = case$basis, bernstein.basis = case$bernstein
    )
    native <- evaluate(x, y, bw, objective = "ls")$objective
    hat <- suppressWarnings(npreghat(bws = bw, txdat = x, output = "matrix"))
    fitted <- drop(hat %*% y)
    reference <- mean(((y - fitted) / (1.0 - diag(hat)))^2)

    expect_equal(
      native, reference, tolerance = 2e-10,
      info = paste(case$basis, case$bernstein, case$kernel)
    )
  }
})

test_that("ridged fixed LP CV agrees with the exact delete-one owner", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  old <- options(np.messages = FALSE, np.tree = FALSE, np.acceleration = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260811L)
  n <- 97L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(4 * pi * x$x1) + 0.5 * x$x2 + rnorm(n, sd = 0.15)
  bw <- npregbw(
    xdat = x, ydat = y, regtype = "lp", bwmethod = "cv.ls",
    bwtype = "fixed", ckertype = "uniform", bws = c(0.30, 0.30),
    bandwidth.compute = FALSE, degree = c(2L, 2L),
    degree.select = "manual", basis = "additive",
    bernstein.basis = FALSE
  )
  evaluate <- getFromNamespace(".npregbw_eval_only", "npRmpi")
  native <- evaluate(x, y, bw, objective = "ls")$objective
  loo.apply <- suppressWarnings(npreghat(
    bws = bw, txdat = x, y = y, output = "apply", leave.one.out = TRUE
  ))
  loo.matrix <- suppressWarnings(npreghat(
    bws = bw, txdat = x, output = "matrix", leave.one.out = TRUE
  ))
  ridge.used <- attr(loo.apply, "ridge.used", exact = TRUE)

  expect_gt(sum(ridge.used > 0.0), 0L)
  expect_equal(
    as.double(loo.apply), drop(loo.matrix %*% y), tolerance = 5e-14
  )
  expect_equal(native, mean((y - drop(loo.apply))^2), tolerance = 2e-10)
})
