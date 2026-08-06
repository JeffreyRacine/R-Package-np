request_driven_regression_runtime <- function() {
  if (exists("spawn_mpi_slaves", mode = "function") &&
      !spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")
}

test_that("npreg request states preserve retained scalar and LP outputs", {
  request_driven_regression_runtime()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026080601L)
  n <- 96L
  x <- data.frame(x1 = runif(n, -0.8, 0.8),
                  x2 = runif(n, -0.7, 0.9))
  y <- sin(1.7 * x$x1) + 0.35 * x$x2 + rnorm(n, sd = 0.08)
  ex <- x[seq_len(23L), , drop = FALSE]

  bw0 <- npregbw(
    xdat = x["x1"], ydat = y, bws = 0.24,
    bandwidth.compute = FALSE, regtype = "lp", degree = 0L
  )
  full0 <- npreg(bws = bw0, txdat = x["x1"], tydat = y,
                 errors = TRUE)
  lean0 <- npreg(bws = bw0, txdat = x["x1"], tydat = y,
                 errors = FALSE)
  expect_identical(lean0$mean, full0$mean)
  expect_null(lean0$merr)
  expect_false(lean0$errors)
  expect_error(se(lean0), "errors=TRUE", fixed = TRUE)

  bw <- npregbw(
    xdat = x, ydat = y, bws = c(0.32, 0.36),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(2L, 1L),
    degree.select = "manual", basis = "glp", bernstein.basis = FALSE
  )
  full <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                errors = TRUE, gradients = TRUE)
  mean.only <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                     errors = FALSE, gradients = FALSE)
  mean.se <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                   errors = TRUE, gradients = FALSE)
  mean.grad <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                     errors = FALSE, gradients = TRUE)

  expect_equal(mean.only$mean, full$mean, tolerance = 2e-12)
  expect_identical(mean.se$mean, full$mean)
  expect_identical(mean.se$merr, full$merr)
  expect_equal(mean.grad$mean, full$mean, tolerance = 2e-12)
  expect_equal(mean.grad$grad, full$grad, tolerance = 2e-11)
  expect_null(mean.only$merr)
  expect_null(mean.only$grad)
  expect_null(mean.only$gerr)
  expect_null(mean.grad$merr)
  expect_null(mean.grad$gerr)
  expect_error(gradients(mean.grad, errors = TRUE),
               "gradients=TRUE and errors=TRUE", fixed = TRUE)
})

test_that("beta request suppression agrees with the independent LP apply owner", {
  request_driven_regression_runtime()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026080602L)
  n <- 96L
  x <- data.frame(x1 = runif(n, -0.8, 0.8),
                  x2 = runif(n, -0.7, 0.9))
  y <- sin(1.7 * x$x1) + 0.35 * x$x2 + rnorm(n, sd = 0.08)
  ex <- x[seq_len(19L), , drop = FALSE]
  bw <- npregbw(
    xdat = x, ydat = y, bws = c(0.18, 0.20),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(2L, 1L),
    degree.select = "manual", basis = "glp", bernstein.basis = TRUE,
    ckertype = "beta", ckerbound = "range",
    ckerlb = c(-1, -1), ckerub = c(1, 1)
  )

  lean <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                errors = FALSE, gradients = TRUE)
  oracle <- function(s) drop(.npreghat_exact_lp_apply_from_regression_core(
    bws = bw, txdat = x, y = matrix(y, ncol = 1L), exdat = ex,
    basis = "glp", degree = c(2L, 1L), bernstein.basis = TRUE,
    s = as.integer(s), return.hat = FALSE
  ))

  expect_equal(lean$mean, oracle(c(0L, 0L)), tolerance = 2e-10)
  expect_equal(lean$grad[, 1L], oracle(c(1L, 0L)), tolerance = 2e-8)
  expect_equal(lean$grad[, 2L], oracle(c(0L, 1L)), tolerance = 2e-8)
})

test_that("all-categorical request suppression honors compression on and off", {
  request_driven_regression_runtime()
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026080603L)
  n <- 90L
  x <- data.frame(u = factor(sample(letters[1:4], n, TRUE)),
                  o = ordered(sample(1:3, n, TRUE)))
  y <- rnorm(n)
  bw <- npregbw(xdat = x, ydat = y, bws = c(0.2, 0.4),
                bandwidth.compute = FALSE, regtype = "lc")

  for (compress in c(FALSE, TRUE)) {
    options(np.categorical.compress = compress)
    full <- npreg(bws = bw, txdat = x, tydat = y, errors = TRUE)
    lean <- npreg(bws = bw, txdat = x, tydat = y, errors = FALSE)
    expect_identical(lean$mean, full$mean)
    expect_null(lean$merr)
  }
})

test_that("prediction and request boundary contracts are explicit", {
  request_driven_regression_runtime()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026080604L)
  x <- data.frame(x = runif(80L))
  y <- sin(2 * pi * x$x) + rnorm(80L, sd = 0.1)
  bw <- npregbw(xdat = x, ydat = y, bws = 0.18,
                bandwidth.compute = FALSE, regtype = "lp", degree = 1L)
  fit <- npreg(bws = bw, txdat = x, tydat = y, errors = FALSE)

  expect_error(npreg(bws = bw, txdat = x, tydat = y, errors = NA),
               "errors")
  expect_identical(predict(fit, newdata = x[1:7, , drop = FALSE]),
                   fitted(npreg(bws = bw, txdat = x, tydat = y,
                                exdat = x[1:7, , drop = FALSE],
                                errors = FALSE)))
  pred.se <- predict(fit, newdata = x[1:7, , drop = FALSE], se.fit = TRUE)
  expect_named(pred.se, c("fit", "se.fit", "df", "residual.scale"))
  expect_length(pred.se$fit, 7L)
  expect_length(pred.se$se.fit, 7L)
})

test_that("request contract has one closed native state space", {
  header <- readLines(test_path("..", "..", "src", "regression_contract.h"),
                      warn = FALSE)
  source <- readLines(test_path("..", "..", "src", "np.c"), warn = FALSE)

  expect_true(any(grepl("NP_REGRESSION_OUTPUT_MEAN_ONLY = 0", header,
                        fixed = TRUE)))
  expect_true(any(grepl("NP_REGRESSION_OUTPUT_FULL = 3", header,
                        fixed = TRUE)))
  expect_true(any(grepl("np_regression_output_request_valid", source,
                        fixed = TRUE)))
})
