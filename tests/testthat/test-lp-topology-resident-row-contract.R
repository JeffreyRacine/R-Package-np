test_that("fixed LP topology objectives agree with independent hat oracles", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.acceleration = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260811L)
  n <- 97L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(4 * pi * x$x1) + 0.5 * x$x2 + rnorm(n, sd = 0.15)
  evaluate <- getFromNamespace(".npregbw_eval_only", "np")
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
    list(basis = "additive", bernstein = FALSE,
         degree = c(2L, 2L), kernel = "uniform"),
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
