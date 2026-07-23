test_that("LP CVAIC agrees with the independent exact-hat definition", {
  old_options <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(202607245L)
  n <- 193L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  signal <- sin(2 * pi * x$x1) * cos(2 * pi * x$x2) + 0.25 * x$x1 * x$x2
  y <- signal + rnorm(n, sd = 0.35 * sd(signal))
  evaluate_objective <- getFromNamespace(".npregbw_eval_only", "npRmpi")
  local_regression <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")

  cases <- list(
    list(type = "fixed", bws = c(0.22, 0.26), degree = c(1L, 1L),
         basis = "glp", bernstein = FALSE),
    list(type = "fixed", bws = c(0.28, 0.32), degree = c(2L, 2L),
         basis = "glp", bernstein = TRUE),
    list(type = "fixed", bws = c(0.30, 0.34), degree = c(2L, 3L),
         basis = "tensor", bernstein = FALSE),
    list(type = "generalized_nn", bws = c(65, 75), degree = c(1L, 1L),
         basis = "glp", bernstein = FALSE),
    list(type = "generalized_nn", bws = c(75, 85), degree = c(2L, 2L),
         basis = "glp", bernstein = FALSE),
    list(type = "adaptive_nn", bws = c(75, 85), degree = c(2L, 2L),
         basis = "glp", bernstein = FALSE)
  )

  for (case in cases) {
    bw <- local_regression(npregbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      bwmethod = "cv.aic",
      bwtype = case$type,
      ckertype = "gaussian",
      ckerorder = 2L,
      bws = case$bws,
      bandwidth.compute = FALSE,
      degree = case$degree,
      degree.select = "manual",
      basis = case$basis,
      bernstein.basis = case$bernstein
    ))

    native <- local_regression(
      evaluate_objective(x, y, bw, objective = "ls")
    )$objective
    hat <- local_regression(
      npreghat(bws = bw, txdat = x, output = "matrix")
    )
    fitted <- drop(hat %*% y)
    trace_hat <- sum(diag(hat))
    reference <- log(mean((y - fitted)^2)) +
      (1.0 + trace_hat / n) / (1.0 - (trace_hat + 2.0) / n)

    scale <- max(1.0, abs(reference))
    expect_lt(
      abs(native - reference),
      1e-8 * scale,
      label = paste(
        "type", case$type,
        "degree", paste(case$degree, collapse = ","),
        "basis", case$basis,
        "bernstein", case$bernstein
      )
    )
  }
})
