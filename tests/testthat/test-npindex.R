test_that("npindex basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 100
  x1 <- runif(n)
  x2 <- runif(n)
  # Single index model: y = g(x1 + x2) + e
  y <- (x1 + x2)^2 + rnorm(n, sd=0.1)
  
  mydat <- data.frame(y, x1, x2)
  bw <- npindexbw(
    y ~ x1 + x2,
    data = mydat,
    bws = c(1, 0.35, 0.45),
    bandwidth.compute = FALSE,
    method = "ichimura"
  )
  
  model <- npindex(bws=bw)
  
  expect_s3_class(model, "singleindex")
  expect_type(predict(model), "double")
  expect_output(summary(model))
})

test_that("npindex public adaptive-nn lc route does not collapse to fixed semantics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(314161)
  n <- 70L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.06)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(20), , drop = FALSE]

  bw.fixed <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 1, 0.85),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed"
  )
  bw.adaptive <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 1, 0.85),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "adaptive_nn"
  )

  fit.fixed <- npindex(
    bws = bw.fixed,
    txdat = tx,
    tydat = y,
    exdat = ex,
    gradients = FALSE
  )
  fit.adaptive <- npindex(
    bws = bw.adaptive,
    txdat = tx,
    tydat = y,
    exdat = ex,
    gradients = FALSE
  )

  expect_gt(max(abs(as.vector(fit.fixed$mean) - as.vector(fit.adaptive$mean))), 1e-6)
})
