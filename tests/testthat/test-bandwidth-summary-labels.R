library(npRmpi)

test_that("fixed mixed-data bandwidth summaries label continuous scale factors correctly", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(123)
  n <- 40L
  xdat <- data.frame(
    x = runif(n),
    z1 = runif(n),
    z2 = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  ydat <- xdat$x + xdat$z1 + as.numeric(xdat$z2 == "b") + rnorm(n)

  bw <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lc",
    bwmethod = "cv.aic",
    nmulti = 0
  )

  s <- paste(capture.output(summary(bw)), collapse = "\n")

  expect_match(s, "Exp\\. Var\\. Name: x\\s+Bandwidth: .*Scale Factor:")
  expect_match(s, "Exp\\. Var\\. Name: z1\\s+Bandwidth: .*Scale Factor:")
  expect_match(s, "Exp\\. Var\\. Name: z2\\s+Bandwidth: .*Lambda Max:")
})
