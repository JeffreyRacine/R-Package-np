test_that("npplreg lp degree metadata is preserved in plbandwidth summary", {
  if (!spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 120
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  z1 <- rbinom(n, 1, 0.5)
  z2 <- rnorm(n)
  y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

  for (deg in c(0L, 2L)) {
    fit <- npplreg(
      y ~ x1 + factor(x2) | factor(z1) + z2,
      regtype = "lp",
      degree = deg
    )

    expect_identical(as.integer(fit$bw$degree), as.integer(deg))
    expect_identical(as.integer(fit$bw$bw$yzbw$degree), as.integer(deg))

    s <- capture.output(summary(fit$bw))
    reg_line <- s[grep("Regression Type:", s)][1]
    expect_true(grepl(paste0("degree = ", deg), reg_line))
  }
})
