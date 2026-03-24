test_that("session explicit conditional bandwidth fits stay local and succeed", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")

  old.autodispatch <- getOption("npRmpi.autodispatch")
  old.messages <- getOption("np.messages")
  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)
  on.exit(options(npRmpi.autodispatch = old.autodispatch, np.messages = old.messages), add = TRUE)

  set.seed(20260324)
  n <- 24L
  x <- runif(n)
  y <- x + rnorm(n, sd = 0.15)
  dat <- data.frame(x = x, y = y)

  bw.cd <- npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lc",
    bws = c(0.8, 0.5),
    bandwidth.compute = FALSE
  )
  fit.cd <- npcdens(bws = bw.cd)
  expect_s3_class(fit.cd, "condensity")
  expect_true(length(fitted(fit.cd)) > 0L)

  bw.cdf <- npcdistbw(
    y ~ x,
    data = dat,
    regtype = "lc",
    bws = c(0.8, 0.5),
    bandwidth.compute = FALSE
  )
  fit.cdf <- npcdist(bws = bw.cdf)
  expect_s3_class(fit.cdf, "condistribution")
  expect_true(length(fitted(fit.cdf)) > 0L)
})
