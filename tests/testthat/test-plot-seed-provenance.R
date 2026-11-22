test_that("seed provenance distinguishes an explicit value from an injected default", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  withr::local_seed(510)
  withr::local_options(list(np.messages = FALSE, np.plot.progress = FALSE))
  d <- data.frame(x = runif(72L))
  b <- npudensbw(~x, data = d, bws = .17, bandwidth.compute = FALSE)
  fit <- npudens(b, tdat = d)
  counter <- 0L
  for (obj in list(b, fit)) {
    before <- .Random.seed
    a <- plot(obj, errors = "asymptotic", band = "simultaneous", output = "data", neval = 5L)
    b42 <- plot(obj, errors = "asymptotic", band = "simultaneous", output = "data", neval = 5L,
      random.seed = { counter <- counter + 1L; 42L })
    expect_false("random.seed" %in% a$request$supplied)
    expect_true("random.seed" %in% b42$request$supplied)
    expect_identical(a$panels, b42$panels)
    expect_identical(before, .Random.seed)
  }
  expect_identical(counter, 2L)
})
