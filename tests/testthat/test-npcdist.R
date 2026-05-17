test_that("npcdist basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  bw <- npcdistbw(xdat=faithful$eruptions, ydat=faithful$waiting, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  fit <- npcdist(bws=bw)
  
  expect_s3_class(fit, "condistribution")
  expect_type(predict(fit), "double")
  expect_equal(length(predict(fit)), nrow(faithful))
  
  expect_output(summary(fit))
})

test_that("npcdist formula path handles non-syntactic variable names", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  dat <- data.frame(
    check.names = FALSE,
    "y var" = runif(8),
    "x var" = runif(8)
  )

  bw <- npcdistbw(
    `y var` ~ `x var`,
    data = dat,
    bws = c(0.5, 0.5),
    bandwidth.compute = FALSE
  )
  fit <- npcdist(bws = bw)

  expect_s3_class(bw, "condbandwidth")
  expect_s3_class(fit, "condistribution")
  expect_identical(bw$variableNames$response, "y var")
  expect_identical(bw$variableNames$terms, "x var")
  expect_identical(length(fitted(fit)), nrow(dat))
})
