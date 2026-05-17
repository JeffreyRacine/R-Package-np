test_that("npcdens basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  # waiting is Y, eruptions is X
  bw <- npcdensbw(xdat=faithful$eruptions, ydat=faithful$waiting, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  fit <- npcdens(bws=bw)
  
  expect_s3_class(fit, "condensity")
  expect_type(predict(fit), "double")
  expect_equal(length(predict(fit)), nrow(faithful))
  
  expect_output(summary(fit))
})

test_that("npcdens works with formula", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("Italy")
  Italy_sub <- Italy[1:50, ]
  bw <- npcdensbw(gdp~ordered(year), data=Italy_sub, 
                                bws=c(0.5, 1.0), bandwidth.compute=FALSE)
  
  fit <- npcdens(bws=bw)
  expect_s3_class(fit, "condensity")
})

test_that("npcdens formula path handles non-syntactic variable names", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  dat <- data.frame(
    check.names = FALSE,
    "y var" = rnorm(8),
    "x var" = runif(8)
  )

  bw <- npcdensbw(
    `y var` ~ `x var`,
    data = dat,
    bws = c(0.5, 0.5),
    bandwidth.compute = FALSE
  )
  fit <- npcdens(bws = bw)

  expect_s3_class(bw, "conbandwidth")
  expect_s3_class(fit, "condensity")
  expect_identical(bw$variableNames$response, "y var")
  expect_identical(bw$variableNames$terms, "x var")
  expect_identical(length(fitted(fit)), nrow(dat))
})
