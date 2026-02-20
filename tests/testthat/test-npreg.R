test_that("npreg basic functionality works", {
  # # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("cps71")
  cps71_sub <- cps71[1:100, ]
  # Local constant
  bw_lc <- npregbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                                 bws=1.0, bandwidth.compute=FALSE, regtype="lc")
  model_lc <- npreg(bws=bw_lc)
  
  expect_s3_class(model_lc, "npregression")
  expect_type(predict(model_lc), "double")
  expect_equal(length(predict(model_lc)), 100)
  
  # Local linear
  bw_ll <- npregbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                                 bws=1.0, bandwidth.compute=FALSE, regtype="ll")
  model_ll <- npreg(bws=bw_ll)
  
  expect_s3_class(model_ll, "npregression")
  expect_type(predict(model_ll), "double")
})

test_that("npreg works with formula and mixed data", {
  # # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("Italy")
  Italy_sub <- Italy[1:50, ]
  bw <- npregbw(gdp~ordered(year), data=Italy_sub, 
                              bws=0.5, bandwidth.compute=FALSE)
  model <- npreg(bws=bw)
  
  expect_s3_class(model, "npregression")
  expect_output(summary(model))
})