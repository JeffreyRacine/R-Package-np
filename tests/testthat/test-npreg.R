test_that("npreg basic functionality works", {
  data("cps71")
  # Use a small subset for speed
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
  data("Italy")
  Italy_sub <- Italy[1:50, ]
  
  bw <- npregbw(gdp~ordered(year), data=Italy_sub, 
                bws=0.5, bandwidth.compute=FALSE)
  model <- npreg(bws=bw)
  
  expect_s3_class(model, "npregression")
  expect_output(summary(model))
})

test_that("npregbw cv.ls handles all-categorical fixed bandwidths", {
  old_rng <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_rng)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(5000001L)

  n <- 200L
  xbar <- factor(rbinom(n, 1L, 0.5), levels = 0:1)
  xtilde <- factor(rbinom(n, 1L, 0.5), levels = 0:1)
  y <- ifelse(as.integer(as.character(xbar)) == 1L, 1, 0) + rnorm(n)

  bw <- npregbw(
    xdat = data.frame(xbar = xbar, xtilde = xtilde),
    ydat = y,
    bwmethod = "cv.ls",
    bwscaling = FALSE,
    bwtype = "fixed",
    ukertype = "aitchisonaitken"
  )

  expect_equal(unname(bw$bw),
               c(0.016073908490923, 0.499999933944028),
               tolerance = 1e-10)
  expect_equal(unname(bw$fval), 1.0135341078806, tolerance = 1e-10)
  expect_true(is.finite(unname(bw$num.feval)))
  expect_gt(unname(bw$num.feval), 0)
})
