test_that("multivariate capture preserves fitted rows and variance policy", {
  i <- seq_len(80L)
  x <- data.frame(x=sin(i*sqrt(2)),z=cos(i*sqrt(3)))
  y <- x$x+.4*x$z^2+.1*sin(i*sqrt(5))
  bw <- npregbw(xdat=x,ydat=y,bws=c(.4,.6),bandwidth.compute=FALSE,
    regtype="lp",degree=c(1L,2L),bernstein.basis=FALSE)
  control <- np_bias_control(adaptation="none")
  state <- np:::.np_bias_regression_prepare(bw,x,y,control=control,law="gaussian")
  expect_identical(state$main.capture$version,"fixed-multivariate-fitted-v1")
  expect_null(state$local.inference)
  expect_identical(state$variance,np:::.np_bias_variance(state,state$residual))
  expect_identical(state$settings$iterations,16L)
  expect_identical(state$settings$main$degree,bw$degree)
  grid <- x[c(1L,40L,40L,80L),,drop=FALSE]
  fit <- np:::.np_plot_regression_eval(bw,x,y,grid,gradients=TRUE,need.bias.capture=TRUE)
  targets <- c("mean","d1","d2")
  capture <- np:::.np_bias_capture_prepare(state,fit,grid,targets)
  index <- c(4L,2L,3L,1L)
  map <- expand.grid(row=seq_along(index),target=targets,stringsAsFactors=FALSE)
  value <- capture$rows(index,map)
  expect_equal(value$original_center,c(fit$mean[index],fit$grad[index,1L],fit$grad[index,2L]),tolerance=1e-8)
  changed <- fit;changed$fit.request$ydat[1L] <- changed$fit.request$ydat[1L]+1
  expect_error(np:::.np_bias_capture_prepare(state,changed,grid,targets),"grid or response")
  expect_error(np:::.np_bias_capture_validate(capture,state,grid[4:1,,drop=FALSE],targets),"does not match")
  reg <- np:::npValidatedConditionalRegSpec(bw,where="test",ncon.field="ncon")
  expect_error(np:::.np_bias_main_admit(state,matrix(0,1L,80L),1e-6,reg,captured=TRUE),"regularized main rows")
  for(law in c("gaussian","mammen")) {
    expect_true(np:::.np_bias_main_capture_supported(state,law))
    set.seed(42)
    one <- np:::.np_bias_inference(state,grid[1L,,drop=FALSE],targets,B=31L,law=law,
      fit=np:::.np_plot_regression_eval(bw,x,y,grid[1L,,drop=FALSE],gradients=TRUE,need.bias.capture=TRUE))
    expect_identical(one$calibration$limits$pointwise,one$calibration$limits$simultaneous)
    expect_identical(one$calibration$limits$pointwise,one$calibration$limits$bonferroni)
    expect_identical(one$covariance$variance,"smoothed P0 residual variance")
  }
  for(law in c("rademacher","inid","fixed","geom"))
    expect_false(np:::.np_bias_main_capture_supported(state,law))
  single <- state;single$d <- 1L
  expect_false(np:::.np_bias_main_capture_supported(single,"gaussian"))
  for(field in c("ckertype","ckerorder","ckerbound")) {
    unsupported <- state
    unsupported$bws[[field]] <- switch(field,ckertype="epanechnikov",ckerorder=4L,ckerbound="fixed")
    expect_false(np:::.np_bias_main_capture_supported(unsupported,"gaussian"))
  }
  higher <- state;higher$gradient.order <- c(1L,2L)
  expect_false(np:::.np_bias_main_capture_supported(higher,"gaussian"))
  nn <- state;nn$bws$type <- "generalized_nn"
  expect_false(np:::.np_bias_main_capture_supported(nn,"gaussian"))
  mixed <- state;mixed$factor.cols <- 3L
  expect_false(np:::.np_bias_main_capture_supported(mixed,"gaussian"))
  adapted <- state;adapted$adaptation <- list(active=TRUE)
  expect_false(np:::.np_bias_main_capture_supported(adapted,"gaussian"))
  expect_error(np:::.np_bias_inference(state,grid,targets,B=31L,law="rademacher",main.capture=capture),"Gaussian or Mammen")
})
