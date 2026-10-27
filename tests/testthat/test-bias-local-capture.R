test_that("local calibration neighborhoods preserve units, rows and constant variance", {
  owner <- np:::.np_bias_local_neighborhood
  x <- matrix(seq(-1,1,length.out=24),ncol=1)
  q <- seq(.25,1,length.out=24)
  a <- owner(x,q)
  b <- owner(7+3*x,q)
  ii <- c(1L,12L,24L)
  expect_equal(a$raw(ii),b$raw(ii),tolerance=1e-13)
  expect_equal(drop(a$raw(ii) %*% (4*q)),rep(4,3),tolerance=1e-13)
  permutation <- rev(seq_len(24))
  c <- owner(x[permutation,,drop=FALSE],q[permutation])
  expect_equal(c$raw(match(ii,permutation))[,match(seq_len(24),permutation)],
               a$raw(ii),tolerance=1e-13)
  expect_error(a$raw(seq_len(17)),"invalid local calibration neighborhood tile")
  expect_error(owner(matrix(1,24,1),q),"invalid local calibration neighborhood geometry")
})

test_that("capture response identity follows numeric data and rejects changed values", {
  x <- data.frame(x=seq(-1,1,length.out=24))
  y <- setNames(1+x$x+x$x^2+.1*sin(seq_len(24)),rownames(x))
  b <- npregbw(xdat=x,ydat=y,bws=.4,bandwidth.compute=FALSE,
                regtype="lp",degree=2L,bernstein.basis=TRUE)
  state <- np:::.np_bias_regression_prepare(b,x,y,law="gaussian")
  grid <- x[c(1L,12L,24L),,drop=FALSE]
  fit <- np:::.np_plot_regression_eval(b,x,y,grid,need.bias.capture=TRUE)
  expect_s3_class(np:::.np_bias_capture_prepare(state,fit,grid),"np_bias_capture_prototype")
  changed <- fit; changed$fit.request$ydat[1L] <- changed$fit.request$ydat[1L]+1
  expect_error(np:::.np_bias_capture_prepare(state,changed,grid),"grid or response")
  changed$fit.request$ydat <- rev(fit$fit.request$ydat)
  expect_error(np:::.np_bias_capture_prepare(state,changed,grid),"grid or response")
})


test_that("raw main bases use their captured owner without changing auxiliary policy", {
  x <- data.frame(x=seq(-1,1,length.out=24))
  y <- x$x^2+.1*sin(seq_len(24))
  b <- npregbw(xdat=x,ydat=y,bws=.4,bandwidth.compute=FALSE,
               regtype="lp",degree=2L,bernstein.basis=FALSE)
  state <- np:::.np_bias_regression_prepare(b,x,y,law="gaussian")
  expect_identical(state$local.inference$version,"captured-pooled-positive-v1")
  expect_false(state$settings$main$bernstein.basis)
  expect_true(state$settings$bernstein.basis)
  expect_identical(state$settings$iterations,16L)
  grid <- x[c(1L,12L,24L),,drop=FALSE]
  fit <- np:::.np_plot_regression_eval(b,x,y,grid,gradients=TRUE,need.bias.capture=TRUE)
  capture <- np:::.np_bias_capture_prepare(state,fit,grid,c("mean","d1"))
  rows <- capture$rows(1:3,expand.grid(row=1:3,target=c("mean","d1"),stringsAsFactors=FALSE))
  expect_equal(rows$original_center,c(fit$mean,drop(fit$grad)),tolerance=1e-10)
  expect_true(all(is.finite(rows$se)))
  for (law in c("gaussian","mammen")) {
    expect_true(np:::.np_bias_local_supported(state,law))
    set.seed(42)
    one <- np:::.np_bias_inference(state,grid[2L,,drop=FALSE],"mean",B=63L,law=law,
      fit=np:::.np_plot_regression_eval(b,x,y,grid[2L,,drop=FALSE],need.bias.capture=TRUE))
    limits <- one$calibration$limits
    expect_identical(limits$pointwise,limits$simultaneous)
    expect_identical(limits$pointwise,limits$bonferroni)
  }
  for (law in c("rademacher","inid","fixed","geom"))
    expect_false(np:::.np_bias_local_supported(state,law))
  higher <- state; higher$gradient.order <- 2L
  expect_false(np:::.np_bias_local_supported(higher,"gaussian"))
  mixed <- state; mixed$factor.cols <- 2L
  expect_false(np:::.np_bias_local_supported(mixed,"gaussian"))
  multidimensional <- state; multidimensional$d <- 2L
  expect_false(np:::.np_bias_local_supported(multidimensional,"gaussian"))
  nn <- state; nn$bws$type <- "generalized_nn"
  expect_false(np:::.np_bias_local_supported(nn,"gaussian"))
})
