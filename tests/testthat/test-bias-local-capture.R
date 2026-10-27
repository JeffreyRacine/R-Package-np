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
