test_that("mixed M1 keeps explicit none and the selected calibration law", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(2026102501)
  x <- data.frame(x=runif(80,-1,1),z=factor(rep(1:2,each=40)))
  y <- x$x+.2*as.integer(x$z)+rnorm(80,sd=.3)
  b <- npregbw(xdat=x,ydat=y,bws=c(8,.2),regtype="ll",bandwidth.compute=FALSE)
  make <- function(control=np_bias_control(),law="gaussian")
    np:::.np_bias_regression_prepare(b,x,y,control,adapt=TRUE,law=law)
  m <- make(); n <- make(np_bias_control(adaptation="none"))
  expect_true(m$settings$adaptation$active)
  expect_identical(m$settings$iterations,16L)
  expect_identical(m$settings$main$bandwidth,b[["bandwidth"]][["x"]])
  expect_null(n$mixed)
  expect_null(make(np_bias_control(p0.scale=4))$mixed)
  for(law in c("inid","rademacher","fixed","geom")) {
    v <- make(law=law)
    expect_null(v$mixed)
    expect_identical(v$residual,n$residual)
    expect_identical(v$variance,n$variance)
    expect_error(make(np_bias_control(adaptation="local-to-global"),law),"Gaussian or Mammen")
  }
})

test_that("M1 nested pilot has the specified polynomial global limit", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(2026102502)
  x <- data.frame(x=runif(80,-1,1),z=factor(rep(1:2,each=40)))
  b <- npregbw(xdat=x,ydat=x$x,bws=c(1e7,.2),regtype="ll",bandwidth.compute=FALSE)
  s <- np:::.np_bias_regression_prepare(b,x,x$x,np_bias_control(),adapt=TRUE,law="gaussian")
  # Independent OLS projection in each profile: no reference to M1's QR path.
  for(block in s$mixed$blocks) {
    X <- cbind(1,block$xcon[,1]); e <- block$x[2,,drop=FALSE]
    expected <- drop(c(1,e$x)%*%solve(crossprod(X),t(X)))
    H <- np:::.np_bias_mixed_pilot(block,e,block$geo$b0)$H
    expect_equal(drop(H),unname(expected),tolerance=1e-10)
    expect_equal(sum(H),1,tolerance=1e-12)
  }
})
