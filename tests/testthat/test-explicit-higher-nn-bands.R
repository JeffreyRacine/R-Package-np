test_that("explicit higher Epanechnikov NN bands retain the fitted hat and the requested target", {
  skip_on_cran()
  opt <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(opt),add=TRUE)
  x <- data.frame(x=seq(-.9,.98,length.out=100))
  y <- sin(3*x$x)+.1*cos(seq_len(nrow(x)))
  b <- npregbw(xdat=x,ydat=y,bws=26,regtype="lp",degree=3L,basis="glp",
    bernstein.basis=TRUE,bwtype="generalized_nn",ckertype="epanechnikov",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y)
  a <- list(x=fit,errors="bootstrap",center="bias-corrected",neval=3L,B=99L,
    random.seed=146L,alpha=.1,output="data",gradients=TRUE,band="simultaneous")
  expect_error(do.call(plot,a))
  for (construction in c("symmetric-sd","basic-local-envelope")) {
    z <- suppressWarnings(do.call(plot,c(a,list(interval=construction))))
    m <- z$inference[[1L]]
    expect_identical(m$resampling_law,"mammen")
    expect_identical(m$components$simultaneous$interval,construction)
  }
  expect_false(.np_bias_explicit_nn_target(utils::modifyList(b,list(type="adaptive_nn")),TRUE))
  expect_false(.np_bias_explicit_nn_target(utils::modifyList(b,list(bernstein.basis=FALSE)),FALSE))
  st <- .np_bias_regression_prepare(b,x,y,adapt=TRUE,gradient.order=1L,law="mammen",
    calibration="absolute",explicit.nn=TRUE)
  g <- x[c(1L,35L,100L),,drop=FALSE]
  rows <- .np_bias_regression_rows(st,g,data.frame(row=seq_len(nrow(g)),target="d1"))
  H <- npreghat(b,txdat=x,exdat=g,s=1L,output="matrix")
  expect_equal(as.double(rows$original),as.double(H),tolerance=1e-10)
  expect_error(.np_bias_main_admit(st,matrix(1),1e-12,list(degree.engine=3L)),"zero-ridge")
  expect_error(.np_bias_regression_rows(st,data.frame(x=1),data.frame(row=1L,target="d1")),"observed training range")
})
