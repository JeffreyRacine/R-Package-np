test_that("sampled cumulative and total quadrature have distinct contracts", {
  env <- environment(npuniden.boundary)
  cumulative <- get("integrate.trapezoidal", env)
  total <- get(".np_quadrature_total", env)
  prepare <- get(".np_quadrature_prepare", env)
  apply <- get(".np_quadrature_cumulative", env)
  x <- c(0,.1,.5,1)
  expected <- c(0,.0005,.0525,.365)
  expect_equal(cumulative(x,x^2),expected,tolerance=1e-15)
  expect_equal(total(x,x^2),.365,tolerance=1e-15)
  expect_equal(cumulative(rep(x,each=2),rep(x^2,each=2)),
               rep(cumulative(x,x^2),each=2),tolerance=0)
  permutation <- c(4,2,1,3,2)
  expect_equal(cumulative(x[permutation],x[permutation]^2),
               expected[permutation],tolerance=1e-15)
  expect_error(cumulative(c(0,0,1),c(0,1,1)),"conflicting quadrature")
  expect_identical(cumulative(numeric(),numeric()),numeric())
  expect_equal(cumulative(3,5),0)
  expect_equal(cumulative(c(3,3),c(5,5)),c(0,0))
  expect_equal(cumulative(c(2,4),c(1,3)),c(0,4))
  expect_equal(total(c(2,4),c(1,3)),4)
  expect_error(cumulative(c(0,Inf),c(0,1)),"finite numeric")
  expect_error(cumulative(c(0,1),1),"match the coordinates")
  geometry <- prepare(rev(x))
  expect_equal(apply(rev(x^2),geometry),rev(expected),tolerance=1e-15)
  expect_equal(apply(rep(2,4),geometry),rev(2*x),tolerance=1e-15)
  # Exact former total on an admissible uniform grid, including entropy.
  grid <- seq(-4,5,length.out=257)
  value <- exp(-grid^2/3)*(1+cos(grid)^2)
  step <- diff(grid)[1]
  dy <- diff(value)
  old.total <- tail(cumsum(diff(grid)*(head(value,-1)+tail(value,-1))/2),1) -
    step^2/12*(tail(dy,1)/step-dy[1]/step)
  expect_identical(unname(total(grid,value)),unname(old.total))
  expect_identical(total(rep(grid,each=2),rep(value,each=2)),total(grid,value))
  expect_equal(cumulative(grid,value)[1],0,tolerance=0)
  # Translation and scale do not turn an irregular grid into a uniform one.
  expect_false(prepare(1e12+x)$uniform)
  expect_false(prepare(1e-12*x)$uniform)
  expect_equal(total(1e-12*x,x^2),1e-12*.365,tolerance=1e-15)
})

test_that("beta kernels use dimensionless h and published beta2 regions", {
  X <- seq(.05,.95,length.out=15)
  Y <- c(0,.1,.2-.Machine$double.eps,.2,.5,.8,.9,1)
  for(h in c(.1,.25-.Machine$double.eps,.25)) {
    oracle <- vapply(Y,function(x) {
      rho <- function(z) 2*h^2+2.5-sqrt(4*h^4+6*h^2+2.25-z^2-z/h)
      if(x<2*h) mean(dbeta(X,rho(x),(1-x)/h))
      else if(x>1-2*h) mean(dbeta(X,x/h,rho(1-x)))
      else mean(dbeta(X,x/h,(1-x)/h))
    },0.)
    fit <- npuniden.boundary(X,Y,h=h,a=0,b=1,kertype="beta2")
    expect_equal(fit$f,oracle,tolerance=1e-14)
    for(width in c(.05,5)) {
      shifted <- npuniden.boundary(3+width*X,3+width*Y,h=h,a=3,b=3+width,
                                   kertype="beta2")
      expect_equal(width*shifted$f,fit$f,tolerance=2e-13)
      expect_equal(width*shifted$sd.f,fit$sd.f,tolerance=2e-13)
    }
  }
  for(h in c(.25+.Machine$double.eps,.3,Inf)) {
    expect_error(npuniden.boundary(X,Y,h=h,a=0,b=1,kertype="beta2"),
                 "0 < h <= 1/4",fixed=TRUE)
    expect_error(npuniden.boundary(X,Y,grid=c(.1,h),a=0,b=1,kertype="beta2"),
                 "0 < h <= 1/4",fixed=TRUE)
  }
  expect_error(npuniden.boundary(X,Y,h=.1,a=0,b=Inf,kertype="beta1"),
               "finite bounds")
  for(h in c(.1,1)) {
    fit <- npuniden.boundary(X,Y,h=h,a=0,b=1,kertype="beta1")
    shifted <- npuniden.boundary(2+.05*X,2+.05*Y,h=h,a=2,b=2.05,kertype="beta1")
    expect_equal(.05*shifted$f,fit$f,tolerance=2e-13)
  }
  # Skewed, rather than near-uniform, data exercise an interior CV optimum.
  X <- qbeta(seq(.02,.98,length.out=25),2,7)
  for(method in c("cv.ls","cv.ml")) {
    fit <- npuniden.boundary(X,Y,a=0,b=1,kertype="beta2",bwmethod=method)
    expect_true(is.finite(fit$h) && fit$h>0 && fit$h<=.25)
    explicit <- npuniden.boundary(X,Y,grid=c(.01,.05,.2,.25),a=0,b=1,
                                 kertype="beta2",bwmethod=method)
    expect_true(is.finite(explicit$h) && explicit$h>0 && explicit$h<=.25)
  }
})

test_that("sampled quadrature preserves training multiplicity", {
  X <- c(.1,.15,.3,.6,.9)
  Y <- c(.2,.4,.8)
  repeated <- c(X,X[1:2])
  fit <- npuniden.boundary(repeated,Y,h=.2,a=0,b=1)
  expected <- vapply(Y,function(y)
    mean(dnorm((y-repeated)/.2)/(.2*(pnorm((1-y)/.2)-pnorm(-y/.2)))),0.)
  expect_identical(fit$f,expected)
  expect_false(isTRUE(all.equal(fit$f,npuniden.boundary(X,Y,h=.2,a=0,b=1)$f)))
  # An inactive QP has identical ordinates at repeated integration nodes.
  fit <- npuniden.sc(X,Y=c(Y,Y),h=.3,lb=0,ub=10,integral.equal=TRUE)
  expect_true(fit$solve.QP)
  expect_equal(fit$F[1:3],fit$F[4:6],tolerance=0)
  expect_equal(fit$F.sc[1:3],fit$F.sc[4:6],tolerance=0)
})
