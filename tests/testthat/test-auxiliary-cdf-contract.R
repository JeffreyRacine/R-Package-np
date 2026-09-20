test_that("continuous integrals retain absolute origins and infinite tails", {
  q <- .np_density_integral(function(x) 2*x,c(.8,0,.2,.8,1),0,1,c(.2,.8),.1)
  expect_equal(q$F,c(.64,0,.04,.64,1),tolerance=1e-12)
  expect_equal(q$total,1,tolerance=1e-12)
  for(bounds in list(c(-Inf,Inf),c(0,Inf),c(-Inf,0))) {
    a <- bounds[1];b <- bounds[2]
    y <- seq(max(-3,a),min(3,b),length.out=21)
    x <- y[seq(2,20,by=2)]
    q <- .np_density_integral(dnorm,y,a,b,x,.3)
    expect_equal(q$F,pnorm(y)-pnorm(a),tolerance=2e-8)
    expect_equal(q$total,pnorm(b)-pnorm(a),tolerance=2e-8)
  }
  # Wide empty gaps and a narrow distant mode must not hide probability mass.
  x <- c(-100,100);h <- .001
  f <- function(y) (dnorm((y-x[1])/h)+dnorm((y-x[2])/h))/(2*h)
  y <- c(-100.001,-100,0,100,100.001)
  q <- .np_density_integral(f,y,-101,101,x,h)
  expect_equal(q$F,(pnorm((y-x[1])/h)+pnorm((y-x[2])/h))/2,tolerance=2e-8)
  expect_equal(q$total,1,tolerance=2e-8)
  expect_error(.np_density_integral(function(x) rep(NaN,length(x)),.5,0,1,.5,.2),
               "non-finite continuous")
  expect_error(.np_density_integral(function(x)stop("must not evaluate"),.5,0,1,
    seq(0,1,length.out=5200),1e-7),"evaluation budget")
  shifted <- .np_density_integral(function(x)dnorm((x-7)/3)/3,
    c(4,7,10),-Inf,Inf,c(4,10),.9)
  expect_equal(shifted$F,pnorm(c(-1,0,1)),tolerance=2e-8)
})

test_that("positive quadrature has nonnegative mass and one coherent primitive", {
  y <- seq(0,1,length.out=1001)
  f <- function(x) pmax(0,1-4*x)
  q <- .np_density_integral(f,y,0,1,c(.1,.2),.1,positive=TRUE)
  oracle <- ifelse(y<.25,y-2*y^2,.125)
  expect_equal(q$F,oracle,tolerance=2e-8)
  expect_equal(q$total,.125,tolerance=2e-8)
  expect_gte(min(diff(q$F)),-1e-14)
  expect_gte(min(q$F),-1e-14)
  expect_lte(max(q$F),q$total+1e-14)
  expect_error(.np_density_integral(function(x) -x,.5,0,1,.5,.2,positive=TRUE),
               "negative integrand")
})

test_that("boundary CDFs integrate the retained density, not the requested grid", {
  X <- seq(.05,.95,length.out=40);h <- .12;Y <- c(.75,.25,.5,.25)
  density <- function(y) vapply(y,function(z)
    mean(dnorm((z-X)/h))/(h*(pnorm((1-z)/h)-pnorm(-z/h))),0.0)
  oracle <- vapply(Y,function(y) integrate(density,0,y,rel.tol=1e-11)$value,0.0)
  raw <- npuniden.boundary(X,Y,h=h,a=0,b=1)
  proper <- npuniden.boundary(X,Y,h=h,a=0,b=1,proper=TRUE)
  mass <- integrate(density,0,1,rel.tol=1e-11)$value
  expect_equal(raw$F,oracle,tolerance=2e-8)
  expect_equal(raw$f,density(Y),tolerance=1e-14)
  expect_equal(proper$F,oracle/mass,tolerance=2e-8)
  expect_equal(proper$f,density(Y)/mass,tolerance=2e-8)
  expect_identical(raw$F[2],raw$F[4])
  single <- npuniden.boundary(X,Y=.25,h=h,a=0,b=1)
  expect_identical(single$F,unname(raw$F[2]))
  expect_equal(raw$sd.F,sqrt(abs(oracle*(1-oracle)/length(X))),tolerance=2e-8)
  unlimited <- npuniden.boundary(X,Y,h=h,a=-Inf,b=Inf)
  expect_equal(unlimited$F,vapply(Y,function(y)mean(pnorm((y-X)/h)),0.0),tolerance=0)
})

test_that("proper floating-boundary density does not acquire an infinite constant tail", {
  X <- c(.95,1,1.05);Y <- c(0,.2,.5,1,2,10,100,1000)
  fit <- npuniden.boundary(X,Y,h=1,a=0,kertype="fbl",proper=TRUE)
  expect_identical(unname(fit$f[6:8]),c(0,0,0))
  expect_identical(unname(fit$f[1]),0)
  expect_gte(min(fit$f),0)
  expect_gte(min(diff(fit$F)),-1e-14)
  expect_equal(fit$F[6:8],c(1,1,1),tolerance=2e-8)
  expect_identical(unname(fit$F[1]),0)
})

test_that("shape CDF integration retains the fitted QP coefficients and density scale", {
  X <- seq(.05,.95,length.out=12)^2;Y <- c(.13,.43,.83);h <- .3
  # Capture a copy of the actual QP solution without replacing its owner.
  captured <- new.env(parent=emptyenv())
  original <- npuniden.sc
  scope <- new.env(parent=environment(original))
  scope$solve.QP <- function(...) {
    out <- quadprog::solve.QP(...)
    captured$solution <- out$solution
    out
  }
  environment(original) <- scope
  for(bounds in list(c(0,1),c(-Inf,Inf),c(0,Inf),c(-Inf,1)))
    for(shape in c("mono.incr","log-concave")) for(normalize in c(FALSE,TRUE)) {
    a <- bounds[1];b <- bounds[2]
    fit <- original(X,Y,h=h,a=a,b=b,constraint=shape,integral.equal=normalize)
    expect_true(fit$solve.QP)
    solution <- captured$solution
    kernel <- function(z) dnorm((z-X)/h)/(h*(pnorm((b-z)/h)-pnorm((a-z)/h)))
    corr <- if(normalize) fit$f.sc.integral/fit$f.integral else 1
    f <- function(y) vapply(y,function(z)mean(kernel(z)),0.0)
    fs <- function(y) vapply(y,function(z) {
      k <- kernel(z); correction <- sum(k*solution)
      (if(shape=="log-concave") exp(log(mean(k))+correction) else mean(k)+correction)/corr
    },0.0)
    expect_equal(fit$f.sc,fs(Y),tolerance=1e-12)
    expect_equal(fit$F,vapply(Y,function(y)integrate(f,a,y,rel.tol=1e-11)$value,0.0),
                 tolerance=2e-8)
    expect_equal(fit$F.sc,vapply(Y,function(y)integrate(fs,a,y,rel.tol=1e-11)$value,0.0),
                 tolerance=2e-8)
  }
})

test_that("each boundary kernel uses the same continuous CDF contract", {
  # Keep the actual unchanged kernel, but replace neither its values nor the
  # integration owner. QUADPACK below is independent of the package integrator.
  owner <- npuniden.boundary
  parts <- as.list(body(owner))
  location <- which(vapply(parts,function(z) is.call(z) &&
    identical(z[[1]],as.name("<-")) && identical(z[[2]],as.name("int.kernel.squared")),FALSE))
  expect_length(location,1)
  body(owner) <- as.call(c(parts[seq_len(location-1L)],list(quote(return(environment())))))
  X <- seq(.03,.96,length.out=20)^1.3;Y <- c(.13,.31,.72,.91);h <- .18
  for(type in c("gaussian1","gaussian2","beta1","beta2","fb","fbl","fbu","gamma","rigaussian")) {
    args <- list(X=X,Y=Y,h=h,a=0,b=1,kertype=type)
    state <- do.call(owner,args);a <- state$a;b <- state$b
    density <- function(y) vapply(y,function(z)mean(state$kernel(z,X,h,a,b)),0.0)
    splits <- sort(unique(c(a,b,.5,a+h,b-h,X-h,X,X+h)))
    splits <- splits[!is.na(splits) & splits>=a & splits<=b]
    integral <- function(upper,positive) {
      knots <- sort(unique(c(a,upper,splits[splits>a & splits<upper])))
      f <- if(positive)function(y)pmax(density(y),0) else density
      sum(vapply(seq_len(length(knots)-1L),function(i)
        integrate(f,knots[i],knots[i+1L],rel.tol=1e-10,abs.tol=1e-11,
                  subdivisions=500L)$value,0.0))
    }
    for(proper in c(FALSE,TRUE)) {
      fit <- do.call(npuniden.boundary,c(args,list(proper=proper)))
      mass <- if(proper)integral(b,TRUE) else 1
      expect_equal(fit$F,vapply(Y,integral,0.0,positive=proper)/mass,tolerance=2e-8,
                   info=paste(type,proper))
      expect_equal(fit$f,(if(proper)pmax(density(Y),0) else density(Y))/mass,
                   tolerance=2e-8,info=paste(type,proper))
    }
  }
  expect_error(npuniden.boundary(c(0,0),Y=.2,h=.1,a=0,kertype="gamma",proper=TRUE),
               "finite positive whole-support mass")
  set.seed(825);seed <- .Random.seed
  invisible(npuniden.boundary(X,Y,h=.1,a=0,b=1,proper=TRUE))
  expect_identical(seed,.Random.seed)
})
