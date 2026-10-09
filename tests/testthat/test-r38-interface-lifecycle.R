test_that("single-index pass policies expose the actual evaluation row map", {
  withr::local_options(np.messages=FALSE, na.action=getOption("na.action"))
  set.seed(3850);d<-data.frame(x=runif(120),z=rnorm(120));d$y<-sin(d$x+.5*d$z)+rnorm(120,sd=.3)
  b<-npindexbw(y~x+z,data=d,bws=c(1,.5,.4),bandwidth.compute=FALSE)
  e<-d[1:25,];e$x[c(2,7)]<-NA;e$y[c(3,9)]<-NA
  for(ambient in c("na.omit","na.exclude")) for(action in list(na.pass,NULL)) {
    options(na.action=ambient)
    actual<-npindex(b,newdata=e,se=TRUE,gradients=TRUE,na.action=action)
    expected<-npindex(b,newdata=e[-c(2,7),],se=TRUE,gradients=TRUE,na.action=na.fail)
    expect_equal(actual$rows.omit,c(2L,7L));expect_equal(actual$nobs.omit,2L)
    expect_s3_class(actual$omit,"omit");expect_length(fitted(actual),23L)
    expect_equal(fitted(actual),fitted(expected),tolerance=0)
    expect_equal(se(actual),se(expected),tolerance=0)
    expect_equal(gradients(actual),gradients(expected),tolerance=0)
    expect_equal(actual$diagnostics.nobs,21L)
  }
})

test_that("formula index calls retain user scope without a generic activation", {
  withr::local_options(np.messages=FALSE)
  set.seed(3813);d<-data.frame(x=runif(60),z=runif(60));d$y<-d$x+.5*d$z
  fit<-npindex(y~x+z,data=d,bws=c(1,.5,.2),bandwidth.compute=FALSE,se=FALSE)
  expect_identical(environment(fit$call),environment(fit$bws$formula))
  expect_identical(fit$ynames,"y");expect_identical(fit$bws$ynames,"y")
  expect_identical(fit$bws$varnames$y,"y")
  make<-function(){h<-function(x)x^2;npindex(y~h(x)+z,data=d,bws=c(1,.5,.2),bandwidth.compute=FALSE,se=FALSE)}
  local.fit<-unserialize(serialize(make(),NULL))
  expect_true(exists("h",envir=environment(local.fit$call),inherits=FALSE))
  expect_equal(predict(local.fit,newdata=d[1:5,]),predict(local.fit,exdat=data.frame('h(x)'=d$x[1:5]^2,z=d$z[1:5],check.names=FALSE)),tolerance=0)
})

test_that("native vector names support named newdata without renaming frames", {
  withr::local_options(np.messages=FALSE)
  set.seed(387);x<-runif(60);y<-sin(x)+rnorm(60,sd=.1)
  for(wrapped in c(FALSE,TRUE)) {
    dens<-if(wrapped)function(...)npudens(...)else npudens
    reg<-if(wrapped)function(...)npreg(...)else npreg
    cdens<-if(wrapped)function(...)npcdens(...)else npcdens
    f<-dens(tdat=x,bws=.3);expect_identical(f$bws$xnames,"x")
    expect_equal(predict(f,newdata=data.frame(x=.5)),predict(f,edat=.5),tolerance=0)
    f<-reg(txdat=x,tydat=y,bws=.3);expect_identical(f$bws$xnames,"x");expect_identical(f$bws$ynames,"y")
    expect_equal(predict(f,newdata=data.frame(x=.5)),predict(f,exdat=.5),tolerance=0)
    f<-cdens(txdat=x,tydat=y,bws=c(.3,.3));expect_identical(f$bws$xnames,"x");expect_identical(f$bws$ynames,"y")
    expect_equal(predict(f,newdata=data.frame(x=.5,y=.4)),predict(f,exdat=.5,eydat=.4),tolerance=0)
  }
  # Column names from an explicit frame/matrix have precedence over its expression.
  f<-npreg(txdat=data.frame(txdat=x),tydat=y,bws=.3);expect_identical(f$bws$xnames,"txdat")
  f<-npreg(txdat=cbind(named=x),tydat=y,bws=.3);expect_identical(f$bws$xnames,"named")
})

test_that("forwarded CMS subsets retain their caller-local values", {
  withr::local_options(np.messages=FALSE)
  set.seed(1);d<-data.frame(x=runif(70));d$y<-sin(d$x)+rnorm(70,sd=.3)
  selected<-d[d$x>.4,];m<-lm(y~x,data=selected,x=TRUE,y=TRUE)
  wrap<-function(...)npcmstest(y~x,data=d,model=m,B=9,bws=.3,bandwidth.compute=FALSE,...)
  actual<-(function(){keep<-d$x>.4;wrap(subset=keep)})()
  expected<-npcmstest(xdat=selected["x"],ydat=selected$y,model=m,B=9,bws=.3,bandwidth.compute=FALSE)
  expect_equal(actual[c("Jn","In","P")],expected[c("Jn","In","P")],tolerance=0)
  qm<-quantreg::rq(y~x,data=selected,tau=.5,model=TRUE)
  wrap<-function(...)npqcmstest(y~x,data=d,model=qm,B=9,bws=.3,bandwidth.compute=FALSE,...)
  actual<-(function(){keep<-d$x>.4;wrap(subset=keep)})()
  expected<-npqcmstest(xdat=selected["x"],ydat=selected$y,model=qm,B=9,bws=.3,bandwidth.compute=FALSE)
  expect_equal(actual[c("Jn","In","P")],expected[c("Jn","In","P")],tolerance=0)
})

test_that("native descriptions survive inactive forwarding closures without forcing", {
  withr::local_options(np.messages=FALSE)
  set.seed(3889);x<-runif(50);y<-sin(x)
  maker<-function(...)function()npreg(...)
  f<-maker(txdat=x,tydat=y,bws=.3)()
  expect_identical(f$bws$xnames,"x");expect_identical(f$bws$ynames,"y")
  expect_identical(fitted(f),fitted(npreg(txdat=x,tydat=y,bws=.3)))
  hold<-function(...)environment();n<-0L;e<-hold({n<-n+1L;x})
  expect_identical(.np_formula_dot_expression(quote(..1),e),quote({n<-n+1L;x}))
  expect_identical(n,0L)
})

test_that("composite native vector labels agree with their children and predictions", {
  withr::local_options(np.messages=FALSE)
  set.seed(38292);x<-runif(60);z<-runif(60);y<-sin(x)+z+rnorm(60,sd=.2)
  for(wrapped in c(FALSE,TRUE)) {
    pl<-if(wrapped)function(...)npplreg(...)else npplreg
    lsq<-if(wrapped)function(...)nplsqreg(...)else nplsqreg
    f<-pl(txdat=x,tydat=y,tzdat=z,bws=matrix(.3,2,1),bandwidth.compute=FALSE)
    expect_identical(f$bws$xnames,"x");expect_identical(f$bws$ynames,"y");expect_identical(f$bws$znames,"z")
    expect_identical(f$bws$bw$yzbw$ynames,"y")
    expect_identical(f$bws$bw[[2L]]$ynames,"x")
    expect_true(all(vapply(f$bws$bw,function(b)identical(b$xnames,"z"),logical(1))))
    expect_equal(predict(f,newdata=data.frame(x=.5,z=.5)),predict(f,exdat=.5,ezdat=.5),tolerance=0)
    for(tau in list(.5,c(.25,.5))) {
      f<-lsq(txdat=x,tydat=y,bws=.3,tau=tau,delta=.5,scale=rep(1,60),bandwidth.compute=FALSE)
      expect_identical(f$bws$xnames,"x");expect_identical(f$bws$ynames,"y")
      children<-if(is.null(f$bws$tau.bws))list(f$bws)else f$bws$tau.bws
      expect_true(all(vapply(children,function(b)identical(b$reg.bws$xnames,"x")&&identical(b$reg.bws$ynames,"y"),logical(1))))
      expect_equal(predict(f,newdata=data.frame(x=.5)),predict(f,exdat=.5),tolerance=0)
      expect_equal(fitted(lsq(bws=f$bws,txdat=x,tydat=y,tau=tau)),fitted(f),tolerance=0)
    }
  }
})
