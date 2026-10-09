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
