# Unique NA preflight contract; broad ten-family coverage remains external.
test_that("stored native fits refuse ignored NA controls before forcing them", {
  withr::local_options(np.messages = FALSE)
  x <- seq(-1,1,length.out=25L)
  bw <- npudensbw(dat=x,bws=.3,bandwidth.compute=FALSE)
  visits <- 0L
  expect_error(npudens(bw,tdat=x,na.action={visits<-visits+1L;na.fail}),
               "na.action override is not supported")
  expect_identical(visits,0L)
  expect_error(npudens(bw,na.action={visits<-visits+1L;na.fail}),
               "na.action override is not supported")
  expect_identical(visits,0L)
  expect_length(fitted(npudens(bw,tdat=x)),25L)
})
test_that("formula replacement owns NA once and retained data keep their policy", {
  withr::local_options(np.messages=FALSE)
  d <- data.frame(x=seq(-1,1,length.out=25L));d$y<-sin(d$x)
  b <- npregbw(y~x,data=d,bws=.3,bandwidth.compute=FALSE)
  visits<-0L;policy<-function(x){visits<<-visits+1L;na.omit(x)}
  missing<-d;missing$x[4L]<-NA_real_
  a<-npreg(b,data=missing,na.action=policy)
  expect_identical(visits,1L)
  expect_equal(fitted(a),fitted(npreg(b,data=missing[-4L,])),tolerance=0)
  expect_error(npreg(b,na.action={visits<-visits+1L;na.fail}),
               "na.action override is not supported")
  expect_identical(visits,1L)
  expect_equal(fitted(npreg(b)),fitted(npreg(b,data=d)),tolerance=0)
  e<-data.frame(x=c(-.2,NA_real_,.5))
  z<-npreg(b,newdata=e,na.action=na.exclude)
  expect_length(fitted(z),3L);expect_true(is.na(fitted(z)[2L]))
})
