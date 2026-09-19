test_that("LSQ and kernel evaluation share scalar numeric formula preparation", {
  old <- options(np.messages=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(190925)
  d <- data.frame(x=runif(40),y=rnorm(40))
  square <- function(x) x^2
  for (tau in list(.3,c(.7,.3))) {
    a <- nplsqreg(bws=y~I(x^2), data=d,scale=rep(1,40),bw=.5,delta=.4,
                  tau=tau,bandwidth.compute=FALSE)
    b <- nplsqreg(bws=y~square(x), data=d,scale=rep(1,40),bw=.5,delta=.4,
                  tau=tau,bandwidth.compute=FALSE)
    for (se.fit in c(FALSE,TRUE)) {
      pa <- predict(a,newdata=d[1:5,],se.fit=se.fit)
      pb <- predict(b,newdata=d[1:5,],se.fit=se.fit)
      expect_equal(pa,pb,tolerance=1e-12)
    }
    ae <- nplsqreg(bws=y~I(x^2),data=d,newdata=d[1:5,],scale=rep(1,40),
                   bw=.5,delta=.4,tau=tau,bandwidth.compute=FALSE)
    expect_equal(fitted(ae),predict(b,newdata=d[1:5,]),tolerance=1e-12)
  }
  a <- npksum(y~I(x^2),data=d,newdata=d[1:5,],bws=.5)
  b <- npksum(y~square(x),data=d,newdata=d[1:5,],bws=.5)
  expect_equal(a$ksum,b$ksum,tolerance=1e-12)
  calls <- new.env(parent=emptyenv()); calls$n <- 0L
  once <- function(x) {calls$n <- calls$n+1L; I(x^2)}
  f <- nplsqreg(bws=y~once(x),data=d,scale=rep(1,40),bw=.5,delta=.4,
                tau=c(.7,.3),bandwidth.compute=FALSE)
  calls$n <- 0L
  predict(f,newdata=d[1:5,])
  expect_identical(calls$n,1L)
})

