test_that("quantile stopping is origin-free and scales in response units", {
  inv <- getFromNamespace(".npqreg_invert_selected_cdf", "npRmpi")
  b <- list(regtype.engine="lc", basis.engine="glp", degree.engine=0L,
            bernstein.basis.engine=FALSE, xncon=1L)
  x <- data.frame(x=0)
  solve <- function(shift, scale, tau=.61) {
    calls <- 0L
    cdf <- function(bws,xdat,ydat,exdat,ycand) {
      calls <<- calls+1L
      ((ycand-shift)/scale+1)/2
    }
    value <- inv(b,x,data.frame(y=shift+scale*c(-1,1)),x,tau=tau,
                  tol=1e-6,small=scale*1e-9,itmax=100L,cdf.values=cdf)
    list(value=(as.numeric(value)-shift)/scale, calls=calls)
  }
  a <- solve(0,1)
  expect_lt(abs(a$value-.22),1.001e-6)
  for (scale in c(1/8,1,8)) for (shift in c(-16384,0,16384)) {
    z <- solve(shift,scale)
    expect_identical(z$calls,a$calls)
    expect_equal(z$value,a$value,tolerance=0)
  }
})

test_that("quantile refinement terminates at adjacent finite endpoints", {
  inv <- getFromNamespace(".npqreg_invert_selected_cdf", "npRmpi")
  b <- list(regtype.engine="lc", basis.engine="glp", degree.engine=0L,
            bernstein.basis.engine=FALSE, xncon=1L)
  x <- data.frame(x=0)
  for (ends in list(c(1,1+.Machine$double.eps),
                    c(2^1023,2^1023*(1+.Machine$double.eps)),
                    c(-1.6e308,1.6e308))) {
    calls <- 0L
    cdf <- function(bws,xdat,ydat,exdat,ycand) {
      calls <<- calls+1L
      if (max(abs(ends)) > 1e308) (ycand/1.6e308+1)/2 else
        as.numeric(ycand >= ends[2L])
    }
    z <- inv(b,x,data.frame(y=ends),x,tau=.6,tol=1e-6,
              small=.Machine$double.xmin,itmax=100L,cdf.values=cdf)
    expect_true(is.finite(z))
    expect_true(z>=ends[1L] && z<=ends[2L])
    if (ends[1L]>0) expect_identical(calls,2L) else
      expect_lt(abs(as.numeric(z)/1.6e308-.2),2e-6)
  }
})

test_that("translated public quantiles agree with a Gaussian mixture root", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  x <- seq(-1,1,length.out=17L)
  y <- round(8*(x+.3*cos(seq_along(x))))/8
  ex <- c(-.625,.25,.75); taus <- c(.25,.6)
  roots <- vapply(taus,function(tau) vapply(ex,function(at) {
    w <- dnorm((at-x)/.6); w <- w/sum(w)
    f <- function(q) sum(w*pnorm((q-y)/.4))-tau
    if(f(min(y))>=0) min(y) else if(f(max(y))<0) max(y) else
      uniroot(f,range(y),tol=1e-12)$root
  },numeric(1L)),numeric(length(ex)))
  for (shift in c(0,16384)) {
    d <- data.frame(x=x,y=y+shift)
    b <- npcdistbw(y~x,data=d,bws=c(.4,.6),bandwidth.compute=FALSE)
    z <- npqreg(b,exdat=data.frame(x=ex),tau=taus,
                tol=1e-6,small=1e-9,gradients=TRUE,se=TRUE)
    expect_equal(unname(z$quantile)-shift,roots,tolerance=3e-6)
    expect_true(all(is.finite(z$quantgrad)))
    expect_true(all(is.finite(z$quanterr)))
    expect_identical(predict(z,exdat=data.frame(x=ex)),fitted(z))
    if(shift==0) original <- z else {
      expect_equal(unname(z$quantile)-shift,unname(original$quantile),tolerance=1e-10)
      expect_equal(z$quantgrad,original$quantgrad,tolerance=1e-9)
      expect_equal(z$quanterr,original$quanterr,tolerance=1e-9)
    }
  }
})
