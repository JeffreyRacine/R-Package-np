test_that("joint native inputs reject unequal rows before computation", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=12))
  y <- sin(x$x)
  cl <- factor(rep(letters[1:2],6))
  q <- npcdistbw(xdat=x,ydat=y,bws=c(.4,.3),bandwidth.compute=FALSE)
  c <- npcdensbw(xdat=x,ydat=cl,bws=c(.2,.3),bandwidth.compute=FALSE)
  r <- npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE)
  s <- npscoefbw(xdat=x,ydat=y,zdat=x,bws=.3,bandwidth.compute=FALSE)
  set.seed(618); seed <- .Random.seed
  for(n in c(6L,24L)) {
    xx <- x[rep(seq_len(12),length.out=n),,drop=FALSE]
    yy <- rep(y,length.out=n)
    cc <- rep(cl,length.out=n)
    expect_error(npqreg(q,txdat=x,tydat=yy),"same number of rows")
    expect_error(npconmode(c,txdat=x,tydat=cc),"same number of rows")
    expect_error(npconmode(c,txdat=x,tydat=cl,exdat=x,eydat=cc),"same number of rows")
    expect_error(npsigtest(r,xdat=x,ydat=yy,B=9),"same number of rows")
    expect_error(npscoefbw(xdat=x,ydat=y,zdat=xx,bws=.3,bandwidth.compute=FALSE),"same number of rows")
    expect_error(npscoefbw(s,xdat=x,ydat=y,zdat=xx,bandwidth.compute=FALSE),"same number of rows")
  }
  expect_identical(.Random.seed,seed)
})
