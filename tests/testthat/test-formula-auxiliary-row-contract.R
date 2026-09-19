test_that("kernel-sum weights follow the formula sample once", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  d <- data.frame(x=seq(-1,1,length.out=18),y=sin(seq_len(18)))
  w <- cbind(seq_len(18)/18, seq_len(18)^2/18^2)
  e <- data.frame(x=c(-.5,.5))
  for (idx in list(rev(seq_len(18)), c(2L,4L,6L,8L,10L,12L))) {
    a <- npksum(y~x,data=d,subset=idx,newdata=e,bws=.4,weights=w)
    b <- npksum(txdat=d[idx,"x",drop=FALSE],tydat=d$y[idx],exdat=e,bws=.4,weights=w[idx,,drop=FALSE])
    expect_equal(a$ksum,b$ksum)
  }
  d$x[4] <- NA
  a <- npksum(y~x,data=d,newdata=e,bws=.4,weights=w)
  b <- npksum(txdat=d[-4,"x",drop=FALSE],tydat=d$y[-4],exdat=e,bws=.4,weights=w[-4,,drop=FALSE])
  expect_equal(a$ksum,b$ksum)
  # Time-indexed auxiliaries intersect with the formula, before selection.
  y <- ts(seq_len(18)); weights <- ts(matrix(seq_len(18)^2,ncol=1))
  keep <- rev(seq_len(17))
  a <- npksum(y~lag(y,-1),subset=keep,bws=.4,weights=weights)
  aligned <- ts.intersect(y,lag(y,-1),weights)
  b <- npksum(txdat=data.frame(x=aligned[keep,2]),tydat=aligned[keep,1],
               weights=matrix(aligned[keep,3],ncol=1),bws=.4)
  expect_equal(as.numeric(a$ksum),as.numeric(b$ksum))
})

test_that("location-scale formula auxiliaries follow subset and NA exclusion", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  d <- data.frame(x=seq(-1,1,length.out=18),y=sin(seq_len(18)))
  sig <- seq(.2,1,length.out=18)
  p <- rev(seq_len(18))
  args <- list(tau=.75,delta=.7,bw=.4,scale=sig,regtype="ll",bandwidth.compute=FALSE)
  a <- do.call(nplsqregbw,c(list(bws=y~x,data=d,subset=p),args))
  b <- nplsqregbw(xdat=d[p,"x",drop=FALSE],ydat=d$y[p],bws=.4,scale=sig[p],
                 tau=.75,delta=.7,regtype="ll",bandwidth.compute=FALSE)
  expect_equal(a$qdat,b$qdat)
  d$x[4] <- NA
  a <- nplsqreg(bws=y~x,data=d,scale=sig,bw=.4,bandwidth.compute=FALSE,
                 tau=.75,delta=.7,regtype="ll",na.action=na.exclude)
  b <- nplsqreg(txdat=d[-4,"x",drop=FALSE],tydat=d$y[-4],scale=sig[-4],
                 bws=.4,bandwidth.compute=FALSE,tau=.75,delta=.7,regtype="ll")
  expect_equal(as.numeric(fitted(a)[-4]),as.numeric(fitted(b)))
  expect_true(is.na(fitted(a)[4]))
})
