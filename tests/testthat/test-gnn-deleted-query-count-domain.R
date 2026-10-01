# Heavy tree/kernel cases are retained in six singleton test files.

test_that("GNN boundary radii agree with independent Gaussian distance ranks", {
  old <- options(np.messages=FALSE, np.extendednn=TRUE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  x <- c(.03,.09,.18,.29,.43,.57,.71,.84,.97)
  d <- data.frame(x=x); y <- sin(9*x)+x; n <- length(x)
  for (k in c(n-3L,n-2L,n-1L,n+2L)) {
    kk <- min(k,n-2L)
    h <- vapply(seq_len(n), function(i) sort(abs(x[-i]-x[i]))[kk]*k/kk, 0.)
    w <- vapply(seq_len(n), function(i) dnorm((x-x[i])/h[i])/h[i], numeric(n))
    diag(w) <- 0
    expected <- colSums(w*y)/colSums(w)
    b <- npregbw(xdat=d,ydat=y,bws=k,bwtype="generalized_nn",
      regtype="lc",bandwidth.compute=FALSE)
    expect_equal(as.numeric(npreghat(b,txdat=d,y=y,output="apply",
      leave.one.out=TRUE)), expected,tolerance=2e-13)
    expect_equal(as.numeric(.npregbw_eval_only(d,y,b,
      invalid.penalty="dbmax")$objective),mean((y-expected)^2),tolerance=2e-13)
  }
})

test_that("failed beta density geometry is not reported as objective zero", {
  old <- options(np.messages=FALSE, np.extendednn=FALSE)
  on.exit(options(old), add=TRUE)
  d <- data.frame(x=c(rep(.25,6),.5,.7,.9))
  for(type in c("generalized_nn","adaptive_nn")) {
    b <- npudensbw(d,bws=2,bwtype=type,bwmethod="cv.ml",ckertype="beta",
      ckerbound="fixed",ckerlb=0,ckerub=1,bandwidth.compute=FALSE)
    expect_identical(as.numeric(npudensbw.bandwidth(d,bws=b,eval.only=TRUE,
      nmulti=1L,invalid.penalty="dbmax")$fval),-.Machine$double.xmax)
  }
})
