test_that("conditional beta CVLS quadrature includes exact support endpoints", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920313)
  n <- 18L
  x <- data.frame(x=runif(n,.03,.97))
  y <- data.frame(y=runif(n,.03,.97))
  for (q in c(11L,31L,101L)) {
    grid <- seq(0,1,length.out=q)
    weights <- rep(1/(q-1),q)
    weights[c(1,q)] <- weights[c(1,q)]/2
    oracle <- mean(vapply(seq_len(n),function(i) {
      wx <- dbeta(x$x[-i],1+x$x[i]/.3^2,1+(1-x$x[i])/.3^2)
      wx <- wx/sum(wx)
      f <- vapply(grid,function(u)
        sum(wx*dbeta(y$y[-i],1+u/.2^2,1+(1-u)/.2^2)),numeric(1L))
      obs <- sum(wx*dbeta(y$y[-i],1+y$y[i]/.2^2,1+(1-y$y[i])/.2^2))
      2*obs-sum(weights*f^2)
    },numeric(1L)))
    for (scale in c(1,3)) {
      xx <- data.frame(x=5*x$x+7)
      yy <- data.frame(y=scale*y$y-2)
      bw <- npcdensbw(xdat=xx,ydat=yy,bws=c(.2*scale,1.5),
        bandwidth.compute=FALSE,bwmethod="cv.ls",regtype="lc",
        cxkertype="beta",cykertype="beta",cxkerbound="fixed",cxkerlb=7,
        cxkerub=12,cykerbound="fixed",cykerlb=-2,cykerub=scale-2,
        cvls.quadrature.grid="uniform",cvls.quadrature.points=c(q,15L))
      expect_equal(.npcdensbw_eval_only(xx,yy,bw)$objective,oracle/scale,
                   tolerance=1e-11)
    }
  }
})
