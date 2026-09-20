test_that("beta regression CVLS uses literal delete-one NN geometry", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(92039)
  x <- data.frame(x=runif(21,.02,.98),c=factor(rep(letters[1:3],7)))
  y <- sin(4*x$x)+rnorm(21,sd=.1)
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    for(reg in c("lc","ll","lp")) {
      args <- list(xdat=x,ydat=y,bws=c(if(type=="fixed") .4 else 9,.2),
        bwtype=type,bandwidth.compute=FALSE,regtype=reg,ckertype="beta",
        ckerbound="fixed",ckerlb=0,ckerub=1,bwmethod="cv.ls")
      if(reg=="lp") args$degree <- 2L
      bw <- do.call(npregbw,args)
      actual <- .npregbw_eval_only(x,y,bw)$objective
      pred <- vapply(seq_len(nrow(x)),function(i)
        fitted(npreg(bws=bw,txdat=x[-i,,drop=FALSE],tydat=y[-i],
                     exdat=x[i,,drop=FALSE],se=FALSE)),numeric(1))
      expect_equal(actual,mean((y-pred)^2),tolerance=1e-10,
                   info=paste(type,reg))
    }
  }
})

test_that("beta CVLS agrees with independent beta2 distance arithmetic", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(.02,.06,.13,.23,.38,.51,.67,.8,.89,.98)
  y <- sin(5*x)+c(-.1,.1,.04,-.08,.03,.07,-.04,.11,-.08,.02)
  n <- length(x); k <- 3L
  for(type in c("generalized_nn","adaptive_nn")) {
    pred <- vapply(seq_len(n),function(i) {
      donors <- setdiff(seq_len(n),i)
      h <- if(type=="generalized_nn") rep(sort(abs(x[donors]-x[i]))[k],n-1)
      else vapply(donors,function(j)
        sort(abs(x[setdiff(seq_len(n),c(i,j))]-x[j]))[k],numeric(1))
      w <- dbeta(x[donors],1+x[i]/h^2,1+(1-x[i])/h^2)
      sum(w*y[donors])/sum(w)
    },numeric(1))
    bw <- npregbw(xdat=data.frame(x=x),ydat=y,bws=k,bwtype=type,
      ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1,
      regtype="lc",bandwidth.compute=FALSE)
    expect_equal(.npregbw_eval_only(data.frame(x=x),y,bw)$objective,
                 mean((y-pred)^2),tolerance=1e-12)
  }
})
