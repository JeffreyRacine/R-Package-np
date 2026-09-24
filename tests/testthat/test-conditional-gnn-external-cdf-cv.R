test_that("external conditional-CDF CV keeps query and deleted donor identities distinct", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=TRUE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(-.8,-.8,-.25,0,.1,.31,.55,.82,1.05))
  y <- data.frame(y=c(-.4,-.4,-.08,.12,.21,.34,.48,.73,.9))
  for (grid in list(y[9:1,,drop=FALSE],data.frame(y=c(-.7,.05,.7)))) {
    k <- 3L;n <- nrow(x);expected <- 0
    for(i in seq_len(n)) {
      keep <- setdiff(seq_len(n),i)
      hx <- sort(abs(x$x[i]-x$x[keep]))[k]
      wx <- dnorm((x$x[i]-x$x[keep])/hx);wx <- wx/sum(wx)
      for(j in seq_len(nrow(grid))) {
        hy <- sort(abs(grid$y[j]-y$y[keep]))[k]
        pred <- sum(wx*pnorm((grid$y[j]-y$y[keep])/hy))
        expected <- expected+(as.numeric(y$y[i]<=grid$y[j])-pred)^2/(n*nrow(grid))
      }
    }
    b <- npcdistbw(xdat=x,ydat=y,bws=c(k,k),regtype="lc",
      bwtype="generalized_nn",bandwidth.compute=FALSE)
    for(tree in c(FALSE,TRUE)) {
      options(np.tree=tree)
      got <- getFromNamespace(".npcdistbw_eval_only","np")(x,y,bws=b,gydat=grid)$objective
      expect_equal(got,expected,tolerance=2e-10)
    }
  }
})

test_that("external CDF folds cover mixed roles and the high-dimensional row sibling", {
  old <- options(np.messages=FALSE,np.tree=TRUE)
  on.exit(options(old),add=TRUE)
  set.seed(9218);n <- 9L
  x <- data.frame(x=runif(n),group=factor(rep(1:3,3)))
  mixed <- data.frame(y=runif(n),rank=ordered(rep(1:3,3)))
  high <- as.data.frame(matrix(runif(n*5),n,5))
  for(label in c("mixed-lc","mixed-lp","high-lc")) {
    y <- if(label=="high-lc")high else mixed
    grid <- y[c(2,7,5),,drop=FALSE]
    reg <- if(label=="mixed-lp")"lp" else "lc"
    args <- list(xdat=x,ydat=y,bws=c(if(label=="high-lc")rep(4,5) else c(4,.2),4,.2),
      bwtype="generalized_nn",regtype=reg,bandwidth.compute=FALSE)
    if(reg=="lp")args$degree <- 2L
    bw <- do.call(npcdistbw,args)
    reference <- 0
    for(i in seq_len(n)) {
      keep <- setdiff(seq_len(n),i)
      pred <- fitted(npcdist(bws=bw,txdat=x[keep,,drop=FALSE],tydat=y[keep,,drop=FALSE],
        exdat=x[rep(i,nrow(grid)),,drop=FALSE],eydat=grid))
      indicator <- vapply(seq_len(nrow(grid)),function(j)
        all(vapply(seq_len(ncol(y)),function(k)y[i,k]<=grid[j,k],FALSE)),FALSE)
      reference <- reference+mean((indicator-pred)^2)/n
    }
    got <- getFromNamespace(".npcdistbw_eval_only","np")(x,y,gydat=grid,bws=bw,
      invalid.penalty="dbmax")$objective
    expect_equal(got,reference,tolerance=2e-10,info=label)
  }
})
