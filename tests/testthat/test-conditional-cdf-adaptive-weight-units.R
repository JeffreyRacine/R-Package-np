for (kernel in c("gaussian", "epanechnikov", "uniform", "bounded", "beta")) {
  test_that(paste("conditional CDF CV uses absolute adaptive X weights:", kernel), {
    skip_on_cran()
    old <- options(np.messages=FALSE, np.tree=FALSE, np.extendednn=FALSE)
    on.exit(options(old), add=TRUE)
    x1 <- rep(c(.1,.3,.7), each=3L)
    y <- data.frame(y=ordered(c(0,0,1, 0,1,1, 1,2,2), levels=0:2))
    grid <- data.frame(y=ordered(c(2,0,1,0), levels=0:2))
    n <- length(x1)
    for (dimensions in 1:2) {
      x <- data.frame(x1=x1)
      if(dimensions==2L) x$x2 <- 1-x1
      # k=4 and its successor have identical positive radii in every donor
      # row. Thus deleting any one other occurrence leaves all radii intact;
      # this fixture isolates weight units from the separate fold repair.
      h <- vapply(seq_len(n), function(j) {
        distances <- sort(abs(x1[-j]-x1[j]))
        stopifnot(distances[4L] > 0, distances[4L] == distances[5L])
        distances[4L]
      }, 0)
      args <- list(xdat=x, ydat=y, bws=c(.3,rep(4,dimensions)),
        bwtype="adaptive_nn", regtype="lc", bandwidth.compute=FALSE,
        cxkertype=if(kernel=="bounded")"gaussian" else kernel,
        oykertype="wangvanryzin")
      if(kernel %in% c("bounded","beta"))
        args <- c(args,list(cxkerbound="fixed",cxkerlb=rep(0,dimensions),
                            cxkerub=rep(1,dimensions)))
      b <- do.call(npcdistbw,args)
      total <- 0
      for(i in seq_len(n)) {
        fit <- fitted(npcdist(bws=b,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
          exdat=x[rep(i,nrow(grid)),,drop=FALSE],eydat=grid))
        if(kernel=="gaussian") {
          weight <- rep(1,n-1L)
          for(l in seq_len(dimensions))
            weight <- weight*dnorm((x[i,l]-x[-i,l])/h[-i])/h[-i]
          expected <- vapply(as.integer(grid$y),function(g) {
            delta <- g-as.integer(y$y[-i])
            ky <- ifelse(delta<0,.5*.3^(-delta),1-.5*.3^(delta+1))
            sum(weight*ky)/sum(weight)
          },0)
          expect_equal(as.numeric(fit),expected,tolerance=2e-13)
        }
        target <- as.integer(y$y[i]) <= as.integer(grid$y)
        total <- total+mean((target-fit)^2)/n
      }
      for(tree in list(FALSE,TRUE,"auto")) {
        options(np.tree=tree)
        actual <- .npcdistbw_eval_only(x,y,gydat=grid,bws=b)$objective
        expect_equal(as.numeric(actual),as.numeric(total),tolerance=2e-12)
      }
      options(np.tree=FALSE)
    }
  })
}
