test_that("bounded conditional NN CVLS integrates the deleted estimator", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,np.extendednn=TRUE)
  on.exit(options(old))
  set.seed(51924)
  n <- 13L
  x <- data.frame(x=runif(n,.02,.98))
  y0 <- data.frame(y=runif(n,.02,.98),v=runif(n,.02,.98))
  q <- 9L
  for(type in c("generalized_nn","adaptive_nn"))
    for(kernel in c("gaussian","beta"))
      for(shape in c("scalar","two","mixed")) {
        y <- if(shape=="scalar") y0[1] else if(shape=="two") y0 else
          data.frame(y=y0$y,o=ordered(rep(1:2,length.out=n)))
        grid <- do.call(expand.grid,lapply(y,function(a) if(is.factor(a))
          factor(levels(a),levels=levels(a),ordered=is.ordered(a)) else seq(0,1,length.out=q)))
        weights <- rep(1,nrow(grid))
        for(j in seq_along(y)) if(!is.factor(y[[j]]))
          weights <- weights*ifelse(grid[[j]] %in% c(0,1),.5,1)/(q-1)
        ny <- sum(!vapply(y,is.factor,logical(1)))
        bw <- npcdensbw(xdat=x,ydat=y,
          bws=c(ifelse(vapply(y,is.factor,logical(1)),.2,8),8),
          bandwidth.compute=FALSE,bwtype=type,bwmethod="cv.ls",
          regtype=if(shape=="two")"ll" else "lc",
          cxkertype=kernel,cykertype=kernel,cxkerbound="fixed",cxkerlb=0,cxkerub=1,
          cykerbound="fixed",cykerlb=rep(0,ny),cykerub=rep(1,ny),
          cvls.quadrature.grid="uniform",cvls.quadrature.points=c(q,q))
        actual <- .npcdensbw_eval_only(x,y,bw)$objective
        evaluate <- function() mean(vapply(seq_len(n),function(i) {
          fg <- fitted(npcdens(bw,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
            exdat=x[rep(i,nrow(grid)),,drop=FALSE],eydat=grid,se=FALSE))
          fo <- fitted(npcdens(bw,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
            exdat=x[i,,drop=FALSE],eydat=y[i,,drop=FALSE],se=FALSE))
          2*fo-sum(weights*fg^2)
        },numeric(1)))
        oracle <- if(exists(".npRmpi_with_local_regression",mode="function"))
          .npRmpi_with_local_regression(evaluate()) else evaluate()
        expect_equal(actual,oracle,tolerance=2e-10,info=paste(type,kernel,shape))
      }
})

test_that("bounded conditional GNN extended counts use the fold size", {
  old <- options(np.messages=FALSE,np.extendednn=TRUE,np.largeh=FALSE)
  on.exit(options(old))
  x <- data.frame(x=seq(.02,.98,length.out=9))
  y <- data.frame(y=c(.2,.3,.12,.5,.8,.62,.9,.72,.4))
  q <- 9L
  grid <- data.frame(y=seq(0,1,length.out=q))
  weights <- c(.5,rep(1,q-2L),.5)/(q-1)
  for(kernel in c("gaussian","beta")) for(k in c(7L,8L,18L)) {
    bw <- npcdensbw(xdat=x,ydat=y,bws=c(k,k),bandwidth.compute=FALSE,
      bwtype="generalized_nn",bwmethod="cv.ls",cxkertype=kernel,cykertype=kernel,
      cxkerbound="fixed",cxkerlb=0,cxkerub=1,cykerbound="fixed",cykerlb=0,cykerub=1,
      cvls.quadrature.grid="uniform",cvls.quadrature.points=c(q,q))
    actual <- .npcdensbw_eval_only(x,y,bw)$objective
    evaluate <- function() mean(vapply(seq_len(nrow(x)),function(i) {
      fg <- fitted(npcdens(bw,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
        exdat=x[rep(i,q),,drop=FALSE],eydat=grid,se=FALSE))
      fo <- fitted(npcdens(bw,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
        exdat=x[i,,drop=FALSE],eydat=y[i,,drop=FALSE],se=FALSE))
      2*fo-sum(weights*fg^2)
    },numeric(1)))
    oracle <- if(exists(".npRmpi_with_local_regression",mode="function"))
      .npRmpi_with_local_regression(evaluate()) else evaluate()
    expect_equal(actual,oracle,tolerance=2e-10)
  }
})
