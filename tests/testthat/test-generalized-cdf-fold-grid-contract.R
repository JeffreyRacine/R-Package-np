.gnn_cdf_fold_objective <- function(dat, state, grid=NULL) {
  getFromNamespace("npudistbw.dbandwidth","np")(
    dat=dat,bws=state,bandwidth.compute=TRUE,eval.only=TRUE,
    do.full.integral=is.null(grid),gdat=grid,nmulti=1L,
    invalid.penalty="dbmax")$fval
}
.gnn_cdf_fold_refit <- function(dat,state,grid=NULL) {
  n <- nrow(dat)
  mean(vapply(seq_len(n),function(i) {
    args <- list(bws=state,tdat=dat[-i,,drop=FALSE])
    if(!is.null(grid))args$edat <- grid
    query <- if(is.null(grid))dat[-i,,drop=FALSE] else grid
    indicator <- rep(TRUE,nrow(query))
    for(d in seq_along(dat))indicator <- indicator & dat[[d]][i]<=query[[d]]
    mean((indicator-as.numeric(fitted(do.call(npudist,args))))^2)
  },numeric(1)))
}
test_that("GNN CDF CV uses literal deleted samples and query identities", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,np.extendednn=TRUE)
  on.exit(options(old),add=TRUE)
  set.seed(712)
  n <- 9L
  for(p in c(1L,2L,3L,5L)) for(kernel in c("gaussian","beta")) {
    dat <- as.data.frame(matrix(runif(n*p,.03,.97),ncol=p))
    dat[2L,] <- dat[1L,]
    grid <- dat[c(2L,5L,7L),,drop=FALSE]
    for(k in c(3L,10L)) for(external in c(FALSE,TRUE)) {
      a <- list(dat=dat,bws=rep(k,p),bandwidth.compute=FALSE,
        bwtype="generalized_nn",ckertype=kernel,
        ckerorder=if(k==3L) 2L else 4L)
      if(kernel=="beta")a <- c(a,list(ckerbound="fixed",ckerlb=rep(0,p),ckerub=rep(1,p)))
      b <- do.call(npudistbw,a)
      query <- if(external)grid else NULL
      expected <- .gnn_cdf_fold_refit(dat,b,query)
      for(tree in c(FALSE,TRUE)) {
        options(np.tree=tree)
        expect_equal(as.numeric(.gnn_cdf_fold_objective(dat,b,query)),expected,
          tolerance=3e-12,info=paste(p,kernel,k,external,tree))
      }
    }
  }
})
test_that("GNN CDF CV matches independent occurrence-distance beta and Gaussian sums", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(.03,.08,.15,.29,.41,.53,.69,.81,.92); n <- length(x); k <- 3L
  for(kernel in c("gaussian","beta")) {
    a <- list(dat=data.frame(x),bws=k,bwtype="generalized_nn",
      ckertype=kernel,bandwidth.compute=FALSE)
    if(kernel=="beta")a <- c(a,list(ckerbound="fixed",ckerlb=0,ckerub=1))
    b <- do.call(npudistbw,a)
    for(external in c(FALSE,TRUE)) {
      query <- if(external)c(0,.12,.46,.78,1) else x
      expected <- mean(vapply(seq_len(n),function(i) {
        js <- if(external)seq_along(query) else setdiff(seq_len(n),i)
        mean(vapply(js,function(j) {
          donors <- setdiff(seq_len(n),i)
          neighbors <- if(external)donors else setdiff(donors,j)
          h <- sort(abs(query[j]-x[neighbors]))[k]
          fit <- if(kernel=="beta")mean(pbeta(query[j],1+x[donors]/h^2,
            1+(1-x[donors])/h^2)) else mean(pnorm((query[j]-x[donors])/h))
          ((x[i]<=query[j])-fit)^2
        },numeric(1)))
      },numeric(1)))
      expect_equal(as.numeric(.gnn_cdf_fold_objective(data.frame(x),b,
        if(external)data.frame(x=query) else NULL)),expected,
        tolerance=if(kernel=="gaussian") 2e-10 else 2e-13)
    }
  }
})
test_that("GNN CDF fold rows retain ordered support and compact-kernel trees", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  dat <- data.frame(x=c(.02,.07,.15,.29,.43,.58,.69,.87,.96),
                    o=ordered(c(1,3,7,1,3,7,1,3,7),levels=c(1,3,7)))
  for(kernel in c("gaussian","epanechnikov","uniform","beta")) {
    a <- list(dat=dat,bws=c(3,.3),bwtype="generalized_nn",
      okertype="racineliyan",ckertype=kernel,bandwidth.compute=FALSE)
    if(kernel=="beta")a <- c(a,list(ckerbound="fixed",ckerlb=0,ckerub=1))
    b <- do.call(npudistbw,a)
    for(external in c(FALSE,TRUE)) {
      query <- if(external)dat[c(2,5,8),,drop=FALSE] else NULL
      expected <- .gnn_cdf_fold_refit(dat,b,query)
      for(tree in c(FALSE,TRUE)) {
        options(np.tree=tree)
        expect_equal(as.numeric(.gnn_cdf_fold_objective(dat,b,query)),expected,
          tolerance=3e-12)
      }
    }
  }
})
