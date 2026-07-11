test_that("bounded integrated kernels have exact endpoints and uniform large-h limits", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  train <- data.frame(y=rep(0.37,4L))
  eval <- data.frame(y=c(0,0.25,0.5,0.75,1))
  specs <- rbind(
    data.frame(kernel="gaussian",order=c(2L,4L,6L,8L)),
    data.frame(kernel="epanechnikov",order=c(2L,4L,6L,8L)),
    data.frame(kernel="uniform",order=2L),
    data.frame(kernel="truncated gaussian",order=2L)
  )
  for(i in seq_len(nrow(specs))){
    args <- list(dat=train,bws=1e16,bandwidth.compute=FALSE,
                 ckertype=as.character(specs$kernel[i]),
                 ckerbound="fixed",ckerlb=0,ckerub=1)
    if(args$ckertype != "uniform") args$ckerorder <- specs$order[i]
    bw <- suppressWarnings(do.call(npudistbw,args))
    got <- npudist(tdat=train,edat=eval,bws=bw)$dist
    expect_identical(got[c(1L,5L)],c(0,1))
    expect_equal(got,eval$y,tolerance=2e-13,
                 info=paste(specs$kernel[i],specs$order[i]))
  }
})

test_that("bounded Gaussian CDF repairs ordinary-bandwidth endpoint inflation", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(42)
  y <- runif(250); a <- min(y); b <- max(y)
  bw <- npudistbw(dat=data.frame(y=y),bws=.1,bandwidth.compute=FALSE,
                  ckerbound="range")
  grid <- data.frame(y=seq(a,b,length.out=101L))
  got <- npudist(tdat=data.frame(y=y),edat=grid,bws=bw)$dist
  expect_identical(got[c(1L,101L)],c(0,1))
  expect_true(all(is.finite(got)))
  expect_true(all(diff(got)>=-2e-14))
  expect_true(all(got>=-2e-14 & got<=1+2e-14))
})

test_that("bounded distribution endpoints propagate across bandwidth modes", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(123)
  x <- runif(80); y <- runif(80); a <- min(y); b <- max(y)
  ey <- data.frame(y=c(a,(a+b)/2,b))
  for(type in c("fixed","generalized_nn","adaptive_nn")){
    ubws <- if(type=="fixed") .2 else 7L
    ubw <- npudistbw(dat=data.frame(y=y),bws=ubws,bwtype=type,
                     bandwidth.compute=FALSE,ckerbound="range")
    ud <- npudist(tdat=data.frame(y=y),edat=ey,bws=ubw)$dist
    expect_identical(ud[c(1L,3L)],c(0,1),info=type)

    cbws <- if(type=="fixed") c(.2,.2) else c(7L,7L)
    cbw <- npcdistbw(xdat=data.frame(x=x),ydat=data.frame(y=y),bws=cbws,
                     bwtype=type,bandwidth.compute=FALSE,
                     cxkerbound="range",cykerbound="range",regtype="lc")
    cd <- npcdist(bws=cbw,exdat=data.frame(x=rep(.5,3L)),eydat=ey)$condist
    expect_equal(cd[c(1L,3L)],c(0,1),tolerance=2e-13,info=type)
    expect_true(all(is.finite(cd)))
  }
})
