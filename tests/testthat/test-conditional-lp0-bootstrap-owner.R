test_that("LP0 conditional bootstrap uses the same count estimator as LC", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(.03,.97,length.out=24))
  y <- data.frame(y=.05+.9*((seq_len(24)*7)%%25)/25)
  ex <- data.frame(x=c(.21,.53,.76))
  ey <- data.frame(y=c(.33,.64,.49))
  counts <- cbind(rep(1,24),rep(c(0,2,1),8))
  for(cdf in c(FALSE,TRUE)) for(pair in c("gg","bb","gb","bg"))
    for(type in c("fixed","generalized_nn","adaptive_nn")) {
      kernel.x <- if(substr(pair,1,1)=="b") "beta" else "gaussian"
      kernel.y <- if(substr(pair,2,2)=="b") "beta" else "gaussian"
      maker <- if(cdf) npcdistbw else npcdensbw
      fitter <- if(cdf) npcdist else npcdens
      args <- list(xdat=x,ydat=y,bws=if(type=="fixed")c(.23,.31)else c(9,10),
        bwtype=type,bandwidth.compute=FALSE,cxkertype=kernel.x,cykertype=kernel.y)
      if(kernel.x=="beta") args <- c(args,list(cxkerbound="fixed",cxkerlb=0,cxkerub=1))
      if(kernel.y=="beta") args <- c(args,list(cykerbound="fixed",cykerlb=0,cykerub=1))
      a <- do.call(maker,c(args,list(regtype="lc")))
      b <- do.call(maker,c(args,list(regtype="lp",degree=0)))
      expect_identical(.np_con_xregtype(a),"lc")
      expect_identical(.np_con_xregtype(b),"lc")
      got <- .np_inid_boot_from_ksum_conditional(x,y,ex,ey,b,2L,cdf,counts=counts)
      lc <- .np_inid_boot_from_ksum_conditional(x,y,ex,ey,a,2L,cdf,counts=counts)
      expect_equal(got,lc,tolerance=2e-11)
      for(j in 1:2) {
        take <- rep(seq_len(24),counts[,j])
        oracle <- fitter(b,txdat=x[take,,drop=FALSE],tydat=y[take,,drop=FALSE],
                           exdat=ex,eydat=ey,se=FALSE)
        expect_equal(as.double(got$t[j,]),as.double(fitted(oracle)),tolerance=2e-11)
      }
    }
})

test_that("fixed LP0 reaches scalar operators and positive degree retains its owner", {
  x <- data.frame(x=seq(-1,1,length.out=12)); y <- data.frame(y=cos(seq_len(12)))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(.4,.4),bandwidth.compute=FALSE,regtype="lp",degree=0)
  owner <- .np_inid_boot_from_ksum_conditional
  env <- new.env(parent=environment(owner))
  env$.np_ksum_conditional_operator_fixed <- function(...) {
    state$scalar <- state$scalar+1L
    list(num=matrix(1,2,12),den=matrix(2,2,12))
  }
  env$.np_inid_boot_from_conditional_localpoly_fixed <- function(...) "positive-owner"
  environment(owner) <- env
  state <- new.env(parent=emptyenv());state$scalar <- 0L
  z <- owner(x,y,x[1:2,,drop=FALSE],y[1:2,,drop=FALSE],b,1L,FALSE,
             counts=matrix(1,12,1))
  expect_identical(state$scalar,1L)
  expect_equal(as.double(z$t),c(.5,.5))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(.4,.4),bandwidth.compute=FALSE,regtype="lp",degree=1)
  expect_identical(owner(x,y,x[1:2,,drop=FALSE],y[1:2,,drop=FALSE],b,1L,FALSE),
                   "positive-owner")
  expect_identical(state$scalar,1L)
})
