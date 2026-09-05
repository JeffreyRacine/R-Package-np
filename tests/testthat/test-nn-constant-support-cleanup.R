test_that("NN constant-support rejection leaves all bandwidth owners reusable", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 24L
  x <- data.frame(x=sin(seq_len(n)*sqrt(2))+seq_len(n)/n,
    z=cos(seq_len(n)*sqrt(3))+seq_len(n)/(2*n))
  y <- sin(seq_len(n)/3)+seq_len(n)/90
  run <- function(owner,type,side=NULL) {
    tx <- x
    ty <- y
    if(identical(side,"X")) tx$z <- 1
    if(identical(side,"Y")) ty[] <- 1
    set.seed(9L)
    if(owner=="regression")
      return(npregbw(xdat=tx,ydat=ty,bwtype=type,regtype="lc",nomad=FALSE,
        bws=c(8,8),nmulti=1L,itmax=2L,random.seed=9L))
    if(owner %in% c("density","distribution")) {
      fn <- if(owner=="density") npudensbw else npudistbw
      return(fn(dat=tx,bwtype=type,bws=c(8,8),nmulti=1L,itmax=2L))
    }
    fn <- if(owner=="conditional density") npcdensbw else npcdistbw
    fn(xdat=tx,ydat=ty,bwtype=type,regtype="lc",nomad=FALSE,
      bws=c(8,8,8),nmulti=1L,itmax=2L,random.seed=9L)
  }
  for(owner in c("regression","density","distribution",
                 "conditional density","conditional distribution"))
    for(type in c("generalized_nn","adaptive_nn")) {
      before <- run(owner,type)
      for(side in if(startsWith(owner,"conditional")) c("X","Y") else "X")
        for(i in 1:2)
          # Existing zero-scale warnings precede the native support error.
          expect_error(suppressWarnings(run(owner,type,side)),
            "nonfixed nearest-neighbour bandwidths require",fixed=TRUE)
      after <- run(owner,type)
      for(field in c("bw","xbw","ybw","fval","ifval","fval.history"))
        expect_identical(after[[field]],before[[field]])
    }
})
