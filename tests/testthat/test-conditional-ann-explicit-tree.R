test_that("explicit conditional ANN X trees preserve exact fold objectives", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  ns <- environment(npcdensbw)
  for(n in c(63L,64L,65L)) for(family in c("cv.ml","cv.ls","cdf")) {
    set.seed(926L+n)
    x <- data.frame(x=runif(n,-1,1))
    y <- data.frame(y=sin(x$x)+rnorm(n))
    constructor <- if(family=="cdf") get("npcdistbw",ns) else get("npcdensbw",ns)
    evaluator <- get(if(family=="cdf") ".npcdistbw_eval_only" else ".npcdensbw_eval_only",ns)
    b <- constructor(xdat=x,ydat=y,bws=rep(floor(.75*n),2),
      bwtype="adaptive_nn",bwmethod=if(family=="cv.ml") "cv.ml" else "cv.ls",
      regtype="lp",degree=2L,bernstein.basis=TRUE,cxkertype="epanechnikov",
      bandwidth.compute=FALSE)
    value <- vapply(list(FALSE,TRUE,"auto"),function(tree) {
      options(np.tree=tree)
      args <- list(xdat=x,ydat=y,bws=b)
      if(family=="cdf") args$do.full.integral <- TRUE
      do.call(evaluator,args)$objective
    },0.0)
    expect_true(all(is.finite(value)),info=paste(n,family))
    expect_true(max(abs(value-value[1]))<2e-10,info=paste(n,family))
    expect_identical(value[1],value[3])
  }
})

test_that("compact ANN AUTO preserves the same request as explicit trees", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(60927);n<-65L;x<-data.frame(x=runif(n,-1,1));y<-data.frame(y=sin(x$x)+rnorm(n))
  ns<-environment(npcdensbw)
  for(family in c("cv.ml","cv.ls","cdf")) {
    ctor<-if(family=="cdf")npcdistbw else npcdensbw
    b<-ctor(xdat=x,ydat=y,bws=c(40,48),bwtype="adaptive_nn",
      bwmethod=if(family=="cv.ml")"cv.ml"else"cv.ls",regtype="lp",degree=2L,
      bernstein.basis=TRUE,cxkertype="epanechnikov",cykertype="epanechnikov",
      bandwidth.compute=FALSE)
    evaluator<-get(if(family=="cdf")".npcdistbw_eval_only"else".npcdensbw_eval_only",ns)
    selector<-get(if(family=="cdf")".npcdistbw_tree_code"else".npcdensbw_tree_code",ns)
    value<-vapply(list(FALSE,TRUE,"auto"),function(tree){
      options(np.tree=tree)
      expect_identical(selector(b,2L,0L),get(if(identical(tree,FALSE))"DO_TREE_NO"else"DO_TREE_YES",ns))
      a<-list(xdat=x,ydat=y,bws=b);if(family=="cdf")a$do.full.integral<-TRUE
      do.call(evaluator,a)$objective
    },0.0)
    expect_equal(value[1],value[2],tolerance=2e-10)
    expect_identical(value[2],value[3])
  }
})

test_that("ANN tree requests are not narrowed by predictor count or basis width", {
  old <- options(np.messages=FALSE,np.tree=TRUE)
  on.exit(options(old),add=TRUE)
  ns <- environment(npcdensbw)
  set.seed(2926);n<-257L
  x<-as.data.frame(matrix(runif(n*4,-1,1),n,4));names(x)<-paste0("x",1:4)
  y<-data.frame(y=x[[1]]+rnorm(n))
  for(family in c("cv.ml","cv.ls","cdf")) {
    constructor<-if(family=="cdf")get("npcdistbw",ns)else get("npcdensbw",ns)
    b<-constructor(xdat=x,ydat=y,bws=rep(190,5),bandwidth.compute=FALSE,
      bwtype="adaptive_nn",bwmethod=if(family=="cv.ml")"cv.ml"else"cv.ls",
      regtype="lp",degree=rep(2L,4),basis="tensor",bernstein.basis=TRUE,
      cxkertype="epanechnikov")
    selector<-get(if(family=="cdf")".npcdistbw_tree_code"else".npcdensbw_tree_code",ns)
    args<-list(bws=b,ncon=5L,ncat=0L)
    if(family=="cdf")args$cv.context<-TRUE
    expect_identical(do.call(selector,args),get("DO_TREE_YES",ns))
  }
})
