test_that("finite-category density SE uses complete contribution variance", {
  old <- options(np.messages=FALSE, np.largelambda=FALSE,
                 np.tree=FALSE, np.categorical.compress=FALSE)
  on.exit(options(old),add=TRUE)
  oracle <- function(dat,edat,h,uk,ok) {
    z <- matrix(1,nrow(dat),nrow(edat))
    for(j in seq_along(dat)) {
      tr <- as.integer(dat[[j]]); ev <- as.integer(edat[[j]])
      d <- abs(outer(tr,ev,"-")); lambda <- h[j]
      if(!is.ordered(dat[[j]])) {
        nc <- nlevels(dat[[j]])
        k <- if(uk=="aitchisonaitken")
          ifelse(d==0,1-lambda,lambda/(nc-1)) else
          ifelse(d==0,1,lambda)/(1+(nc-1)*lambda)
      } else if(ok=="wangvanryzin") {
        k <- ifelse(d==0,1-lambda,.5*(1-lambda)*lambda^d)
      } else if(ok=="liracine") {
        k <- lambda^d*(1-lambda)/(1+lambda)
      } else {
        den <- vapply(tr,function(t)
          sum(lambda^abs(t-seq_len(nlevels(dat[[j]])))),0.)
        k <- lambda^d/den
      }
      z <- z*k
    }
    list(mean=colMeans(z),se=vapply(seq_len(ncol(z)),function(j) {
      v <- z[,j]-z[1,j]
      sqrt(mean((v-mean(v))^2)/nrow(dat))
    },0.))
  }
  for(tree in c(FALSE,TRUE))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
  for(kind in c("unordered","mixed","wangvanryzin","liracine","racineliyan")) {
  options(np.tree=tree,np.categorical.compress=tree)
    n <- 257L
    a <- rep(c(1,1,1,2,2,3,1),length.out=n)
    b <- rep(c(1,1,2,1,2,3,2),length.out=n)
    dat <- data.frame(
      a=if(kind %in% c("unordered","mixed")) factor(a,levels=1:3) else ordered(a,levels=1:3),
      b=if(kind=="unordered") factor(b,levels=1:3) else ordered(b,levels=1:3))
    edat <- dat[rep(c(1L,3L,4L,6L),3L),,drop=FALSE]
    uk <- if(kind=="mixed") "liracine" else "aitchisonaitken"
    ok <- if(kind %in% c("unordered","mixed")) "racineliyan" else kind
    h <- c(.2,.35)
    bw <- npudensbw(dat=dat,bws=h,bandwidth.compute=FALSE,
                    bwtype=type,ukertype=uk,okertype=ok)
    f <- npudens(bws=bw,tdat=dat,edat=edat)
    ref <- oracle(dat,edat,h,uk,ok)
    expect_equal(as.numeric(fitted(f)),ref$mean,tolerance=1e-12)
    expect_equal(as.numeric(se(f)),ref$se,tolerance=1e-12)
  }

  for(tree in c(FALSE,TRUE))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
  for(h in c(0,.5-1e-9,.5)) {
  options(np.tree=tree,np.categorical.compress=tree)
    dat <- data.frame(a=factor(rep(c("a","a","b"),length.out=257L)))
    edat <- dat[rep(1:3,4L),,drop=FALSE]
    bw <- npudensbw(dat=dat,bws=h,bandwidth.compute=FALSE,bwtype=type)
    f <- npudens(bws=bw,tdat=dat,edat=edat)
    ref <- oracle(dat,edat,h,"aitchisonaitken","racineliyan")
    expect_equal(as.numeric(se(f)),ref$se,tolerance=1e-12)
    if(h==.5) expect_true(all(se(f)==0))
    if(h==0) expect_identical(as.numeric(se(f)),
      sqrt(as.numeric(fitted(f))*(1-as.numeric(fitted(f)))/nrow(dat)))
  }
for(tree in c(FALSE,TRUE))
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    options(np.tree=tree,np.categorical.compress=tree)
    dat <- data.frame(a=ordered(rep(c("a","a","b","c"),length.out=257L)))
    edat <- dat[rep(1:4,3L),,drop=FALSE]
    bw <- npudensbw(dat=dat,bws=1,bandwidth.compute=FALSE,
                    bwtype=type,okertype="racineliyan")
    f <- npudens(bws=bw,tdat=dat,edat=edat)
    expect_true(all(se(f)==0))
    expect_equal(as.numeric(fitted(f)),rep(1/3,nrow(edat)),tolerance=1e-13)
  }
})
