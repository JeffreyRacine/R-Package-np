test_that("pure-ordered CDF SE uses joint contribution variance", {
  old <- options(np.messages=FALSE,np.largelambda=FALSE,
                 np.tree=FALSE,np.categorical.compress=FALSE)
  on.exit(options(old),add=TRUE)
  contribution <- function(tr,ev,h,kind,nlevels) {
    if(kind=="wangvanryzin") {
      # Independent PMF addition around the symmetry anchor. The assembled
      # W1 producer corrects the formerly decreasing upper-tail contribution.
      return(vapply(ev,function(threshold) vapply(tr,function(center) {
        pmf <- function(z) ifelse(z==center,1-h,
                                 .5*(1-h)*h^abs(z-center))
        anchor <- .5+.5*pmf(center)
        if(threshold==center) return(anchor)
        if(threshold>center)
          return(anchor+sum(pmf(seq.int(center+1L,threshold))))
        anchor-sum(pmf(seq.int(threshold+1L,center)))
      },0.),numeric(length(tr))))
    }
    vapply(ev,function(e) vapply(tr,function(t) {
      support <- if(kind=="racineliyan") seq_len(nlevels) else (-120L):120L
      mass <- h^abs(support-t)
      if(kind=="racineliyan") sum(mass[support<=e])/sum(mass)
      else sum(mass[support<=e])*(1-h)/(1+h)
    },0.),numeric(length(tr)))
  }
  for(tree in c(FALSE,TRUE))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
  for(kind in c("wangvanryzin","liracine","racineliyan")) {
    options(np.tree=tree,np.categorical.compress=tree)
    dat <- data.frame(
      a=ordered(rep(c(1,1,1,2,2,3,1),length.out=257L),levels=1:3),
      b=ordered(rep(c(1,1,2,1,2,3,2),length.out=257L),levels=1:3))
    edat <- dat[rep(c(1L,3L,4L,6L),3L),,drop=FALSE]
    h <- c(.2,.35)
    z <- contribution(as.integer(dat$a),as.integer(edat$a),h[1],kind,3)*
         contribution(as.integer(dat$b),as.integer(edat$b),h[2],kind,3)
    target <- vapply(seq_len(ncol(z)),function(j) {
      v <- z[,j]-z[1,j]
      sqrt(mean((v-mean(v))^2)/nrow(dat))
    },0.)
    b <- npudistbw(dat=dat,bws=h,bandwidth.compute=FALSE,bwtype=type,okertype=kind)
    f <- npudist(bws=b,tdat=dat,edat=edat)
    expect_equal(as.numeric(fitted(f)),colMeans(z),tolerance=1e-12)
    expect_equal(as.numeric(se(f)),target,tolerance=1e-12)
    expect_equal(as.numeric(predict(f,edat=edat)),as.numeric(fitted(f)),tolerance=0)
  }
  for(tree in c(FALSE,TRUE))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
  for(h in c(0,1)) {
    options(np.tree=tree,np.categorical.compress=tree)
    dat <- data.frame(a=ordered(rep(c(1,1,2,3),length.out=257L),levels=1:3))
    edat <- dat[rep(1:4,3L),,drop=FALSE]
    b <- npudistbw(dat=dat,bws=h,bandwidth.compute=FALSE,bwtype=type,okertype="racineliyan")
    f <- npudist(bws=b,tdat=dat,edat=edat)
    if(h==1) expect_true(all(se(f)==0))
    else expect_identical(as.numeric(se(f)),
      sqrt(as.numeric(fitted(f))*(1-as.numeric(fitted(f)))/nrow(dat)))
  }
})
