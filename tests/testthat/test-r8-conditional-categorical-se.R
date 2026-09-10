r8_category_factor <- function(train,query,lambda,cdf=FALSE) {
  vapply(as.integer(train),function(donor) {
    weight <- lambda^abs(donor-seq_len(nlevels(train)))
    if(cdf)sum(weight[seq_len(as.integer(query))])/sum(weight) else
      weight[as.integer(query)]/sum(weight)
  },numeric(1))
}

test_that("conditional finite categories retain leading and ratio moments", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 73L; x <- (seq_len(n)/79)^1.1; y <- .15+.45*x+sin(seq_len(n))*.09
  a <- factor(rep(c("a","b","a"),length.out=n))
  b <- ordered(ifelse(seq_len(n)%%5<3,as.character(a),"c"))
  X <- data.frame(x=x,a=a); Y <- data.frame(y=y,b=b)
  EX <- data.frame(x=c(.21,.43,.71),a=factor(c("a","b","a"),levels=levels(a)))
  EY <- data.frame(y=c(.25,.37,.49),b=ordered(c("a","b","c"),levels=levels(b)))
  rk <- 1/(2*sqrt(pi))
  for(type in c("fixed","generalized_nn")) for(kind in c("density","ratio","cdf")) {
    hX <- if(type=="fixed") .18 else 17L
    hY <- if(type=="fixed") .23 else 15L
    for(ly in c(0,.31,1)) {
      yy <- if(kind=="ratio")data.frame(b=b) else Y
      ey <- if(kind=="ratio")EY["b"] else EY
      bwfn <- if(kind=="cdf")npcdistbw else npcdensbw
      fn <- if(kind=="cdf")npcdist else npcdens
      bw <- bwfn(xdat=X,ydat=yy,
        bws=c(if(kind=="ratio")ly else c(hY,ly),hX,.22),
        bwtype=type,oykertype="racineliyan",regtype="lc",bandwidth.compute=FALSE)
      expect_equal(unname(bw$xbw),c(hX,.22))
      f <- fn(bws=bw,txdat=X,tydat=yy,exdat=EX,eydat=ey,gradients=TRUE,se=TRUE)
      radiusX <- if(type=="fixed")rep(hX,3) else
        vapply(EX$x,function(v)sort(abs(x-v))[[hX]],numeric(1))
      radiusY <- if(type=="fixed")rep(hY,3) else
        vapply(EY$y,function(v)sort(abs(y-v))[[hY]],numeric(1))
      reference <- vapply(seq_len(3),function(j) {
        KX <- dnorm((EX$x[j]-x)/radiusX[j])/radiusX[j]
        L <- ifelse(a==EX$a[j],1-.22,.22)
        B <- r8_category_factor(b,EY$b[j],ly,kind=="cdf")
        TY <- if(kind=="ratio")rep(1,n) else if(kind=="density")
          dnorm((EY$y[j]-y)/radiusY[j])/radiusY[j] else
          pnorm((EY$y[j]-y)/radiusY[j]*sqrt(2)*.7071067810)
        den <- sum(KX*L); mean <- sum(KX*L*TY*B)/den
        S <- if(kind=="density")sum(KX*TY*L^2*B^2) else
          sum(KX*L^2*((TY*B-mean)^2+TY*(1-TY)*B^2))
        variance <- S*rk/(radiusX[j]*den^2)
        if(kind=="density") variance <- variance*rk/radiusY[j]
        c(mean,variance)
      },numeric(2))
      expect_equal(as.numeric(fitted(f)),unname(reference[1,]),tolerance=1e-11)
      expect_equal(f$conderr^2,unname(reference[2,]),tolerance=1e-11)
      expect_equal(f$congerr[,1]^2,
        unname(reference[2,])/(2*radiusX^2),tolerance=1e-11)
      if(kind=="ratio" && ly==1)expect_lt(max(f$conderr),1e-14)
    }
  }
})

test_that("pure-category ratio covariance respects profile counts and NN labels", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 193L
  x <- ordered(rep(c("a","b","a","c"),length.out=n))
  y <- ordered(ifelse(seq_len(n)%%7<4,as.character(x),"b"),levels=levels(x))
  ## Keep repeated profiles within each fixed-bandwidth MPI owner block.
  ex <- rep(x[1:6],12); ey <- rep(y[4:9],12)
  for(cdf in c(FALSE,TRUE)) for(type in c("fixed","generalized_nn","adaptive_nn")) {
    bwfn <- if(cdf)npcdistbw else npcdensbw; fn <- if(cdf)npcdist else npcdens
    bw <- bwfn(xdat=x,ydat=y,bws=c(.31,.41),bwtype=type,regtype="lc",
      oxkertype="racineliyan",oykertype="racineliyan",bandwidth.compute=FALSE)
    reference <- vapply(seq_along(ex),function(j) {
      L <- r8_category_factor(x,ex[j],.41)
      B <- r8_category_factor(y,ey[j],.31,cdf)
      den <- sum(L); mean <- sum(L*B)/den
      c(mean,sum(L^2*(B-mean)^2)/den^2)
    },numeric(2))
    for(tree in c(FALSE,TRUE)) {
      options(np.tree=tree)
      f <- fn(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey,se=TRUE)
      expect_equal(as.numeric(fitted(f)),unname(reference[1,]),tolerance=1e-12)
      expect_equal(f$conderr^2,unname(reference[2,]),tolerance=1e-12)
    }
  }
})

test_that("fully smoothed X cancels from a finite response-kernel mean", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- factor(rep(c("a","b","a"),length.out=61L))
  y <- factor(rep(c("a","a","b","c","a"),length.out=61L))
  for(ly in c(0,.2,2/3)) {
    bw <- npcdensbw(xdat=x,ydat=y,bws=c(ly,.5),bandwidth.compute=FALSE)
    f <- npcdens(bws=bw,txdat=x,tydat=y,exdat=x[1:2],eydat=y[c(1,4)],se=TRUE)
    ref <- vapply(y[c(1,4)],function(v) {
      B <- ifelse(y==v,1-ly,ly/2)
      c(mean(B),mean((B-mean(B))^2)/length(B))
    },numeric(2))
    expect_equal(f$condens,unname(ref[1,]),tolerance=1e-13)
    expect_equal(f$conderr^2,unname(ref[2,]),tolerance=1e-13)
  }
})

test_that("conditional CDF bounds and empty support retain point policy", {
  old <- options(np.messages=FALSE,np.tree=TRUE)
  on.exit(options(old),add=TRUE)
  n <- 81L; x <- seq(.1,.9,length.out=n); y <- .2+.6*x
  z <- ordered(rep(letters[1:3],length.out=n))
  X <- data.frame(x=x,z=z); Y <- data.frame(y=y,z=z)
  EX <- data.frame(x=c(.3,.5,.7),z=z[1:3])
  EY <- data.frame(y=c(.25,.45,.65),z=z[1:3])
  h <- .19; lx <- .31; ly <- .41
  bw <- npcdistbw(xdat=X,ydat=Y,bws=c(h,ly,h,lx),regtype="lc",
    oxkertype="racineliyan",oykertype="racineliyan",
    cykerbound="fixed",cykerlb=-.2,cykerub=1.2,bandwidth.compute=FALSE)
  f <- npcdist(bws=bw,txdat=X,tydat=Y,exdat=EX,eydat=EY,se=TRUE)
  reference <- vapply(seq_len(3),function(j) {
    K <- dnorm((EX$x[j]-x)/h)/h
    A <- r8_category_factor(z,EX$z[j],lx)
    B <- r8_category_factor(z,EY$z[j],ly,TRUE)
    J <- (pnorm((EY$y[j]-y)/h)-pnorm((-.2-y)/h))/
      (pnorm((1.2-y)/h)-pnorm((-.2-y)/h))
    D <- sum(K*A); m <- sum(K*A*J*B)/D
    c(m,sum(K*A^2*((J*B-m)^2+J*(1-J)*B^2))/(2*sqrt(pi)*h*D^2))
  },numeric(2))
  expect_equal(f$condist,unname(reference[1,]),tolerance=1e-9)
  expect_equal(f$conderr^2,unname(reference[2,]),tolerance=1e-9)
  bw <- npcdistbw(xdat=X,ydat=Y,bws=c(h,1,h,lx),regtype="lc",
    cxkertype="epanechnikov",cykertype="uniform",
    oxkertype="racineliyan",oykertype="racineliyan",bandwidth.compute=FALSE)
  EX$x <- c(.5,.5,10); EY$y <- c(-10,10,.5); EY$z <- z[c(1,1,1)]
  f <- suppressWarnings(npcdist(bws=bw,txdat=X,tydat=Y,exdat=EX,eydat=EY,se=TRUE))
  expect_equal(f$condist[1:2],c(0,1/3),tolerance=1e-14)
  expect_lt(max(f$conderr[1:2]),1e-14)
  expect_true(is.na(f$condist[3]))
  expect_true(is.na(f$conderr[3]))
})
