test_that("mixed density errors retain finite joint category factors", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 81L; x <- (seq_len(n)/83)^1.2
  a <- factor(rep(c("a","b","a"),length.out=n))
  b <- factor(ifelse(seq_len(n)%%5<3,as.character(a),"c"))
  d <- data.frame(x=x,a=a,b=b)
  e <- data.frame(x=c(.061,.417,.913),a=factor(c("a","b","a"),levels=levels(a)),
                  b=factor(c("b","c","a"),levels=levels(b)))
  rk <- 1/(2*sqrt(pi))
  for(type in c("fixed","generalized_nn")) for(lam in c(0,.22,.5)) {
    h <- if(type=="fixed") .17 else 17L
    bw <- npudensbw(dat=d,bws=c(h,lam,.28),bwtype=type,bandwidth.compute=FALSE)
    f <- npudens(bws=bw,tdat=d,edat=e,se=TRUE)
    off <- npudens(bws=bw,tdat=d,edat=e,se=FALSE)
    expect_identical(f$dens,off$dens)
    radius <- if(type=="fixed") rep(h,3) else
      vapply(e$x,function(v)sort(abs(x-v))[[h]],numeric(1))
    oracle <- vapply(seq_len(nrow(e)),function(j) {
      L <- ifelse(a==e$a[j],1-lam,lam)*
        ifelse(b==e$b[j],1-.28,.28/(nlevels(b)-1))
      K <- dnorm((e$x[j]-x)/radius[j])/radius[j]
      c(mean=mean(K*L),variance=rk*sum(K*L^2)/(n^2*radius[j]))
    },numeric(2))
    expect_equal(f$dens,unname(oracle[1,]),tolerance=1e-12)
    expect_equal(f$derr^2,unname(oracle[2,]),tolerance=1e-12)
  }
})

test_that("mixed ordered CDF errors keep the leading continuous indicator target", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 73L; x <- (seq_len(n)/79)^1.1
  z <- ordered(rep(c("a","b","c","a","a"),length.out=n))
  d <- data.frame(x=x,z=z)
  e <- data.frame(x=c(.13,.43,.87),z=ordered(c("a","b","c"),levels=levels(z)))
  for(type in c("fixed","generalized_nn")) for(lam in c(0,.37,1)) {
    h <- if(type=="fixed") .19 else 15L
    bw <- npudistbw(dat=d,bws=c(h,lam),bwtype=type,okertype="racineliyan",
                     bandwidth.compute=FALSE)
    f <- npudist(bws=bw,tdat=d,edat=e,se=TRUE)
    off <- npudist(bws=bw,tdat=d,edat=e,se=FALSE)
    expect_identical(f$dist,off$dist)
    radius <- if(type=="fixed") rep(h,3) else
      vapply(e$x,function(v)sort(abs(x-v))[[h]],numeric(1))
    oracle <- vapply(seq_len(nrow(e)),function(j) {
      L <- vapply(as.integer(z),function(donor)
        sum(lam^abs(donor-seq_len(as.integer(e$z[j]))))/
          sum(lam^abs(donor-seq_len(nlevels(z)))),numeric(1))
      ## Retain the published Gaussian integral's historical erf multiplier.
      F <- pnorm((e$x[j]-x)/radius[j]*sqrt(2)*.7071067810)
      c(mean=mean(F*L),variance=(mean(F*L^2)-mean(F*L)^2)/n)
    },numeric(2))
    expect_equal(f$dist,unname(oracle[1,]),tolerance=1e-12)
    expect_equal(f$derr^2,unname(oracle[2,]),tolerance=1e-12)
  }
})

test_that("complete categorical smoothing has the correct variance scaling", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- seq(.01,.99,length.out=67); z <- factor(rep(letters[1:3],length.out=67))
  d <- data.frame(x=x,z=z); e <- data.frame(x=c(.2,.6),z=z[1:2])
  plain <- npudens(bws=npudensbw(dat=x,bws=.2,bandwidth.compute=FALSE),
                    tdat=x,edat=e$x,se=TRUE)
  mixed <- npudens(bws=npudensbw(dat=d,bws=c(.2,2/3),bandwidth.compute=FALSE),
                    tdat=d,edat=e,se=TRUE)
  expect_equal(mixed$dens,plain$dens/3,tolerance=1e-12)
  expect_equal(mixed$derr,plain$derr/3,tolerance=1e-12)
})

test_that("compact CDF endpoints do not manufacture categorical uncertainty", {
  old <- options(np.messages=FALSE,np.tree=TRUE)
  on.exit(options(old),add=TRUE)
  d <- data.frame(x=seq(.01,.99,length.out=193),
                  z=ordered(rep(letters[1:3],length.out=193)))
  e <- data.frame(x=c(-10,10),z=ordered(c("a","a"),levels=levels(d$z)))
  for(type in c("fixed","generalized_nn")) {
    h <- if(type=="fixed") .2 else 15L
    bw <- npudistbw(dat=d,bws=c(h,1),bwtype=type,ckertype="uniform",
                     okertype="racineliyan",bandwidth.compute=FALSE)
    ## The fixed query is exactly outside support. GNN radii expand with
    ## the query and are checked only against their actual integral target.
    f <- npudist(bws=bw,tdat=d,edat=e,se=TRUE)
    radius <- if(type=="fixed") rep(h,2) else
      vapply(e$x,function(v)sort(abs(d$x-v))[[h]],numeric(1))
    target <- vapply(seq_len(2),function(j) {
      F <- pmin(1,pmax(0,((e$x[j]-d$x)/radius[j]+1)/2))
      sqrt((mean(F)/9-mean(F/3)^2)/nrow(d))
    },numeric(1))
    expect_equal(f$derr,target,tolerance=1e-12)
    if(type=="fixed") {
      expect_equal(f$dist,c(0,1/3),tolerance=1e-14)
      expect_identical(f$derr,c(0,0))
    }
  }
})

test_that("compact CDF tree and cached joint profiles retain zero donors", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d <- data.frame(x=seq(.01,.99,length.out=193),
    a=ordered(rep(c("a","b","a"),length.out=193)),
    b=ordered(rep(c("a","b","c","a","a"),length.out=193)))
  e <- d[c(3,24,67,129,187),]; e$x <- c(.017,.117,.337,.667,.947)
  cat.cdf <- function(train,query,lambda) vapply(as.integer(train),function(donor)
    sum(lambda^abs(donor-seq_len(as.integer(query))))/
      sum(lambda^abs(donor-seq_len(nlevels(train)))),numeric(1))
  for(type in c("fixed","generalized_nn")) {
    h <- if(type=="fixed") .035 else 13L
    bw <- npudistbw(dat=d,bws=c(h,1,.31),bwtype=type,
      ckertype="epanechnikov",okertype="racineliyan",bandwidth.compute=FALSE)
    radius <- if(type=="fixed")rep(h,nrow(e)) else
      vapply(e$x,function(x)sort(abs(d$x-x))[[h]],numeric(1))
    reference <- vapply(seq_len(nrow(e)),function(j) {
      z <- (e$x[j]-d$x)/radius[j]
      ## The public integral retains these historical decimal coefficients.
      F <- ifelse(z < -sqrt(5),0,ifelse(z > sqrt(5),1,
        z*(.3354101967-.02236067978*z*z)+.5))
      L <- cat.cdf(d$a,e$a[j],1)*cat.cdf(d$b,e$b[j],.31)
      c(mean(F*L),(mean(F*L^2)-mean(F*L)^2)/nrow(d))
    },numeric(2))
    for(tree in c(FALSE,TRUE)) {
      options(np.tree=tree)
      f <- npudist(bws=bw,tdat=d,edat=e,se=TRUE)
      expect_equal(f$dist,unname(reference[1,]),tolerance=1e-12)
      expect_equal(f$derr^2,unname(reference[2,]),tolerance=1e-12)
    }
  }
})
