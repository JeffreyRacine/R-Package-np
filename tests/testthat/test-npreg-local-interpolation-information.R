test_that("local square supported polynomial designs report unidentified uncertainty", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  for(p in 2:3) for(kernel in c("uniform","epanechnikov"))
   for(bernstein in c(FALSE,TRUE)) {
    x<-data.frame(x=rep(c(0,4,8),each=p)+rep(if(p==2)c(0,1)else c(0,.3,1),3))
    y<-sin(seq_len(nrow(x)))
    b<-npregbw(xdat=x,ydat=y,bws=1.1,regtype="lp",degree=p-1,
      bernstein.basis=bernstein,ckertype=kernel,bandwidth.compute=FALSE)
    K<-npksum(txdat=x,bws=1.1,ckertype=kernel,return.kernel.weights=TRUE)$kw
    expect_identical(as.integer(colSums(K!=0)),rep(p,nrow(x)))
    W<-outer(x$x,0:(p-1),`^`)
    expect_true(all(vapply(seq_len(nrow(x)),function(i)
      qr(W[K[,i]!=0,,drop=FALSE])$rank==p,logical(1L))))
    H<-npreghat(b,txdat=x)
    expect_true(all(attr(H,"ridge.used")==0))
    expect_warning(fit<-npreg(b,txdat=x,tydat=y,se=TRUE,gradients=TRUE),
                        "leaves no residual information")
    expect_true(all(is.na(se(fit))))
    expect_true(all(is.na(gradients(fit,se=TRUE))))
    expect_equal(fitted(fit),y,tolerance=1e-10)
    perm<-rev(seq_len(nrow(x)))
    expect_warning(other<-npreg(b,txdat=x[perm,,drop=FALSE],tydat=y[perm],
        exdat=x,se=TRUE,gradients=TRUE),"leaves no residual information")
    expect_identical(is.na(se(other)),is.na(se(fit)))
    expect_equal(fitted(other),fitted(fit),tolerance=1e-10)
    expect_equal(gradients(other),gradients(fit),tolerance=1e-10)
  }
})

test_that("extra donors and positive ridge do not acquire the square certificate", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  for(shape in c("extra","duplicate")) {
    localx<-if(shape=="extra")c(0,.2,.6,1)else c(0,0,1)
    x<-data.frame(x=rep(c(0,4,8),each=length(localx))+rep(localx,3))
    y<-sin(seq_len(nrow(x)))
    b<-npregbw(xdat=x,ydat=y,bws=1.1,regtype="lp",degree=2,
      ckertype="uniform",bandwidth.compute=FALSE)
    H<-npreghat(b,txdat=x)
    if(shape=="duplicate")expect_true(all(attr(H,"ridge.used")>0))else
      expect_true(all(attr(H,"ridge.used")==0))
    fit<-npreg(b,txdat=x,tydat=y,se=TRUE,gradients=TRUE)
    expect_true(all(is.finite(se(fit))))
    expect_true(all(is.finite(gradients(fit,se=TRUE))))
    expect_equal(drop(H%*%y),fitted(fit),tolerance=1e-10)
  }
})

test_that("compact NN uncertainty preserves observation permutation identities", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(294);n<-27L
  x<-data.frame(x=rnorm(n),z=runif(n,-1,1));y<-sin(x$x)+x$z+rnorm(n)
  e<-x[c(4,8,15,21),];ix<-sample(n);ie<-c(4,2,1,3)
  # Retain the documented uncertainty and inherited uniform-order advisories.
  uniform.advisories<-0L
  fitfun<-function(...)withCallingHandlers(npreg(...),warning=function(w){
    if(grepl("leaves no residual information",conditionMessage(w),fixed=TRUE))
      invokeRestart("muffleWarning")
    if(grepl("ignoring kernel order specified with uniform kernel type",
             conditionMessage(w),fixed=TRUE)){
      uniform.advisories<<-uniform.advisories+1L
      invokeRestart("muffleWarning")
    }
  })
  for(bwt in c("fixed","generalized_nn","adaptive_nn"))for(reg in c("lc","ll","lp")){
    b<-npregbw(xdat=x,ydat=y,bws=if(bwt=="fixed")c(.8,.6)else c(13,15),
      bandwidth.compute=FALSE,bwtype=bwt,regtype=reg,
      degree=if(reg=="lp")c(2,1)else NULL,ckertype="uniform")
    a<-fitfun(b,txdat=x,tydat=y,exdat=e,gradients=TRUE,se=TRUE)
    t<-fitfun(b,txdat=x[ix,],tydat=y[ix],exdat=e,gradients=TRUE,se=TRUE)
    v<-fitfun(b,txdat=x,tydat=y,exdat=e[ie,],gradients=TRUE,se=TRUE)
    expect_equal(fitted(a),fitted(t),tolerance=1e-10)
    expect_equal(se(a),se(t),tolerance=1e-10)
    expect_equal(gradients(a,se=TRUE),gradients(t,se=TRUE),tolerance=1e-10)
    expect_equal(se(a),se(v)[order(ie)],tolerance=1e-10)
  }
  expect_identical(uniform.advisories,0L)
})
