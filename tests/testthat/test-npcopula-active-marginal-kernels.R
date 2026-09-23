test_that("mixed beta copulas project only active marginal kernel metadata", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  ns <- asNamespace("npRmpi")
  args <- get(".npcopula_marginal_bw_args",ns)
  mbw <- get(".npcopula_marginal_bw",ns)
  d <- data.frame(o=ordered(rep(1:3,8)),x=seq(.1,.9,length.out=24))
  for(target in c("density","distribution")) for(ok in c("liracine","racineliyan")) {
    constructor <- if(target=="density") npudensbw else npudistbw
    b <- do.call(constructor,list(dat=d,bws=c(.2,.15),ckertype="beta",
          okertype=ok,ckerbound="fixed",ckerlb=0,ckerub=1,bandwidth.compute=FALSE))
    snapshot <- serialize(b,NULL)
    ca <- args(b,d,1L,target)
    expect_null(ca$ckertype); expect_null(ca$ckerorder)
    expect_null(ca$ckerbound); expect_null(ca$ckerlb); expect_null(ca$ckerub)
    expect_identical(ca$okertype,ok)
    cc <- args(b,d,2L,target)
    expect_identical(cc$ckertype,"beta")
    expect_identical(cc$ckerlb,b$ckerlb[2L])
    expect_identical(cc$ckerub,b$ckerub[2L])
    for(k in c(FALSE,TRUE)) {
      child <- mbw(b,d,1L,target,kbandwidth=k)
      expect_identical(child$ncon,0L)
      actual <- get(".np_ksum_unconditional_eval_exact",ns)(
        xdat=d[1L],exdat=d[1L],bws=if(k)child else get("kbandwidth",ns)(child),
        operator=if(target=="density")"normal" else "integral")
      cb <- do.call(constructor,list(dat=d[1L],bws=.2,okertype=ok,
                                     bandwidth.compute=FALSE))
      ref <- fitted(if(target=="density")npudens(cb) else npudist(cb))
      expect_equal(as.numeric(actual),as.numeric(ref),tolerance=1e-13)
    }
    a <- npcopula(b,data=d,se=TRUE)
    expect_identical(serialize(b,NULL),snapshot)
    expect_identical(a$bws$ckertype,"beta")
    expect_true(all(is.finite(se(a))))
    joint <- fitted(if(target=="density")npudens(b) else npudist(b))
    fx <- npudistbw(dat=d[2L],bws=.15,ckertype="beta",ckerbound="fixed",
                   ckerlb=0,ckerub=1,bandwidth.compute=FALSE)
    fo <- npudistbw(dat=d[1L],bws=.2,okertype=ok,bandwidth.compute=FALSE)
    expect_equal(a$u1,as.numeric(fitted(npudist(fo))),tolerance=1e-13)
    expect_equal(a$u2,as.numeric(fitted(npudist(fx))),tolerance=1e-13)
    if(target=="density") {
      dx <- npudensbw(dat=d[2L],bws=.15,ckertype="beta",ckerbound="fixed",
                     ckerlb=0,ckerub=1,bandwidth.compute=FALSE)
      do <- npudensbw(dat=d[1L],bws=.2,okertype=ok,bandwidth.compute=FALSE)
      joint <- joint/fitted(npudens(do))/fitted(npudens(dx))
    }
    expect_equal(as.numeric(fitted(a)),as.numeric(joint),tolerance=1e-13)
    expect_identical(predict(a),fitted(a))
    reverse <- do.call(constructor,list(dat=d[2:1],bws=c(.15,.2),ckertype="beta",
      okertype=ok,ckerbound="fixed",ckerlb=0,ckerub=1,bandwidth.compute=FALSE))
    z <- npcopula(reverse,data=d[2:1],se=TRUE)
    expect_equal(as.numeric(fitted(z)),as.numeric(fitted(a)),tolerance=1e-13)
    expect_equal(as.numeric(se(z)),as.numeric(se(a)),tolerance=1e-13)
  }
})

test_that("mixed beta copula grids support uncertainty and bootstrap replay", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  d <- data.frame(o=ordered(rep(1:3,8)),x=seq(.1,.9,length.out=24))
  u <- data.frame(o=c(.3,.7),x=c(.35,.65))
  for(target in c("density","distribution")) {
    constructor <- if(target=="density")npudensbw else npudistbw
    b <- do.call(constructor,list(dat=d,bws=c(.2,.15),ckertype="beta",
      ckerbound="range",bandwidth.compute=FALSE))
    fit <- npcopula(b,data=d,u=u,n.quasi.inv=20L,se=TRUE)
    expect_true(all(is.finite(fitted(fit))))
    expect_true(all(is.finite(se(fit))))
    expect_identical(predict(fit,u=u,n.quasi.inv=20L),fitted(fit))
    set.seed(937)
    out <- plot(fit,output="data",errors="bootstrap",B=3L,band="pmzsd")
    expect_equal(nrow(out),4L)
    expect_true(all(is.finite(as.matrix(out[c("center","lower","upper")]))))
    expect_true(all(out$lower<=out$upper))
  }
})
