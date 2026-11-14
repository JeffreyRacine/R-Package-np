test_that("conditional empty-window disclosure preserves the requested family", {
  withr::local_options(np.messages=FALSE)
  set.seed(9781);d<-data.frame(x=runif(150));d$y<-d$x+rexp(150)
  b<-npcdensbw(xdat=d["x"],ydat=d["y"],bws=c(.25,.3),
               bandwidth.compute=FALSE,cxkertype="uniform",cykertype="uniform")
  f<-npcdens(b,txdat=d["x"],tydat=d["y"]);saved<-f
  total.empty<-0L
  for(surface in c(TRUE,FALSE)) {
    p<-suppressWarnings(plot(f,errors="asymptotic",band="all",neval=12L,
                              perspective=surface,output="data"))
    for(j in seq_along(p$panels)) {
      g<-p$panels[[j]]$grid
      nx<-vapply(g[[1]],function(z)sum(abs(d$x-z)<b$xbw),0L)
      nj<-mapply(function(x,y)sum(abs(d$x-x)<b$xbw & abs(d$y-y)<b$ybw),g[[1]],g[[2]])
      empty<-nx>0 & nj==0
      total.empty<-total.empty+sum(empty)
      expect_true(any(!empty))
      point<-as.data.frame(p,panel=j,band="pointwise")
      sim<-as.data.frame(p,panel=j,band="simultaneous")
      expect_equal(point$estimate[nx>0],(nj/(2*b$ybw*nx))[nx>0],tolerance=1e-12)
      expect_true(all(is.na(point$lower[empty])))
      expect_true(all(point$interval_status[empty]=="undefined"))
      expect_true(all(is.finite(point$lower[!empty])))
      expect_true(all(is.na(sim$lower)))
      expect_equal(nrow(point),nrow(g))
    }
  }
  expect_gt(total.empty,0L)
  expect_identical(f,saved)
})

test_that("empty-window admission excludes distinct inference contracts", {
  withr::local_options(np.messages=FALSE)
  d<-data.frame(x=seq(0,1,length.out=30),y=seq(0,1,length.out=30))
  b<-npcdensbw(xdat=d["x"],ydat=d["y"],bws=c(.2,.2),
               bandwidth.compute=FALSE,cxkertype="uniform",cykertype="uniform")
  admit<-getFromNamespace(".np_plot_conditional_empty_mask",'np')
  ey<-data.frame(y=c(-.1,0,.5,1,1.1));zero<-rep(0,5)
  call<-function(bw=b,center="estimate",gr=FALSE,proper=FALSE,method="asymptotic")
    suppressWarnings(admit(bw,zero,ey,center,gr,proper,method))
  expect_identical(call(),rep(TRUE,5))
  for(mode in c("bias-corrected","none"))
    expect_false(any(call(center=mode)))
  expect_false(any(call(gr=TRUE)));expect_false(any(call(proper=TRUE)))
  expect_false(any(call(method="none")))
  bound<-b;bound$cykerbound<-"fixed";bound$cykerlb<-0;bound$cykerub<-1
  expect_identical(call(bound),c(FALSE,TRUE,TRUE,TRUE,FALSE))
  nn<-b;nn$type<-"generalized_nn"
  expect_false(any(call(nn,method="bootstrap")))
  signed<-b;signed$cykertype<-"epanechnikov";signed$cykerorder<-4L
  expect_false(any(call(signed)))
})
