test_that("native vector fits expose caller variable labels", {
 if(environmentName(environment(npudist))=="npRmpi"){
  skip_on_cran();if(!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES","1"))))stop("MPI unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
 }
 old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
 set.seed(831);n<-90L;x1<-runif(n);x2<-rnorm(n);yy<-x2+sin(x1)+rnorm(n,sd=.2)
 build<-function(frames=FALSE) {
  if(frames){x1<-data.frame(x1=x1);x2<-data.frame(x2=x2)}
  list(pl=npplreg(tydat=yy,txdat=x2,tzdat=x1,bws=matrix(.3,2,1),bandwidth.compute=FALSE),
    cd=npcdens(txdat=x1,tydat=yy,bws=c(.3,.3),bandwidth.compute=FALSE),
    cc=npcdist(txdat=x1,tydat=yy,bws=c(.3,.3),bandwidth.compute=FALSE),
    sc=npscoef(tydat=yy,txdat=x2,tzdat=x1,bws=.3,bandwidth.compute=FALSE,betas=TRUE),
    qr=npqreg(txdat=x1,tydat=yy,bws=c(.3,.3),bandwidth.compute=FALSE,tau=.5),
    ls=nplsqreg(tydat=yy,txdat=x1,bws=.3,bandwidth.compute=FALSE,tau=.5))
 }
 a<-build();b<-build(TRUE)
 labels<-function(z)list(pl=c(z$pl$data.xnames,z$pl$data.znames,names(coef(z$pl))),
  cd=c(names(z$cd$xeval),names(z$cd$yeval)),cc=c(names(z$cc$xeval),names(z$cc$yeval)),
  sc=c(names(z$sc$eval$exdat),names(z$sc$eval$ezdat)),qr=names(z$qr$xeval),
  ls=c(z$ls$reg.bws$xnames,z$ls$reg.bws$ynames,z$ls$fit$xnames,z$ls$fit$ynames))
 wanted<-list(pl=c("x2","x1","x2"),cd=c("x1","yy"),cc=c("x1","yy"),sc=c("x2","x1"),qr="x1",ls=c("x1","yy","x1","yy"))
 la<-labels(a);lb<-labels(b)
 for(nm in names(a)){
  expect_identical(la[[nm]],wanted[[nm]],info=nm)
  expect_identical(la[[nm]],lb[[nm]],info=paste("frame",nm))
  expect_equal(unname(fitted(a[[nm]])),unname(fitted(b[[nm]])),tolerance=1e-12,info=nm)
 }
 expect_true(is.finite(coef(a$pl)["x2"]))
 if(identical(Sys.getenv("N5_FULL"),"true")){
   for(nm in names(a)){
     restored<-unserialize(serialize(a[[nm]],NULL));expect_equal(unname(fitted(restored)),unname(fitted(a[[nm]])),tolerance=0,info=nm)
     out<-plot(a[[nm]],errors="none",neval=5L,output="data")
     expect_false(any(grepl("npRmpi_autod",unlist(lapply(out$panels,function(p)p$variables$name)))))
     expect_false(any(grepl("npRmpi_autod",capture.output(summary(a[[nm]])))))
   }
 }
 if(nzchar(Sys.getenv("N5_SNAPSHOT")))saveRDS(a,Sys.getenv("N5_SNAPSHOT"))
})

test_that("explicit token-shaped names and leased bandwidth identity are preserved", {
 if(environmentName(environment(npudist))!="npRmpi") return(invisible(NULL))
 skip_on_cran();if(!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES","1"))))stop("MPI unavailable")
 on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
 options(np.messages=FALSE);set.seed(832);x<-runif(60);y<-sin(x)+rnorm(60,sd=.1)
 d<-data.frame(x=x);names(d)<-".__npRmpi_autod_xdat_2"
 for(f in list(npregbw,npindexbw)) {
  b<-suppressWarnings(f(xdat=d,ydat=y,bws=if(identical(f,npindexbw)) c(.3,1) else .3,bandwidth.compute=FALSE))
  expect_identical(b$xnames,names(d))
 }
 b<-npcdensbw(xdat=d,ydat=data.frame(y=y),bws=c(.3,.3),bandwidth.compute=FALSE)
 fit<-npcdens(bws=b,txdat=d,tydat=data.frame(y=y))
 expect_identical(names(fit$xeval),names(d));expect_identical(fit$bws$xnames,names(d))
 r<-npscoefbw(xdat=x,zdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE)
 frozen<-serialize(r,NULL);a<-npscoef(bws=r,txdat=x,tzdat=x,tydat=y)
 c<-npscoef(bws=r,txdat=x,tzdat=x,tydat=y)
 expect_identical(serialize(r,NULL),frozen);expect_identical(fitted(a),fitted(c))
 bad<-r;bad$bw<-bad$bw*1.01;bad$bandwidth$z<-bad$bandwidth$z*1.01
 # Existing contract: a mutated object is not reused through its stale lease;
 # the ordinary materialization path remains available for the edited value.
 resolve<-getFromNamespace(".npRmpi_lease_resolve","npRmpi")
 expect_false(is.null(resolve(r,"npscoef","bws")))
 expect_null(resolve(bad,"npscoef","bws"))
 got<-npscoef(bws=bad,txdat=x,tzdat=x,tydat=y)
 ref<-npscoef(bws=getFromNamespace(".npRmpi_autodispatch_untag","npRmpi")(bad),txdat=x,tzdat=x,tydat=y)
 expect_identical(fitted(got),fitted(ref))
 expect_false(identical(fitted(got),fitted(a)))
 expect_identical(serialize(r,NULL),frozen)
})


test_that("descriptive frame finalization preserves exact row-name representation", {
  helper <- getFromNamespace(".npRmpi_native_identity_apply", "npRmpi")
  for (rows in list(NULL, c("row-a", "row-b", "row-c"))) {
    frame <- data.frame(value = 1:3, row.names = rows)
    names(frame) <- ".__npRmpi_autod_txdat_2"
    original.rows <- attr(frame, "row.names", exact = TRUE)
    got <- helper(list(xeval = frame), c(".__npRmpi_autod_txdat_2" = "value"))
    expect_identical(attr(got$xeval, "row.names", exact = TRUE), original.rows)
    expect_identical(got$xeval[[1L]], frame[[1L]])
    expect_identical(names(got$xeval), "value")
  }
})
