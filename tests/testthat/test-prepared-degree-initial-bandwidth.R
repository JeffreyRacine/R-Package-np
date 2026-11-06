test_that("invalid decoded degree-search start is refused before worker work", {
 skip_on_cran();if(!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES","1"))))stop("MPI unavailable")
 on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
 options(np.messages=FALSE);remote<-getFromNamespace("mpi.remote.exec","npRmpi")
 set.seed(8802);d<-data.frame(x=rep(1,70),y=rnorm(70))
 dd<-data.frame(x=seq(0,1,length.out=40),y=sin(seq(0,1,length.out=40)))
 b<-npregbw(y~x,data=dd,bws=.2,bandwidth.compute=FALSE);cold<-fitted(npreg(bws=b));pids<-remote(Sys.getpid())
 e<-tryCatch(suppressWarnings(npregbw(y~x,data=d,regtype="lp",degree.select="coordinate",degree.min=0L,degree.max=2L,nmulti=1L)),error=identity)
 expect_true(inherits(e,"error"));expect_match(conditionMessage(e),"initial bandwidth returned invalid stored bandwidths")
 expect_identical(fitted(npreg(bws=b)),cold);expect_identical(remote(Sys.getpid()),pids)
 if(identical(Sys.getenv("N5_FULL"),"true")) {
  set.seed(8810);z<-data.frame(x=runif(60));z$y<-sin(z$x)+rnorm(60,sd=.1)
  a<-npregbw(y~x,data=z,regtype="lp",degree.select="coordinate",degree.min=0L,degree.max=1L,nmulti=1L)
  expect_true(all(is.finite(a$bw)));expect_true(is.finite(a$fval))
  expect_identical(fitted(npreg(bws=b)),cold)
  if(nzchar(Sys.getenv("N5_SNAPSHOT")))saveRDS(list(bw=a$bw,degree=a$degree,fval=a$fval,rng=.Random.seed),Sys.getenv("N5_SNAPSHOT"))
 }
})
