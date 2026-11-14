# CF-172: the response evaluation node is not another deleted donor.
r77_data <- function() list(
 x=c(.505490091559477,.589119883659296,.790222305464558,.715090542472899,
     .172868311740458,.92881602099631,.738325426881202,.83093855200801,
     .663035267321393,.393898012689315,.6485223339824,.332087422483601),
 y=c(.536363800793576,.460730826295579,.716847643469465,.623389720122072,
     .426177741326129,.659515603834944,.619751672810095,.62108321779509,
     .577242290220906,.42537193287042,.720935943159341,.415262877711799))

r77_oracle <- function(x,y,k,degree=0L,grid=y,empirical=TRUE,
                       kernel="gaussian",fixed=FALSE) {
 n<-nrow(x)
 radius<-function(distances) { cap<-length(distances)-1L; sort(distances)[min(k,cap)]*max(1,k/cap) }
 K<-switch(kernel,gaussian=dnorm,
   epanechnikov=function(z) pmax(0,3/(4*sqrt(5))*(1-z^2/5)),
   uniform=function(z) .5*(abs(z)<1))
 P<-switch(kernel,gaussian=pnorm,
   epanechnikov=function(z) {u<-pmin(sqrt(5),pmax(-sqrt(5),z));.5+3/(4*sqrt(5))*(u-u^3/15)},
   uniform=function(z) pmin(1,pmax(0,.5+.5*z)))
 mean(vapply(seq_len(n),function(i){
   keep<-setdiff(seq_len(n),i)
   hx<-vapply(seq_len(ncol(x)),function(l) if(fixed) k else radius(abs(x[keep,l]-x[i,l])),0)
   z<-sweep(sweep(as.matrix(x[keep,,drop=FALSE]),2,as.numeric(x[i,]),"-"),2,hx,"/")
   w<-apply(matrix(K(z),nrow=n-1),1,prod)
   if(degree==0L) influence<-w/sum(w) else {
     stopifnot(ncol(x)==1L)
     design<-outer(z[,1],0:degree,"^")
     influence<-as.vector(design%*%solve(crossprod(design,w*design),c(1,rep(0,degree))))*w
   }
   js<-if(empirical) setdiff(seq_len(nrow(grid)),i) else seq_len(nrow(grid))
   mean(vapply(js,function(j){
     hy<-vapply(seq_len(ncol(y)),function(l) if(fixed) k else radius(abs(y[keep,l]-grid[j,l])),0)
     u<-sweep(sweep(-as.matrix(y[keep,,drop=FALSE]),2,as.numeric(grid[j,]),"+"),2,hy,"/")
     pred<-sum(influence*apply(matrix(P(u),nrow=n-1),1,prod))
     (as.numeric(all(y[i,]<=grid[j,]))-pred)^2
   },0))
 },0))
}

test_that("conditional CDF GNN empirical CV deletes only the fold observation", {
 expect_true(spawn_mpi_slaves(1L))
 on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
 old<-options(np.tree=FALSE,np.messages=FALSE,np.extendednn=FALSE)
 on.exit(options(old),add=TRUE)
 f<-r77_data();x<-data.frame(x=f$x);y<-data.frame(y=f$y)
 for(k in c(4L,8L,10L)) for(degree in 0:2) for(empirical in c(FALSE,TRUE)) {
   b<-npcdistbw(xdat=x,ydat=y,bws=c(k,k),bwtype="generalized_nn",
       bandwidth.compute=FALSE,bwscaling=FALSE,regtype="lp",degree=degree,
       bernstein.basis=TRUE)
   grid<-if(empirical)y else y[c(8,2,5),,drop=FALSE]
   got<-getFromNamespace(".npcdistbw_eval_only","npRmpi")(x,y,bws=b,
      gydat=if(empirical)NULL else grid,do.full.integral=empirical)$objective
   expect_lt(abs(got-r77_oracle(x,y,k,degree,grid,empirical)),2e-9)
 }
})

test_that("conditional CDF GNN response mask and row owners retain the deleted-fold target", {
 expect_true(spawn_mpi_slaves(1L))
 on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
 old<-options(np.tree=FALSE,np.messages=FALSE,np.extendednn=FALSE)
 on.exit(options(old),add=TRUE)
 f<-r77_data();x<-data.frame(x=f$x)
 for(d in c(1L,2L,4L,5L)) {
   y<-as.data.frame(vapply(seq_len(d),function(l) f$y[(seq_along(f$y)+2*l-3L)%%length(f$y)+1L],numeric(length(f$y))))
   b<-npcdistbw(xdat=x,ydat=y,bws=rep(8,d+1L),bwtype="generalized_nn",bandwidth.compute=FALSE,bwscaling=FALSE)
   got<-getFromNamespace(".npcdistbw_eval_only","npRmpi")(x,y,bws=b,do.full.integral=TRUE)$objective
   expect_lt(abs(got-r77_oracle(x,y,8)),2e-9)
 }
})


test_that("non-beta conditional CDF recovery uses the deleted domain on every grid", {
  caps <- getFromNamespace(".npcdistbw_nn_recovery_caps", "npRmpi")
  schedule <- getFromNamespace(".np_nn_ordinary_schedule", "npRmpi")
  for(n in c(3L,12L,48L)) for(q in 0:2) {
    setup <- list(nobs=n,cont_flat=seq_len(q+1L))
    template <- list(type="generalized_nn",iycon=rep(TRUE,q),
                     cxkertype="gaussian",cykertype="gaussian")
    for(full in c(FALSE,TRUE)) for(external in c(FALSE,TRUE)) {
      grid <- if(external)data.frame(y=.5)else NULL
      proposal <- caps(template,setup,grid,full)
      incumbent <- caps(template,setup,grid,full,incumbent=TRUE)
      expect_identical(proposal,rep.int(n-2L,q+1L))
      expect_identical(incumbent,rep.int(n-1L,q+1L))
      for(k in unique(c(max(1L,n-3L),n-2L,n-1L,n))) {
        point <- rep(as.double(k),q+1L)
        got <- schedule(point,seq_along(point),proposal,incumbent)
        if(k==n || k==n-2L) expect_length(got,0L) else {
          expect_true(length(got)>0L)
          expect_identical(tail(got,1L)[[1L]],as.double(proposal))
        }
      }
    }
  }
})


test_that("deleted-fold CDF extended scales and radius refusals are protected", {
 expect_true(spawn_mpi_slaves(1L))
 on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
 old<-options(np.tree=FALSE,np.messages=FALSE,np.extendednn=TRUE,np.macMseries.accelerate=FALSE)
 on.exit(options(old),add=TRUE)
 f<-r77_data();x<-data.frame(x=f$x);y<-data.frame(y=f$y)
 for(accelerate in c(FALSE,TRUE))for(kernel in c("gaussian","uniform","epanechnikov"))for(k in c(11L,13L))for(full in c(FALSE,TRUE)) {
  options(np.macMseries.accelerate=accelerate)
  b<-npcdistbw(xdat=x,ydat=y,bws=c(k,k),bwtype="generalized_nn",bandwidth.compute=FALSE,bwscaling=FALSE,cxkertype=kernel,cykertype=kernel)
  grid<-if(full)y else y[c(8,2,5),,drop=FALSE]
  got<-getFromNamespace(".npcdistbw_eval_only","npRmpi")(x,y,bws=b,gydat=if(full)NULL else grid,do.full.integral=full)$objective
  expect_lt(abs(got-r77_oracle(x,y,k,grid=grid,empirical=full,kernel=kernel)),2e-9)
 }
 options(np.extendednn=FALSE)
 for(counts in list(c(1L,8L),c(11L,8L),c(8L,11L))) {
  b<-npcdistbw(xdat=x,ydat=y,bws=counts,bwtype="generalized_nn",bandwidth.compute=FALSE)
  raw<-getFromNamespace(".npcdistbw_eval_only","npRmpi")(x,y,bws=b,do.full.integral=TRUE,invalid.penalty="dbmax")$objective
  expect_identical(raw,.Machine$double.xmax)
 }
})
