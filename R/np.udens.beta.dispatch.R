# Deterministic leaves contain no MPI, sampling, estimator call or retained cache.
.np_udens_beta_region_leaf <- function(task,theta,Veta,alpha) {
  r<-if(ncol(task$grid)==1L).np_udens_beta_region_q1(theta,Veta,task$grid[,1],alpha) else
    .np_udens_beta_region_q2(theta,Veta,task$grid,alpha)
  result<-cbind(r$lower,r$upper)
  if(!identical(dim(result),c(task$bsz,2L))||any(!is.finite(result))||any(result[,1]>result[,2]))
    stop('invalid beta plot region task result',call.=FALSE)
  result
}

.np_udens_beta_refit_leaf <- function(task,x,theta) {
  if(ncol(x)==1L)result<-.np_udens_beta_q1_refits(x[,1],task$counts) else {
    fit<-.np_udens_beta_q2_refits_native(x,task$counts,theta)
    result<-cbind(fit$theta,fit$iterations)
  }
  if(!identical(dim(result),c(task$bsz,as.integer(length(theta)+1L)))||any(!is.finite(result)))
    stop('invalid beta plot refit task result',call.=FALSE)
  result
}

.np_udens_beta_region_dispatch <- function(theta,Veta,grid,alpha,progress) {
  values<-matrix(0,nrow(grid),2L)
  for(ix in .np_udens_beta_blocks(nrow(grid),8L)) {
    progress('nonlinear guide region',min(ix)-1L,nrow(grid))
    task<-list(start=min(ix),bsz=length(ix),grid=grid[ix,,drop=FALSE])
    values[ix,]<-.np_udens_beta_region_leaf(task,theta,Veta,alpha)
    progress('nonlinear guide region',max(ix),nrow(grid))
  }
  list(lower=values[,1],upper=values[,2])
}

.np_udens_beta_refit_dispatch <- function(x,theta,stream,progress) {
  B<-stream$request$B;values<-matrix(0,B,length(theta)+1L)
  batch<-min(32L,max(1L,floor(8*1024^2/(nrow(x)*8))))
  for(ix in .np_udens_beta_blocks(B,batch)) {
    task<-list(start=min(ix),bsz=length(ix),counts=.np_udens_beta_random_read(stream,min(ix),max(ix)))
    values[ix,]<-.np_udens_beta_refit_leaf(task,x,theta)
    progress('bootstrap guide refits',max(ix),B)
  }
  list(theta=values[,seq_along(theta),drop=FALSE],iterations=values[,length(theta)+1L])
}
