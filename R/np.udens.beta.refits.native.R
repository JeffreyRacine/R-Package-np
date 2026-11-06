.np_udens_beta_q2_refits_native <- function(x,counts,theta) {
  q<-.np_udens_beta_series_rule(theta,32L)
  target<-crossprod(counts/nrow(x),.np_udens_beta_basis(x))
  initial<-cbind(.np_udens_beta_q1_refits(x[,1],counts)[,1:2,drop=FALSE],.np_udens_beta_q1_refits(x[,2],counts)[,1:2,drop=FALSE],0)
  axes<-lapply(q$axes,function(a)cbind(a$lv,a$l1v,log(a$prob)))
  result<-.Call(C_np_plot_udens_beta_joint_refits,axes,as.double(q$base),as.double(q$A),target,initial)
  list(theta=result[,1:5,drop=FALSE],iterations=result[,6],rule=q)
}

.np_udens_beta_series_projection <- function(q,h,n,grid,inner) {
  cache<-.np_udens_beta_axis_cache(grid,q$axes,h,n^(-13/54),inner)
  list(cache=cache,rule=q,grid=grid)
}

.np_udens_beta_series_means_native <- function(prep,thetas) {
  q<-prep$rule;cache<-prep$cache
  axes<-lapply(q$axes,function(a)cbind(a$lv,a$l1v,log(a$prob)))
  ans<-.Call(C_np_plot_udens_beta_series_means,cache$parts,cache$index,axes,
    as.double(q$base),as.double(q$A),matrix(as.double(thetas),nrow(thetas)))
  names(ans)<-c('mu','A','terms','tail_bound')
  ans$fitted<-exp(sweep(thetas%*%t(.np_udens_beta_basis(prep$grid)),1,ans$A))
  ans
}
