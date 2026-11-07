.np_udens_beta_covariance_native <- function(prep,B,law,progress=NULL) {
  axes<-lapply(prep$axes,function(a)cbind(a$lv,a$l1v,log(a$prob)))
  .Call(C_np_plot_udens_beta_joint_covariance,prep$cache$parts,prep$cache$index,axes,
    c(prep$theta[5],prep$logz,nrow(prep$x)),prep$ratio,B%*%prep$Iinv,
    as.integer(law),progress)
}

.np_udens_beta_dense_compare <- function(main,refined,name) {
  if(!identical(dim(main),dim(refined)))stop('beta plot covariance dimensions',call.=FALSE)
  error<-scale<-0
  for(ix in .np_udens_beta_blocks(nrow(main))) {
    error<-max(error,max(abs(main[ix,,drop=FALSE]-refined[ix,,drop=FALSE])))
    scale<-max(scale,max(abs(refined[ix,,drop=FALSE])))
  }
  limit<-1e-12+1e-9*scale
  if(!is.finite(error)||error>limit)stop(sprintf('beta plot %s refinement failed: %.17g > %.17g',name,error,limit),call.=FALSE)
  c(error=error,limit=limit)
}
