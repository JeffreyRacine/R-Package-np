.np_udens_beta_joint_block_native <- function(cache,i,j) {
  .Call(C_np_plot_udens_beta_joint_block,cache$parts,lapply(cache$index,as.integer),as.integer(i),as.integer(j))
}

.np_udens_beta_joint_raw_native <- function(cache,axes,gamma,logz,ratio,full) {
  logs<-lapply(axes,function(a)cbind(a$lv,a$l1v,log(a$prob)))
  raw<-.Call(C_np_plot_udens_beta_joint_moments,cache$parts,lapply(cache$index,as.integer),
    logs,as.double(c(gamma,logz)),ratio,as.logical(full))
  lapply(raw,function(z)setNames(z,c('mean','variance','cross','cov','logmass')))
}
