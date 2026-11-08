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

# Transport siblings: the leaves above and every numerical certificate are shared.
.np_udens_beta_mpi_seed_state <- function()
  .np_rng_scope_capture(uniform.only = TRUE)

.np_udens_beta_region_dispatch <- function(theta,Veta,grid,alpha,progress) {
  seed.state<-.np_udens_beta_mpi_seed_state()
  on.exit(.np_rng_scope_restore(seed.state),add=TRUE)
  .npRmpi_bootstrap_fanout_enabled(B=nrow(grid),chunk.size=8L,what='beta nonlinear guide region')
  tasks<-lapply(.np_udens_beta_blocks(nrow(grid),8L),function(ix)
    list(start=min(ix),bsz=length(ix),grid=grid[ix,,drop=FALSE]))
  progress('nonlinear guide region',0L,nrow(grid))
  worker<-function(task,theta,Veta,alpha).np_udens_beta_region_leaf(task,theta,Veta,alpha)
  values<-.npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=2L,
    what='beta nonlinear guide region',profile.where='beta:region',
    progress.label='Beta nonlinear guide region',progress.unit='targets',
    theta=theta,Veta=Veta,alpha=alpha)
  if(!identical(dim(values),c(nrow(grid),2L))||any(!is.finite(values)))
    stop('incomplete beta plot distributed region',call.=FALSE)
  progress('nonlinear guide region',nrow(grid),nrow(grid))
  list(lower=values[,1],upper=values[,2])
}

.np_udens_beta_refit_dispatch <- function(x,theta,stream,progress) {
  seed.state<-.np_udens_beta_mpi_seed_state()
  on.exit(.np_rng_scope_restore(seed.state),add=TRUE)
  B<-stream$request$B;n<-nrow(x)
  .npRmpi_bootstrap_fanout_enabled(n=n,B=B,chunk.size=32L,what='beta bootstrap guide refits')
  workers<-.npRmpi_bootstrap_worker_count();width<-max(1L,workers+1L)
  # Bound coordinator, serialized payloads and rank-local count copies together.
  # The count wave plus the common data stay below 32 MiB, independent of B.
  available<-32*1024^2-4*as.double(length(x))*8*width
  batch<-min(32L,floor(available/(4*as.double(n)*8*width)))
  if(!is.finite(batch)||batch<1L)stop('beta plot MPI refit workspace budget',call.=FALSE)
  batch<-as.integer(batch);values<-matrix(0,B,length(theta)+1L)
  blocks<-.np_udens_beta_blocks(B,batch)
  context<-new.env(parent=emptyenv())
  context$state<-.np_plot_bootstrap_progress_begin(total=B,label='Beta bootstrap guide refits',unit='draws')
  context$done<-0L;context$use.bootstrap.done<-FALSE;complete<-FALSE
  on.exit({if(complete).np_plot_progress_end(context$state) else
    if(is.null(context$state[['plot_context']]))tryCatch(.np_progress_abort(context$state),error=function(e)NULL,interrupt=function(e)NULL)},add=TRUE)
  worker<-function(task,x,theta).np_udens_beta_refit_leaf(task,x,theta)
  for(wave in .np_udens_beta_blocks(length(blocks),width)) {
    ids<-unlist(blocks[wave],use.names=FALSE)
    tasks<-lapply(blocks[wave],function(ix)list(start=min(ix),bsz=length(ix),
      counts=.np_udens_beta_random_read(stream,min(ix),max(ix))))
    context$done<-min(ids)-1L
    fit<-.npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=length(theta)+1L,
      what='beta bootstrap guide refits',profile.where='beta:refits',
      progress.context=if(is.null(context$state))NULL else context,x=x,theta=theta)
    if(!identical(dim(fit),c(length(ids),as.integer(length(theta)+1L)))||any(!is.finite(fit)))
      stop('incomplete beta plot distributed refits',call.=FALSE)
    values[ids,]<-fit;context$done<-max(ids)
    context$state<-.np_plot_progress_tick(context$state,done=max(ids),force=max(ids)==B)
    progress('bootstrap guide refits',max(ids),B)
  }
  complete<-TRUE
  list(theta=values[,seq_along(theta),drop=FALSE],iterations=values[,length(theta)+1L])
}
