# Exact NN bootstrap: counts are held fixed; radii and scales are refitted.
.np_udens_b2_nn_refit <- function(state, index) {
  if(length(index)!=state$n||anyNA(index)||any(index<1L|index>state$n)||any(index!=floor(index)))
    stop("invalid density NN resample indices",call.=FALSE)
  out<-state;out$category.cache<-out$replicate.index<-NULL;out$x<-state$x[index,,drop=FALSE];rownames(out$x)<-NULL
  if(state$r+state$s>0L)for(j in seq_along(out$category))
    out$category[[j]]$training<-state$category[[j]]$training[index]
  if(state$type=="adaptive_nn"&&state$q>0L)out$ann<-.np_udens_b2_ann_geometry(out)
  out
}

.np_udens_b2_nn_boot_chunk <- function(state, grid, base, first, last, B, law, blocklen) {
  indices<-.np_udens_b2_indices(state$n,first,last,B,law,blocklen)
  out<-matrix(0,nrow(grid),length(indices));degenerate<-integer(nrow(grid))
  for(a in seq_along(indices)){
    fit<-.np_udens_b2_nn_refit(state,indices[[a]])
    fit$category.cache<-base$category.cache;fit$replicate.index<-indices[[a]]
    z<-.np_udens_b2_nn_contribution(fit,grid,base$use.corrected)
    delta<-z$center-base$center;variance<-.np_udens_b2_variance(z$factor,law,blocklen)
    if(any(!is.finite(variance))||any(variance<0))stop("nonfinite density NN replicate variance",call.=FALSE)
    zero<-variance==0
    if(any(zero)){
      certified<-.np_udens_b2_block_zero(z$factor,z$factor,blocklen,law)
      if(any(zero&(!certified|delta==0)))stop("undefined density NN bootstrap pivot",call.=FALSE)
      degenerate<-degenerate+as.integer(zero)
    }
    out[,a]<-delta/sqrt(variance)
    if(anyNA(out[,a]))stop("undefined density NN bootstrap pivot",call.=FALSE)
  }
  list(signed=out,degeneracy=degenerate,cancellation=integer(nrow(grid)))
}

.np_udens_b2_nn_bootstrap <- function(state, grid, law="inid",blocklen=NULL,
                                      alpha=c(.05,.1),B=9999L,tile=25L,corrected=TRUE) {
  if(!state$type%in%c("generalized_nn","adaptive_nn"))stop("unadmitted exact density NN bootstrap type",call.=FALSE)
  if(state$q==0L)return(.np_udens_b2_bootstrap(state,grid,law,blocklen,alpha,B,tile,corrected))
  state$cache.enabled<-TRUE
  owner<-new.env(parent=environment(.np_udens_b2_bootstrap))
  owner$.np_udens_bc_grid<-.np_udens_b2_nn_grid
  owner$.np_udens_bc_contribution<-function(state,grid,corrected){
    z<-.np_udens_b2_nn_contribution(state,grid,corrected);z$use.corrected<-corrected
    if(state$r+state$s>0L)z$category.cache<-lapply(state$category,function(c)
      .np_udens_bc_categorical_rows(c,grid[[c$index]]))
    z
  }
  environment(owner$.np_udens_bc_contribution)<-owner
  owner$.np_udens_b2_boot_chunk<-.np_udens_b2_nn_boot_chunk
  dispatch<-.np_udens_b2_boot_dispatch;environment(dispatch)<-owner
  owner$.np_udens_b2_boot_dispatch<-dispatch
  entry<-.np_udens_b2_bootstrap;environment(entry)<-owner
  out<-entry(state,grid,law,blocklen,alpha,B,tile,corrected)
  out$settings$construction<-"density-bc-v2-exact-NN-bootstrap"
  out$settings$type<-state$type;out$settings$counts<-state$requested
  out
}
