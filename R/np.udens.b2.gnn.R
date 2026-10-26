# Private random-radius density prototype. Existing fixed and native owners are untouched.
.np_udens_b2_nn_prepare <- function(bws, x, max.bytes=256*1024^2) {
  type<-bws[["type",exact=TRUE]]
  if(!type%in%c("generalized_nn"))stop("unadmitted density NN type",call.=FALSE)
  numeric<-vapply(x,is.numeric,logical(1L));raw<-.np_physical_bandwidth(bws)
  package<-getNamespaceName(environment(.np_physical_bandwidth))
  rank<-if(any(numeric)).Call("C_np_plot_ann_rank_decode",as.double(raw[numeric]),as.integer(nrow(x)),PACKAGE=package)else NULL
  if(!is.null(rank)&&(any(rank$extended)||any(rank$distance.scale!=1)))
    stop("density NN inference requires interior unscaled counts",call.=FALSE)
  # Reuse the fixed adapter solely for retained categorical/type metadata.
  spec<-bws;spec[["type"]]<-"fixed"
  spec[["bw"]][numeric]<-1;spec[["bandwidth"]][numeric]<-1
  state<-.np_udens_bc_prepare(spec,x,max.bytes)
  state$numeric<-numeric;state$continuous<-NULL;state$type<-type
  state$nn.bws<-bws;state$lookup<-if(is.null(rank))integer()else rank$lookup
  state$requested<-raw[numeric];state$package<-package
  state$construction<-"density-bc-v2-random-radius"
  state$metadata$type<-type;state$metadata$radius<-"recomputed from the sample at fixed counts"
  state$metadata$rank.remainder.order<-if(length(state$lookup))1/state$lookup else numeric()
  state
}

.np_udens_b2_nn_grid <- function(state, grid) {
  if(!is.data.frame(grid)||!nrow(grid)||!identical(names(grid),names(state$x))||anyNA(grid)||
     !all(vapply(grid[state$numeric],is.numeric,logical(1L)))||
     any(!is.finite(as.matrix(grid[state$numeric]))))stop("invalid density NN evaluation frame",call.=FALSE)
  if(state$r+state$s>0L).np_udens_bc_mixed_grid(state,grid)
  invisible(TRUE)
}

.np_udens_b2_gnn_geometry <- function(state, grid) {
  n<-state$n;q<-state$q;m<-nrow(grid);h<-density<-matrix(0,m,q)
  for(j in seq_len(q)){
    k<-state$lookup[j];w<-min(floor(n^(2/3)),floor(min(k,n-k)/2))
    if(w<1L||k-w<1L||k+w>n)stop("density GNN spacing pilot is non-interior",call.=FALSE)
    x<-state$x[state$numeric][[j]];t<-grid[state$numeric][[j]]
    for(a in seq_len(m)){
      d<-abs(x-t[a]);order<-sort.int(d,partial=c(k-w,k,k+w))
      h[a,j]<-order[k];density[a,j]<-(2*w/n)/(order[k+w]-order[k-w])
    }
  }
  if(any(!is.finite(h))||any(h<=0)||any(!is.finite(density))||any(density<=0))
    stop("density GNN radius or spacing pilot is unavailable",call.=FALSE)
  list(radius=h,density.sum=density)
}

.np_udens_b2_nn_rows <- function(state, grid, geometry, corrected=TRUE) {
  m<-nrow(grid);n<-state$n;q<-state$q
  K<-Q<-matrix(1,m,n);radius<-matrix(0,m,n)
  x<-state$x[state$numeric];g<-grid[state$numeric]
  for(j in seq_len(q)){
    h<-geometry$radius[,j];if(state$type=="adaptive_nn")h<-matrix(rep(h,each=m),m,n)
    u<-outer(g[[j]],x[[j]],"-")/h
    K<-K*dnorm(u)/h;radius<-radius+u^2
  }
  Q<-K
  if(state$r+state$s>0L)for(z in state$category){
    v<-.np_udens_bc_categorical_rows(z,grid[[z$index]])
    K<-K*v$original;Q<-Q*if(corrected)v$corrected else v$original
  }
  P<-if(corrected)1+q/2-radius/2 else matrix(1,m,n)
  L<-Q*P
  if(any(!is.finite(L))||any(!is.finite(K)))stop("nonfinite density NN contributions",call.=FALSE)
  list(original=K,corrected=L,Q=Q,P=P)
}

.np_udens_b2_nn_contribution <- function(state, grid, corrected=TRUE) {
  .np_udens_b2_nn_grid(state,grid)
  if(!is.logical(corrected)||length(corrected)!=1L||is.na(corrected))stop("corrected must be TRUE or FALSE",call.=FALSE)
  .np_udens_bc_budget(state,nrow(grid),1L,nrow(grid),FALSE)
  if(state$q==0L)return(.np_udens_bc_contribution(state,grid,corrected))
  geometry<-.np_udens_b2_gnn_geometry(state,grid)
  z<-.np_udens_b2_nn_rows(state,grid,geometry,corrected)
  phi<-z$corrected;x<-state$x[state$numeric];g<-grid[state$numeric]
  for(j in seq_len(state$q)){
    h<-geometry$radius[,j];u<-outer(g[[j]],x[[j]],"-")/h
    D<-z$Q/h*((u^2-1)*z$P+if(corrected)u^2 else 0)
    inside<-abs(outer(g[[j]],x[[j]],"-"))<=h
    phi<-phi-rowMeans(D)/geometry$density.sum[,j]*inside
  }
  center<-z$corrected[,1]+rowMeans(z$corrected-z$corrected[,1])
  F<-(phi-phi[,1]-rowMeans(phi-phi[,1]))/state$n
  if(any(!is.finite(F)))stop("nonfinite density NN influence",call.=FALSE)
  list(original=z$original,corrected=z$corrected,original.center=rowMeans(z$original),
       center=center,factor=F,se=sqrt(rowSums(F^2)),geometry=geometry)
}

.np_udens_b2_nn_gaussian <- function(state,grid,alpha=c(.05,.1),B=9999L,tile=25L,corrected=TRUE) {
  owner<-new.env(parent=environment(.np_udens_bc_gaussian))
  owner$.np_udens_bc_grid<-.np_udens_b2_nn_grid
  owner$.np_udens_bc_contribution<-.np_udens_b2_nn_contribution
  for(nm in c(".np_udens_bc_covariance",".np_udens_bc_gaussian")){
    f<-get(nm,owner,inherits=TRUE);environment(f)<-owner;assign(nm,f,owner)
  }
  entry<-.np_udens_bc_inference;environment(entry)<-owner
  out<-entry(state,grid,law="gaussian",alpha=alpha,B=B,tile=tile,corrected=corrected)
  out$settings$h<-out$settings$auxiliary.h<-NULL
  out$settings$counts<-state$requested;out$settings$covariance<-"complete direct and random-radius observation influence"
  out
}
