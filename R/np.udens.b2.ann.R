# Donor-centered radii use the existing canonical geometry, including identity.
.np_udens_b2_ann_geometry <- function(state) {
  z<-.np_plot_ann_rank_geometry(state$nn.bws,state$x,state$package)
  if(!isTRUE(z$valid))stop("density ANN radius or spacing pilot is unavailable",call.=FALSE)
  radius<-matrix(0,state$n,state$q)
  for(j in seq_len(state$q))radius[z$order[[j]],j]<-z$geometry[[j]]$radius
  list(radius=radius,order=z$order,geometry=z$geometry)
}

.np_udens_b2_ann_influence <- function(D, geometry, j, n) {
  g<-geometry$geometry[[j]];order<-geometry$order[[j]]
  out<-matrix(0,nrow(D),n)
  # Interval starts/stops are native zero-based half-open sorted indices.
  # Group duplicate endpoints before prefix summation, never an n-by-n map.
  for(a in seq_len(nrow(D))){
    w<--D[a,order]/g$density.sum/n
    sums<-rowsum(c(w,-w),c(g$left+1L,g$right+1L),reorder=FALSE)
    delta<-numeric(n+1L);delta[as.integer(rownames(sums))]<-sums[,1L]
    out[a,order]<-cumsum(delta)[seq_len(n)]
  }
  out
}
