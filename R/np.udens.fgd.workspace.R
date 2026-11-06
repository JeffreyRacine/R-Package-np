.np_udens_fgd_workspace <- function(n,grid,theta,full,B=9999L,workers=0L) {
  # The Gaussian adapter uses the same 24-component joint accumulator. The
  # original common forecast bounds its retained quadrature and SD projections.
  base<-.np_udens_beta_workspace(n,grid,theta,full,B,workers)
  m<-as.double(nrow(grid));n<-as.double(n);B<-as.double(B)
  stages<-base$stages
  if(full) {
    # Pilot and empirical laws are sequential, with shared innovations on disk.
    # Charge the input, correlation/factor, eigen and reconstruction work, and
    # empirical feature/product temporaries; no dense matrix survives into MPI.
    stages<-c(stages,fgd_factor=8*(4*m*m+6*m*64+4*B+20*m),
      fgd_empirical=8*(4*m*m+2*m*n+40*m+16*n))
  }
  bytes<-max(stages)+4*1024^2
  if(!is.finite(bytes)||bytes>256*1024^2)
    stop(sprintf('Gaussian plot inference workspace forecast exceeds 256 MiB (%s stage)',names(which.max(stages))),call.=FALSE)
  list(bytes=bytes,stages=stages,cap=256*1024^2)
}

.np_udens_fgd_disk <- function(n,m,q,B,sd,draws) {
  n<-as.double(n);m<-as.double(m);B<-as.double(B)
  bytes<-if(sd)4*n*B+8*m*B+8*m*n+262144 else
    if(draws)8*m*B+16*m*m+if(q==2L)16*m*n+262144 else 262144 else 0
  if(!is.finite(bytes)||bytes>512*1024^2)stop('Gaussian plot inference exceeds owned-disk budget',call.=FALSE)
  bytes
}
