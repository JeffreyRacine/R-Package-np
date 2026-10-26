# Exact circular selected-block law; no residual pilot or resample relabeling.
.np_udens_b2_block_length <- function(blocklen,n) {
  if(!is.numeric(blocklen)||length(blocklen)!=1L||!is.finite(blocklen)||
     blocklen<1||blocklen>n||blocklen!=floor(blocklen))
    stop("density block bootstrap requires an explicit integer block length in [1,n]",call.=FALSE)
  as.integer(blocklen)
}

.np_udens_b2_block_zero <- function(raw, F, blocklen, law) {
  constant <- apply(raw,1L,function(z)all(z==z[1L]))
  if(law!="fixed"||blocklen==1L)return(constant)
  n<-ncol(F);Z<-F-rowMeans(F);certificate<-rep(TRUE,nrow(F))
  for(width in unique(c(blocklen,n%%blocklen))){
    if(!width)next
    sums<-matrix(0,nrow(F),n)
    for(k in seq_len(width))sums<-sums+Z[,((seq_len(n)+k-2L)%%n)+1L,drop=FALSE]
    certificate<-certificate&apply(sums,1L,function(z)all(z==0))
  }
  constant|certificate
}

.np_udens_b2_block_chunk <- function(state,grid,base,first,last,B,law,blocklen) {
  indices<-.np_udens_b2_indices(state$n,first,last,B,law,blocklen)
  out<-matrix(0,nrow(grid),length(indices));degenerate<-integer(nrow(grid))
  for(j in seq_along(indices)){
    idx<-indices[[j]];raw<-base$corrected[,idx,drop=FALSE]
    delta<-rowSums(base$factor[,idx,drop=FALSE])
    center<-raw[,1L]+rowMeans(raw-raw[,1L]);F<-(raw-center)/state$n
    variance<-.np_udens_b2_variance(F,law,blocklen)
    if(any(!is.finite(variance))||any(variance<0)||any(!is.finite(delta)))
      stop("nonfinite or negative density block replicate variance",call.=FALSE)
    zero<-variance==0
    if(any(zero)){
      cert<-.np_udens_b2_block_zero(raw,F,blocklen,law)
      if(any(zero&!cert))stop("uncertified zero density block variance",call.=FALSE)
      if(any(zero&delta==0))stop("undefined zero-over-zero density bootstrap pivot",call.=FALSE)
      degenerate<-degenerate+as.integer(zero)
    }
    out[,j]<-delta/sqrt(variance)
    if(any(!is.finite(out[!zero,j])))stop("uncertified infinite density block pivot",call.=FALSE)
  }
  list(signed=out,degeneracy=degenerate,cancellation=integer(nrow(grid)))
}
