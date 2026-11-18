# O(B) standardized family extrema; draws stay with their existing owner.
.np_band_local_prepare <- function(B, map, alpha, byte.limit=128*1024^2) {
  if(!is.numeric(B)||length(B)!=1L||!is.finite(B)||B<2||B!=floor(B)||B>.Machine$integer.max)
    stop("invalid local-envelope draw count",call.=FALSE)
  if(!is.data.frame(map)||!all(c("row","target")%in%names(map))||!nrow(map)||
     anyNA(map)||anyDuplicated(map$row)||length(unique(map$target))!=1L)
    stop("invalid local-envelope family map",call.=FALSE)
  if(!is.numeric(alpha)||!length(alpha)||any(!is.finite(alpha))||
     any(alpha<=0|alpha>=.5)||anyDuplicated(alpha))
    stop("invalid local-envelope alpha",call.=FALSE)
  if(!is.numeric(byte.limit)||length(byte.limit)!=1L||!is.finite(byte.limit)||
     96*as.double(B)+32*as.double(nrow(map))*length(alpha)>byte.limit)
    stop("local-envelope workspace exceeds declared byte budget before allocation",call.=FALSE)
  list(B=as.integer(B),map=map,alpha=alpha,minimum=rep(Inf,B),
    maximum=rep(-Inf,B),seen=rep(FALSE,nrow(map)),byte.limit=byte.limit)
}

.np_band_local_update <- function(state, signed, scale, map, replication.ids) {
  rr<-match(map$row,state$map$row)
  if(!is.matrix(signed)||!is.numeric(signed)||nrow(signed)!=nrow(map)||ncol(signed)!=state$B||
     any(!is.finite(signed))||anyNA(rr)||anyDuplicated(rr)||
     !identical(map,state$map[rr,,drop=FALSE])||!identical(replication.ids,seq_len(state$B))||
     !is.numeric(scale)||length(scale)!=length(rr)||any(!is.finite(scale))||any(scale<0))
    stop("invalid local-envelope tile, scales, target map or replication identities",call.=FALSE)
  if(8*as.double(length(signed))+96*as.double(state$B)+32*nrow(state$map)*length(state$alpha)>state$byte.limit)
    stop("local-envelope tile exceeds declared byte budget",call.=FALSE)
  if(any(state$seen[rr]))stop("duplicate local-envelope target",call.=FALSE)
  for(j in seq_along(rr)) {
    x<-signed[j,]
    if(scale[j]==0) {
      if(any(x!=0))stop("zero target scale with nonzero local-envelope draws",call.=FALSE)
    } else x<-x/scale[j]
    if(any(!is.finite(x)))stop("nonfinite local-envelope standardized draws",call.=FALSE)
    state$minimum<-pmin(state$minimum,x);state$maximum<-pmax(state$maximum,x)
  }
  state$seen[rr]<-TRUE
  state
}

.np_band_local_finish <- function(state) {
  if(!all(state$seen))stop("missing local-envelope targets",call.=FALSE)
  lo<-sort(state$minimum,method="radix");hi<-sort(state$maximum,method="radix")
  absolute<-sort(pmax(-state$minimum,state$maximum),method="radix")
  k<-1+floor(state$B*state$alpha/2);ka<-ceiling(state$B*(1-state$alpha))
  list(cutoffs=rbind(pmin(lo[k],-absolute[ka]),pmax(hi[state$B+1-k],absolute[ka])),
    map=state$map,alpha=state$alpha,kS=k,kA=ka,
    units="observed-SE standardized signed roots",
    construction="local signed family extrema protected by their absolute maximum")
}

.np_band_local_protect <- function(limits, local, map, scale) {
  if(!identical(names(limits),c("pointwise","simultaneous","bonferroni"))||
     !is.numeric(scale)||length(scale)!=nrow(map)||any(!is.finite(scale))||any(scale<0))
    stop("invalid local-envelope protection input",call.=FALSE)
  seen<-rep(FALSE,nrow(map))
  for(family in local) {
    rr<-which(map$target==family$map$target[1L])
    if(!identical(map[rr,,drop=FALSE],family$map)||any(seen[rr])||
       !identical(dim(family$cutoffs),c(2L,length(family$alpha)))||
       any(!is.finite(family$cutoffs))||any(family$cutoffs[1L,]>family$cutoffs[2L,]))
      stop("invalid local-envelope protection family",call.=FALSE)
    for(band in c("simultaneous","bonferroni")) {
      x<-limits[[band]]
      if(!identical(dim(x),c(nrow(map),2L,length(family$alpha)))||any(!is.finite(x)))
        stop("invalid local-envelope protection limits",call.=FALSE)
      for(a in seq_along(family$alpha)) {
        low<-scale[rr]*family$cutoffs[1L,a];high<-scale[rr]*family$cutoffs[2L,a]
        if(any(!is.finite(c(low,high))))stop("nonfinite local-envelope endpoint offsets",call.=FALSE)
        limits[[band]][rr,1L,a]<-pmin(x[rr,1L,a],low)
        limits[[band]][rr,2L,a]<-pmax(x[rr,2L,a],high)
      }
    }
    seen[rr]<-TRUE
  }
  if(!all(seen))stop("missing local-envelope protection family",call.=FALSE)
  limits
}
