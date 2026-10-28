# Pure finite-array calibration. Family identity is supplied, never inferred.
.np_band_rank_prepare <- function(B, map, alpha, byte.limit=128*1024^2) {
  if(length(B)!=1L || !is.numeric(B) || !is.finite(B) || B<1 || B!=floor(B) || B>.Machine$integer.max)
    stop("invalid rank draw count",call.=FALSE)
  if(!is.data.frame(map) || !all(c("row","target") %in% names(map)) ||
     !nrow(map) || anyNA(map) || anyDuplicated(map$row) || length(unique(map$target))!=1L)
    stop("invalid rank family map",call.=FALSE)
  if(!is.numeric(alpha)||!length(alpha)||any(!is.finite(alpha))||any(alpha<=0|alpha>=.5)||anyDuplicated(alpha))
    stop("invalid rank alpha",call.=FALSE)
  # Account for values, sort copies, indices/tie runs and K even with one target.
  if(!is.numeric(byte.limit)||length(byte.limit)!=1L||!is.finite(byte.limit)||byte.limit<1||80*as.double(B)>byte.limit)
    stop("rank workspace exceeds declared byte budget before allocation",call.=FALSE)
  m<-nrow(map);B<-as.integer(B);a<-length(alpha)
  list(B=B,map=map,alpha=alpha,kP=1+floor(B*alpha/2),kB=1+floor(B*alpha/(2*m)),
    K=rep.int(B,B),seen=rep(FALSE,m),replayed=rep(FALSE,m),
    point=array(0,c(m,2L,a)),bonf=array(0,c(m,2L,a)),sim=array(0,c(m,2L,a)),
    receipt=matrix(0,m,3L),ties=integer(m),byte.limit=byte.limit)
}

.np_band_rank_update <- function(state, signed, map, replication.ids, replay=FALSE) {
  rr<-match(map$row,state$map$row)
  if(!is.matrix(signed)||!is.numeric(signed)||nrow(signed)!=nrow(map)||ncol(signed)!=state$B||
     any(!is.finite(signed))||anyNA(rr)||anyDuplicated(rr)||
     !identical(map,state$map[rr,,drop=FALSE])||!identical(replication.ids,seq_len(state$B)))
    stop("invalid rank tile, target map or replication identities",call.=FALSE)
  if(8*as.double(length(signed))+72*state$B>state$byte.limit)
    stop("rank tile exceeds declared byte budget",call.=FALSE)
  seen<-if(replay)state$replayed else state$seen
  if(any(seen[rr]) || (replay && is.null(state$kS)))stop("duplicate or premature rank tile",call.=FALSE)
  for(j in seq_along(rr)) {
    x<-signed[j,];o<-order(x,method="radix");v<-x[o];i<-rr[j]
    receipt<-c(sum(x),sum(x*x),sum(x*replication.ids))
    if(replay) {
      if(!identical(receipt,state$receipt[i,]))stop("signed replay receipt mismatch",call.=FALSE)
      for(a in seq_along(state$alpha))state$sim[i,,a]<-v[c(state$kS[a],state$B+1-state$kS[a])]
    } else {
      # Each sorted tie run contributes inclusive <= and >= counts in O(B).
      runs<-rle(v)$lengths;ends<-cumsum(runs);starts<-ends-runs+1L
      rank<-integer(state$B);rank[o]<-rep(pmin(ends,state$B+1L-starts),runs)
      state$K<-pmin(state$K,rank);state$ties[i]<-state$B-length(runs)
      state$receipt[i,]<-receipt
      for(a in seq_along(state$alpha)) {
        state$point[i,,a]<-v[c(state$kP[a],state$B+1-state$kP[a])]
        state$bonf[i,,a]<-v[c(state$kB[a],state$B+1-state$kB[a])]
      }
    }
  }
  if(replay)state$replayed[rr]<-TRUE else state$seen[rr]<-TRUE
  state
}

.np_band_rank_cutoff <- function(state) {
  if(!all(state$seen))stop("missing rank targets",call.=FALSE)
  ordered<-sort(state$K,method="radix")
  r<-ceiling(state$B*(1-state$alpha))
  state$kS<-pmin(state$kP,ordered[state$B-r+1L])
  if(any(state$kS<state$kB))stop("rank union-bound invariant failed",call.=FALSE)
  state$content<-vapply(state$kS,function(k)sum(state$K>=k),integer(1))
  state$K<-NULL
  state
}

.np_band_rank_finish <- function(state) {
  if(!all(state$replayed))stop("missing rank replay targets",call.=FALSE)
  list(limits=list(pointwise=state$point,simultaneous=state$sim,bonferroni=state$bonf),
    kP=state$kP,kS=state$kS,kB=state$kB,content=state$content,
    ties=state$ties,map=state$map,orientation="estimate minus truth; root-inverted endpoints")
}

.np_band_rank_invert <- function(center, scale, limits, deterministic=rep(FALSE,length(scale))) {
  if(!is.numeric(center)||!is.numeric(scale)||length(center)!=length(scale)||
     any(!is.finite(center))||any(!is.finite(scale))||any(scale<0))
    stop("invalid rank center or scale",call.=FALSE)
  if(!is.logical(deterministic)||length(deterministic)!=length(scale)||
     anyNA(deterministic)||any(deterministic & scale!=0))
    stop("invalid deterministic target certificate",call.=FALSE)
  lapply(limits,function(x) {
    if(length(dim(x))!=3L||dim(x)[1L]!=length(center)||dim(x)[2L]!=2L||any(!is.finite(x)))
      stop("invalid root endpoint dimensions",call.=FALSE)
    if(any(scale==0 & !deterministic) && any(x[scale==0 & !deterministic,,,drop=FALSE]!=0))
      stop("zero target scale with nonzero signed roots",call.=FALSE)
    out<-x
    for(a in seq_len(dim(x)[3L])) {
      if(any(x[,1L,a]>x[,2L,a]))stop("reversed root endpoints",call.=FALSE)
      out[,1L,a]<-center-scale*x[,2L,a]
      out[,2L,a]<-center-scale*x[,1L,a]
    }
    out
  })
}
