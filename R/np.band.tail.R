# Pure, one-pass empirical tail envelopes over explicit scalar-target families.
.np_band_tail_prepare <- function(B, map, alpha, byte.limit=128*1024^2) {
  if(!is.numeric(B)||length(B)!=1L||!is.finite(B)||B<1||B!=floor(B)||B>.Machine$integer.max)
    stop("invalid tail draw count",call.=FALSE)
  if(!is.data.frame(map)||!all(c("row","target")%in%names(map))||!nrow(map)||
     anyNA(map)||anyDuplicated(map$row)||length(unique(map$target))!=1L)
    stop("invalid tail family map",call.=FALSE)
  if(!is.numeric(alpha)||!length(alpha)||any(!is.finite(alpha))||
     any(alpha<=0|alpha>=.5)||anyDuplicated(alpha))stop("invalid tail alpha",call.=FALSE)
  m<-nrow(map);a<-length(alpha)
  if(!is.numeric(byte.limit)||length(byte.limit)!=1L||!is.finite(byte.limit)||
     96*as.double(B)+48*as.double(m)*a>byte.limit)
    stop("tail workspace exceeds declared byte budget before allocation",call.=FALSE)
  B<-as.integer(B)
  list(B=B,map=map,alpha=alpha,kP=1+floor(B*alpha/2),kB=1+floor(B*alpha/(2*m)),
    minimum=rep(Inf,B),maximum=rep(-Inf,B),seen=rep(FALSE,m),
    point=array(0,c(m,2L,a)),bonf=array(0,c(m,2L,a)),byte.limit=byte.limit)
}

.np_band_tail_update <- function(state, signed, map, replication.ids) {
  rr<-match(map$row,state$map$row)
  if(!is.matrix(signed)||!is.numeric(signed)||nrow(signed)!=nrow(map)||ncol(signed)!=state$B||
     any(!is.finite(signed))||anyNA(rr)||anyDuplicated(rr)||
     !identical(map,state$map[rr,,drop=FALSE])||!identical(replication.ids,seq_len(state$B)))
    stop("invalid tail tile, target map or replication identities",call.=FALSE)
  if(8*as.double(length(signed))+96*as.double(state$B)+48*nrow(state$map)*length(state$alpha)>state$byte.limit)
    stop("tail tile exceeds declared byte budget",call.=FALSE)
  if(any(state$seen[rr]))stop("duplicate tail target",call.=FALSE)
  for(j in seq_along(rr)) {
    x<-signed[j,];v<-sort(x,method="radix");i<-rr[j]
    state$minimum<-pmin(state$minimum,x);state$maximum<-pmax(state$maximum,x)
    for(a in seq_along(state$alpha)) {
      state$point[i,,a]<-v[c(state$kP[a],state$B+1-state$kP[a])]
      state$bonf[i,,a]<-v[c(state$kB[a],state$B+1-state$kB[a])]
    }
  }
  state$seen[rr]<-TRUE
  state
}

.np_band_tail_finish <- function(state) {
  if(!all(state$seen))stop("missing tail targets",call.=FALSE)
  lo<-sort(state$minimum,method="radix");hi<-sort(state$maximum,method="radix")
  sim<-state$point
  for(a in seq_along(state$alpha)) {
    sim[,1L,a]<-lo[state$kP[a]];sim[,2L,a]<-hi[state$B+1-state$kP[a]]
    state$bonf[,1L,a]<-pmin(state$bonf[,1L,a],sim[,1L,a])
    state$bonf[,2L,a]<-pmax(state$bonf[,2L,a],sim[,2L,a])
  }
  list(limits=list(pointwise=state$point,simultaneous=sim,bonferroni=state$bonf),
    kP=state$kP,kS=state$kP,kB=state$kB,map=state$map,
    orientation="estimate minus truth; root-inverted endpoints",
    construction="separate empirical lower and upper family extrema")
}

.np_band_tail_protect <- function(signed, absolute) {
  if(!identical(names(signed),names(absolute)))stop("incompatible tail protection bands",call.=FALSE)
  limits<-signed;expanded<-list()
  for(band in names(signed)) {
    x<-signed[[band]];y<-absolute[[band]]
    if(!identical(dim(x),dim(y))||length(dim(x))!=3L||dim(x)[2L]!=2L||
       any(!is.finite(x))||any(!is.finite(y)))stop("invalid tail protection endpoints",call.=FALSE)
    for(a in seq_len(dim(x)[3L])) {
      if(any(x[,1L,a]>x[,2L,a])||any(y[,1L,a]>y[,2L,a]))
        stop("reversed tail protection endpoints",call.=FALSE)
      limits[[band]][,1L,a]<-pmin(x[,1L,a],y[,1L,a])
      limits[[band]][,2L,a]<-pmax(x[,2L,a],y[,2L,a])
    }
    expanded[[band]]<-limits[[band]]!=x
  }
  list(limits=limits,raw.limits=signed,expanded=expanded,
    protection="interval hull with unchanged empirical absolute-root rectangle")
}
