# Plot-only compact q2 complete ratio-pilot normal variance.
# Internal nine-target reference; public admission is qualified separately.
.np_udens_compact_normal_q2_axis <- function(kernel) {
 function(grid,axes,h,rb=NULL,inner=32L){
  parts<-index<-vector('list',2L)
  for(j in 1:2){
   target<-sort(unique(grid[,j]));index[[j]]<-as.integer(match(grid[,j],target))
   pl<-.np_udens_compact_operator_plan(target,h[j],kernel)
   a<-.np_udens_compact_operator(pl,exp(axes[[j]]$lv))-1
   z<-array(0,c(length(target),ncol(a),24L))
   for(s in 1:4){z[,,s]<-a;z[,,s+4L]<-if(j==1L)a else 1;z[,,9L+4L*(s-1L)]<-if(j==1L)1 else -.25*a}
   parts[[j]]<-z
  }
  list(parts=parts,index=index)
 }
}
.np_udens_compact_normal_q2_prepare <- function(x,h,grid,theta,outer,state,kernel,held.bytes=0){
 .np_udens_compact_normal_q2_validate(x,h,grid,kernel)
 if(length(theta)!=5L||any(!is.finite(theta))||!outer%in%c(16L,32L,64L)||
    length(held.bytes)!=1L||!is.finite(held.bytes)||held.bytes<0)
   stop('invalid compact normal q2 preparation',call.=FALSE)
  n<-nrow(x)
 axes<-lapply(1:2,function(j).np_udens_compact_rule(theta[(2*j-1):(2*j)],outer,.np_udens_compact_operator_plan(sort(unique(grid[,j])),h[j],kernel)))
 nodes<-vapply(axes,`[[`,integer(1),'nodes')
 bytes<-.np_udens_compact_normal_q2_capacity(n,nrow(grid),nodes,held.bytes)
  axis<-.np_udens_compact_normal_q2_axis(kernel);cache<-axis(grid,axes,h);Td<-.np_udens_beta_basis(x);lg<-drop(Td%*%theta)-state$A;rho<-n^(-1/3);mats<-list();shifts<-numeric(2)
 for(j in 1:2){v<-exp(axes[[j]]$lv);a<-v/rho;b<-(1-v)/rho;L<-outer(a,log(x[,j]))+outer(b,log1p(-x[,j]))-lbeta(a+1,b+1);L<-sweep(L,2,.5*lg);shifts[j]<-max(L);mats[[j]]<-exp(L-shifts[j])}
 ratio<-tcrossprod(mats[[1]],mats[[2]])/n;if(any(!is.finite(ratio))||any(ratio<=0))stop('C2 ratio cache range');ratio<-log(ratio)+sum(shifts)
 list(x=x,h=h,grid=grid,theta=theta,inner=outer,outer=outer,state=state,axes=axes,cache=cache,ratio=ratio,Td=Td,Iinv=solve(state$I),axis.cache=axis,workspace=bytes,logz=state$A-axes[[1]]$A-axes[[2]]$A)
}

.np_udens_compact_normal_q2_validate <- function(x,h,grid,kernel) {
  if(!is.matrix(x)||!is.numeric(x)||ncol(x)!=2L||nrow(x)<10L||
     !is.matrix(grid)||!is.numeric(grid)||ncol(grid)!=2L||!nrow(grid)||nrow(grid)>9L||
     !is.numeric(h)||length(h)!=2L||any(!is.finite(c(x,h,grid)))||
     any(x<=0|x>=1)||any(grid<=0|grid>=1)||any(h<=0)||
     length(kernel)!=1L||!kernel%in%c('uniform','epanechnikov'))
    stop('invalid compact normal q2 request',call.=FALSE)
  invisible(NULL)
}

# Research diagonal functional variance; existing native owner unchanged.
.np_udens_compact_normal_q2_block <- function(x,h,grid,kernel,progress=function(...)NULL,fit=NULL,held.bytes=0) {
  if(!is.matrix(x)||!is.numeric(x)||ncol(x)!=2L||nrow(x)<10L||
     !is.matrix(grid)||!is.numeric(grid)||ncol(grid)!=2L||!nrow(grid)||nrow(grid)>9L||
     !is.numeric(h)||length(h)!=2L||any(!is.finite(c(x,h,grid)))||
     any(h<=0)||any(x<=0|x>=1)||any(grid<=0|grid>=1))
    stop('invalid beta normal q2 moment request',call.=FALSE)
  n<-nrow(x);m<-nrow(grid)
  # Cheap input-size refusal precedes guide or kernel preparation.
  if(8*(4*as.double(m)*n+40*as.double(m)^2)>192*1024^2)
    stop('beta normal q2 plot workspace exceeds 192 MiB',call.=FALSE)
  progress('original joint guide',0L,1L)
  if(is.null(fit))fit<-.np_udens_beta_q2_fit(x);theta<-fit$theta
  progress('original joint guide',1L,1L);levels<-vector('list',2L);workspace<-0
  raw<-function(p).np_udens_beta_joint_raw_native(p$cache,p$axes,p$theta[5L],p$logz,p$ratio,FALSE)
  blocks<-.np_udens_beta_blocks(m,16L)
  for(level in 1:2) {
    outer<-c(16L,32L)[level]
    width<-min(16L,m)
    workspace<-max(workspace,.np_udens_beta_normal_q2_budget(n,width,theta,outer))
    progress(paste('joint reference pilot',level),0L,1L)
    first<-blocks[[1L]]
    p<-.np_udens_compact_normal_q2_prepare(x,h,grid[first,,drop=FALSE],theta,outer,fit$state,kernel,held.bytes)
    initial<-raw(p)[[2L]];target<-initial$mean[length(first)+seq_len(5L)]
    .np_udens_compact_normal_q2_capacity(n,length(first),held.bytes=held.bytes+as.double(object.size(p)))
    fq<-.np_udens_beta_normal_q2_fit(target,theta)
    workspace<-max(workspace,.np_udens_beta_normal_q2_budget(n,width,fq$theta,outer))
    pq<-.np_udens_compact_normal_q2_prepare(x,h,grid[first,,drop=FALSE],fq$theta,outer,fq$state,kernel,held.bytes+as.double(object.size(p)))
    workspace<-max(workspace,p$workspace,pq$workspace)
    z<-list(center=numeric(m),se=numeric(m),theta=theta,thetaQ=fq$theta,
      meanQ=numeric(m),D=matrix(0,m,5L),variance.checks=list())
    for(rows in blocks) {
      if(!identical(rows,first)) {
        p$grid<-pq$grid<-grid[rows,,drop=FALSE]
        p$cache<-.np_udens_compact_normal_q2_axis(kernel)(p$grid,p$axes,h)
        pq$cache<-.np_udens_compact_normal_q2_axis(kernel)(pq$grid,pq$axes,h)
      }
      ii<-seq_along(rows);tt<-length(rows)+seq_len(5L)
      qr<-if(identical(rows,first))initial else raw(p)[[2L]]
      .np_udens_beta_check(qr$mean[tt],target,'q2 blocked reference score',1e-10)
      gq<-raw(pq)[[1L]]
      Tg<-.np_udens_beta_basis(grid[rows,,drop=FALSE]);gf<-exp(drop(Tg%*%fq$theta)-fq$state$A)
      G<-sweep(sweep(Tg,2L,fq$state$mean),1L,gf,'*')
      D<-G-gq$cross[ii,,drop=FALSE];J<-D%*%solve(fq$state$I)
      cross<-qr$cross[ii,,drop=FALSE];V<-qr$cross[tt,,drop=FALSE]
      components<-cbind(kernel=qr$variance[ii],cross=2*rowSums(J*cross),
        guide=rowSums((J%*%V)*J))
      variance<-(components[,1L]+components[,2L]+components[,3L])/n
      if(any(!is.finite(variance))||any(variance<0))
        stop('beta normal q2 complete variance is not nonnegative',call.=FALSE)
      base<-.np_udens_beta_q2_moment_block(p,as.integer(ii),FALSE)
      z$center[rows]<-base$center;z$se[rows]<-sqrt(variance)
      z$meanQ[rows]<-qr$mean[ii];z$D[rows,]<-D
      if(is.null(z$reference))z$reference<-list(normal=exp(qr$logmass),V=V,I=fq$state$I,
        score_error=fq$score_error,iterations=fq$iterations)
      else .np_udens_beta_check(V,z$reference$V,'q2 blocked reference covariance')
      z$variance.checks[[length(z$variance.checks)+1L]]<-list(components=components,identities=base$identities)
      progress(paste('complete joint variance',level),max(rows),m)
    }
    levels[[level]]<-z
    rm(p,pq,initial,qr,gq,base)
  }
  fields<-c('center','se','thetaQ','meanQ','D')
  levels[[1L]]$checks<-setNames(lapply(fields,function(k).np_udens_beta_check(
    levels[[1L]][[k]],levels[[2L]][[k]],paste('N3 refinement',k))),fields)
  levels[[1L]]$workspace<-workspace
  levels[[1L]]
}


# Bound live invocation workspace before allocating target tensors or ratio caches.
# nodes are the actual target-knot-split axes, never marginal beta-only proxies.
.np_udens_compact_normal_q2_capacity <- function(n,m,nodes=NULL,held.bytes=0) {
  if(!is.numeric(n)||length(n)!=1L||!is.finite(n)||n<10||n!=floor(n)||
     !is.numeric(m)||length(m)!=1L||!is.finite(m)||m<1||m>9||m!=floor(m)||
     !is.numeric(held.bytes)||length(held.bytes)!=1L||!is.finite(held.bytes)||held.bytes<0)
    stop('invalid compact normal q2 workspace geometry',call.=FALSE)
  stages<-c(guide=28*600000,input=4*as.double(m)*n+40*m*m)
  if(!is.null(nodes)) {
    if(!is.numeric(nodes)||length(nodes)!=2L||any(!is.finite(nodes))||
       any(nodes<2|nodes!=floor(nodes))||any(nodes>1600)||prod(as.double(nodes))>2560000)
      stop('compact normal q2 joint refinement node cap',call.=FALSE)
    stages<-c(stages,preparation=20*as.double(n)*max(nodes)+4*prod(as.double(nodes))+
      96*as.double(m)*sum(nodes)+40*m*m+150*m+60*n)
  }
  bytes<-held.bytes+8*max(stages)
  if(!is.finite(bytes)||bytes>192*1024^2)
    stop('compact normal q2 plot workspace exceeds 192 MiB',call.=FALSE)
  bytes
}

# Target-specific quadrature and reference fitting stay local to each block.
# Only the original x-only guide is shared; no numerical fallback or global cache.
.np_udens_compact_normal_q2_moments <- function(x,h,grid,kernel,
                                               progress=function(...)NULL) {
  if(!is.matrix(grid)||!is.numeric(grid)||ncol(grid)!=2L||!nrow(grid)||
     any(!is.finite(grid))||any(grid<=0|grid>=1))
    stop('invalid compact normal q2 evaluation grid',call.=FALSE)
  m<-nrow(grid);n<-nrow(x)
  .np_udens_compact_normal_q2_validate(x,h,grid[seq_len(min(9L,m)),,drop=FALSE],kernel)
  held<-8*(400*as.double(m)+96*as.double(n)+200000)
  workspace<-.np_udens_compact_normal_q2_capacity(n,min(9L,m),held.bytes=held)
  progress('original joint guide',0L,1L);fit<-.np_udens_beta_q2_fit(x)
  progress('original joint guide',1L,1L)
  out<-list(center=numeric(m),se=numeric(m),theta=fit$theta,meanQ=numeric(m),
    D=matrix(0,m,5L),blocks=list(),workspace=workspace)
  blocks<-.np_udens_beta_blocks(m,9L)
  for(b in seq_along(blocks)) {
    rows<-blocks[[b]]
    notify<-function(stage,done,total)progress(paste('block',b,'of',length(blocks),stage),done,total)
    z<-.np_udens_compact_normal_q2_block(x,h,grid[rows,,drop=FALSE],kernel,notify,fit,held)
    for(k in c('center','se','meanQ'))out[[k]][rows]<-z[[k]]
    out$D[rows,]<-z$D
    out$blocks[[b]]<-list(rows=rows,thetaQ=z$thetaQ,reference=z$reference,
      checks=z$checks,variance=z$variance.checks)
    out$workspace<-max(out$workspace,z$workspace)
    rm(z)
    progress('completed target blocks',max(rows),m)
  }
  out
}
