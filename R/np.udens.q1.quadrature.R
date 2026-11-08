# Plot-only q1 second-moment integration. Original first moments, guide fits,
# corrected centers and the independent outer256 certificate are unchanged.
.np_udens_q1_variance_plan <- function(theta,h,grid) {
  if(length(theta)!=2L || any(!is.finite(theta)) || any(theta<=-1) ||
     length(h)!=1L || !is.finite(h) || h<=0 || !length(grid) ||
     any(!is.finite(grid)) || any(grid<=0|grid>=1))
    stop('invalid q1 plot variance geometry',call.=FALSE)
  a<-theta[1]+1;b<-theta[2]+1;budget<-70+abs(lbeta(a,b))
  TL<-max(40,budget/a);TR<-max(40,budget/b)
  if(max(TL,TR)>1e7)stop('unresolved beta plot guide boundary',call.=FALSE)
  fixed<-c(0,.5,1,2,3,5,10,20,40,80,160,320,640,1280,2560,5120,10240,
    20480,40960,81920,163840,327680,655360,1310720,2621440,5242880)
  # Equal angular spacing resolves narrow interior beta peaks; unit logit
  # panels resolve the distinct endpoint regime. Plan once for the whole family.
  K<-ceiling((pi/2)/(8*h))
  if(K>128L)stop('beta plot marginal quadrature node budget',call.=FALSE)
  angular<-(seq_len(K-1L)/K)*(pi/2)
  features<-2*(log(sin(angular))-log(cos(angular)))
  tt<-qlogis(sort(unique(as.double(grid))))
  endpoint<-sort(unique(c(seq.int(-8L,8L),unlist(lapply(unique(floor(tt)),function(t)t+(-2:3))))))
  extra<-c(features,endpoint);extra<-extra[extra> -TL & extra<TR]
  edges<-sort(unique(c(-TL,-fixed[fixed<TL],fixed[fixed<TR],TR,extra)))
  nodes<-64L*(length(edges)-1L)
  if(nodes>8192L)stop('beta plot marginal quadrature node budget',call.=FALSE)
  list(edges=edges,order=64L,nodes=nodes)
}

.np_udens_q1_variance_rule <- function(theta,plan) {
  edges<-plan$edges;order<-plan$order
  if(!identical(order,64L) || any(!is.finite(edges)) || any(diff(edges)<=0) ||
     length(edges)<2L || plan$nodes!=order*(length(edges)-1L) || plan$nodes>8192L)
    stop('invalid q1 plot variance plan',call.=FALSE)
  gl<-.np_udens_beta_rule(c(1,1),order)
  t<-unlist(lapply(seq_len(length(edges)-1L),function(i)edges[i]+diff(edges[i+0:1])*gl[,1]))
  qw<-unlist(lapply(seq_len(length(edges)-1L),function(i)diff(edges[i+0:1])*gl[,2]))
  softplus<-function(t)pmax(t,0)+log1p(exp(-abs(t)))
  lv<--softplus(-t);l1v<--softplus(t);T<-cbind(lv,l1v)
  lw<-drop(T%*%theta)+lv+l1v+log(qw);A<-.np_udens_beta_logsum(lw)
  p<-exp(lw-A);p<-p/sum(p);mean<-drop(crossprod(p,T));D<-sweep(T,2,mean)
  list(A=A,mean=mean,I=crossprod(D,D*p),lv=lv,l1v=l1v,T=T,prob=p,nodes=length(t))
}

.np_udens_q1_score_columns <- function(grid,coordinates,h,rb,inner,score) {
  sizes<-vapply(coordinates,function(z)length(z$lv),integer(1));ends<-cumsum(sizes)
  total<-sum(as.double(sizes))
  if(any(sizes<1L) || any(vapply(coordinates,function(z)length(z$lv)!=length(z$l1v),logical(1))))
    stop('invalid q1 score columns',call.=FALSE)
  if(!is.finite(total) || 8*total>32*1024^2)stop('beta plot score tile budget',call.=FALSE)
  rows<-max(1L,floor(32*1024^2/(8*total)))
  lv<-unlist(lapply(coordinates,`[[`,'lv'),use.names=FALSE)
  l1v<-unlist(lapply(coordinates,`[[`,'l1v'),use.names=FALSE)
  out<-lapply(sizes,function(n)matrix(0,length(grid),n))
  # Native preparation depends on target/h, not on donor-column ordering.
  # Batching reuses that preparation without changing any column arithmetic.
  for(first in seq.int(1L,length(grid),by=rows)) {
    ix<-seq.int(first,min(length(grid),first+rows-1L))
    values<-score(grid[ix],lv,l1v,h*h,rb,inner)
    for(j in seq_along(sizes))out[[j]][ix,]<-values[,ends[j]-sizes[j]+seq_len(sizes[j]),drop=FALSE]
  }
  out
}

.np_udens_q1_variance_workspace <- function(n,rows,old.nodes,new.nodes) {
  total<-as.double(old.nodes)+new.nodes+n
  tile<-min(rows,max(1,floor(32*1024^2/(8*total))))
  native.tile<-min(tile,8,max(1,floor(32*1024^2/(24*8*total))))
  # Retained output columns, current score tile and 24-component native scratch.
  score<-8*(rows*total+tile*total+24*native.tile*total+8*total)+8*1024^2
  moments<-8*(8*rows*new.nodes+4*rows*n+5*rows^2+32*new.nodes+12*n)
  max(score,moments)
}

.np_udens_q1_ratio_weights <- function(x,theta,oldA,st) {
 n<-length(x);rho<-n^(-2/5);lgx<-drop(cbind(log(x),log1p(-x))%*%theta)-oldA;lp<-numeric(st$nodes)
 for(start in seq.int(1L,st$nodes,by=32L)){
  ix<-seq.int(start,min(st$nodes,start+31L));v<-exp(st$lv[ix]);a<-v/rho;b<-(1-v)/rho
  logs<-outer(a,log(x))+outer(b,log1p(-x))-lbeta(a+1,b+1);logs<-sweep(logs,2,lgx)
  lp[ix]<-apply(logs,1,.np_udens_beta_logsum)-log(n)
 }
 lw<-log(st$prob)+lp;normal<-.np_udens_beta_logsum(lw);pp<-exp(lw-normal);pp<-pp/sum(pp)
 list(prob=pp,normal=exp(normal))
}
