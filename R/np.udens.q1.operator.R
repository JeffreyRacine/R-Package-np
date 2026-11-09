# Resolved plot-only q1 operator moments. Guide fitting and statistical laws
# remain unchanged; both mean/projection and covariance use the same measure.
.np_udens_q1_operator_plan <- function(theta,h,grid) {
  plan<-.np_udens_q1_variance_plan(theta,h,grid)
  plan$refined.nodes<-96L*(length(plan$edges)-1L)
  if(plan$refined.nodes>8192L)stop('beta plot marginal quadrature node budget',call.=FALSE)
  plan
}

.np_udens_q1_operator_rule <- function(theta,plan,outer,keep.log=FALSE) {
  if(length(outer)!=1L || !outer%in%c(128L,256L))
    stop('invalid q1 plot operator refinement',call.=FALSE)
  order<-if(outer==128L)64L else 96L
  plan$order<-as.integer(order);plan$nodes<-order*(length(plan$edges)-1L)
  edges<-plan$edges;order<-plan$order
  if(!order%in%c(64L,96L) || any(!is.finite(edges)) || any(diff(edges)<=0) ||
     length(edges)<2L || plan$nodes!=order*(length(edges)-1L) || plan$nodes>8192L)
    stop('invalid q1 plot operator plan',call.=FALSE)
  gl<-.np_udens_beta_rule(c(1,1),order)
  t<-unlist(lapply(seq_len(length(edges)-1L),function(i)edges[i]+diff(edges[i+0:1])*gl[,1]))
  qw<-unlist(lapply(seq_len(length(edges)-1L),function(i)diff(edges[i+0:1])*gl[,2]))
  softplus<-function(t)pmax(t,0)+log1p(exp(-abs(t)))
  lv<--softplus(-t);l1v<--softplus(t);T<-cbind(lv,l1v)
  lw<-drop(T%*%theta)+lv+l1v+log(qw);A<-.np_udens_beta_logsum(lw)
  p<-exp(lw-A);total<-sum(p);p<-p/total;mean<-drop(crossprod(p,T));D<-sweep(T,2,mean)
  out<-list(A=A,mean=mean,I=crossprod(D,D*p),lv=lv,l1v=l1v,T=T,prob=p,nodes=length(t))
  if(keep.log)out$log.prob<-lw-A-log(total)
  out
}
