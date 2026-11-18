# Invocation-owned beta q2 geometry. Shared guide rules and other kernels retain
# their original quadrature. Each base measure receives its own positive rule.
.np_udens_beta_q2_axes <- function(theta,h,grid,order) {
  if(!is.numeric(theta)||length(theta)!=5L||any(!is.finite(theta))||any(theta[1:4]<=-1)||
     !is.numeric(h)||length(h)!=2L||any(!is.finite(h))||any(h<=0)||
     !is.matrix(grid)||!is.numeric(grid)||ncol(grid)!=2L||!nrow(grid)||
     any(!is.finite(grid))||any(grid<=0|grid>=1)||
     length(order)!=1L||!is.numeric(order)||!is.finite(order)||!order%in%c(32L,36L))
    stop("invalid beta q2 quadrature geometry",call.=FALSE)
  targets<-lapply(1:2,function(j)sort(unique(as.double(grid[,j]))))
  plans<-lapply(1:2,function(j).np_udens_q1_variance_plan(theta[(2*j-1):(2*j)],h[j],targets[[j]]))
  nodes<-vapply(plans,function(p)as.integer(order*(length(p$edges)-1L)),integer(1L))
  if(any(nodes>1600L)||prod(as.double(nodes))>2560000)
    stop("beta plot joint refinement node cap",call.=FALSE)
  gl<-.np_udens_beta_rule(c(1,1),as.integer(order))
  axes<-lapply(1:2,function(j) {
    edges<-plans[[j]]$edges;delta<-diff(edges)
    t<-as.vector(outer(gl[,1],delta)+rep(head(edges,-1L),each=order))
    weights<-as.vector(outer(gl[,2],delta))
    softplus<-function(t)pmax(t,0)+log1p(exp(-abs(t)))
    lv<--softplus(-t);l1v<--softplus(t);T<-cbind(lv,l1v)
    lw<-drop(T%*%theta[(2*j-1):(2*j)])+lv+l1v+log(weights)
    A<-.np_udens_beta_logsum(lw);prob<-exp(lw-A);prob<-prob/sum(prob)
    mu<-drop(crossprod(prob,T));dm<-sweep(T,2,mu)
    list(A=A,mean=mu,I=crossprod(dm,dm*prob),lv=lv,l1v=l1v,T=T,
      prob=prob,nodes=nodes[j],log.weights=log(weights))
  })
  out<-list(theta=theta,h=h,targets=targets,order=as.integer(order),axes=axes)
  .np_udens_beta_q2_axes_check(out,theta,h,grid,order)
  out
}

.np_udens_beta_q2_axes_check <- function(prepared,theta,h,grid,order) {
  if(!is.list(prepared)||!identical(prepared$theta,theta)||!identical(prepared$h,h)||
     !identical(prepared$order,as.integer(order))||length(prepared$axes)!=2L||
     length(prepared$targets)!=2L||!is.matrix(grid)||ncol(grid)!=2L||!nrow(grid)||
     any(!is.finite(grid))||any(grid<=0|grid>=1))
    stop("incompatible beta q2 prepared axes",call.=FALSE)
  nodes<-numeric(2L)
  for(j in 1:2) {
    z<-prepared$axes[[j]]
    if(!is.list(z)||length(z$nodes)!=1L||!is.numeric(z$nodes))
      stop("invalid beta q2 prepared axis",call.=FALSE)
    nodes[j]<-z$nodes
    if(length(nodes[j])!=1L||!is.finite(nodes[j])||nodes[j]<1||nodes[j]!=floor(nodes[j])||
       nodes[j]>1600||!all(grid[,j]%in%prepared$targets[[j]])||
       !is.matrix(z$T)||!identical(dim(z$T),c(as.integer(nodes[j]),2L))||
       length(z$prob)!=nodes[j]||length(z$lv)!=nodes[j]||length(z$l1v)!=nodes[j]||
       length(z$log.weights)!=nodes[j]||any(!is.finite(c(z$T,z$prob,z$log.weights,z$A)))||
       any(z$prob<0)||any(diff(z$lv-z$l1v)<0)||
       any(diff(z$lv-z$l1v)==0 & (diff(z$lv)!=0|diff(z$l1v)!=0))||
       any(z$lv>0)||any(z$l1v>0)||
       !identical(unname(z$T),unname(cbind(z$lv,z$l1v))))
      stop("invalid beta q2 prepared axis",call.=FALSE)
    # Distinct positive quadrature nodes may round to identical paired logs.
    # Retain their separate weights and arithmetic; reject inversions or
    # ambiguous equal logits rather than requiring artificial separation.
    th<-theta[(2*j-1):(2*j)]
    lw<-drop(z$T%*%th)+z$lv+z$l1v+z$log.weights
    A<-.np_udens_beta_logsum(lw);pp<-exp(lw-A);pp<-pp/sum(pp)
    .np_udens_beta_check(z$prob,pp,"q2 prepared positive weights",1e-12)
    ex<-.np_udens_bounded_normal_analytic(th)
    .np_udens_beta_check(z$A,A,"q2 prepared axis mass")
    .np_udens_beta_check(z$A,ex$A,"q2 prepared marginal mass")
    .np_udens_beta_check(drop(crossprod(z$prob,z$T)),z$mean,"q2 prepared axis mean",1e-10)
    .np_udens_beta_check(z$mean,ex$mean,"q2 prepared marginal mean",1e-10)
    dm<-sweep(z$T,2,z$mean)
    .np_udens_beta_check(crossprod(dm,dm*z$prob),z$I,"q2 prepared axis information")
    .np_udens_beta_check(z$I,ex$I,"q2 prepared marginal information")
  }
  if(prod(nodes)>2560000)stop("beta plot joint refinement node cap",call.=FALSE)
  prepared$axes
}

# Full covariance uses one exact preparation at a time. The invocation store
# owns these files and removes them on both success and error.
.np_udens_beta_q2_preparation_write <- function(prep,store,name) {
  bytes<-ceiling(2*as.double(object.size(prep))+65536)
  path<-.np_udens_beta_store_file(store,name,bytes)
  saveRDS(prep,path,compress=FALSE)
  size<-file.info(path)$size
  if(!is.finite(size)||size<=0||size>bytes)
    stop('beta plot preparation store size',call.=FALSE)
  list(path=path,size=size,theta=prep$theta,h=prep$h,grid=prep$grid,
       inner=prep$inner,outer=prep$outer)
}

.np_udens_beta_q2_preparation_read <- function(record) {
  size<-file.info(record$path)$size
  if(!is.finite(size)||!identical(size,record$size))
    stop('beta plot preparation store size',call.=FALSE)
  prep<-readRDS(record$path)
  if(!all(vapply(c('theta','h','grid','inner','outer'),function(k)
    identical(prep[[k]],record[[k]]),logical(1L))))
    stop('beta plot preparation store identity',call.=FALSE)
  prep
}
