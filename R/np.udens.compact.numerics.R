# Fixed-support, order-two compact kernels: plot-only local-cubic correction.
# The fitted estimator and its physical bandwidth are never changed.
.np_udens_compact_operator_plan <- function(target,h,kernel) {
  if(!is.numeric(target)||!length(target)||any(!is.finite(target))||any(target<0|target>1)||
     length(h)!=1L||!is.finite(h)||h<=0||length(kernel)!=1L||!kernel%in%c('uniform','epanechnikov'))
    stop('invalid compact plot operator geometry',call.=FALSE)
  span<-if(kernel=='uniform')1 else sqrt(5);s<-min(h,1)
  poly<-if(kernel=='uniform')c(.5,0,0) else c(3/(4*sqrt(5)),0,-3/(20*sqrt(5))*(s/h)^2)
  coefficients<-matrix(0,length(target),4L);mass<-condition<-numeric(length(target))
  for(i in seq_along(target)) {
    lo<-max(-target[i]/s,-span*h/s);hi<-min((1-target[i])/s,span*h/s)
    moment<-vapply(0:6,function(j)sum(poly*(hi^(j+1:3)-lo^(j+1:3))/(j+1:3)),numeric(1L))
    mass[i]<-moment[1L];M<-outer(0:3,0:3,function(j,k)moment[j+k+1L]/mass[i]);condition[i]<-rcond(M)
    if(!is.finite(condition[i])||condition[i]<1e-12)stop('compact local-polynomial moment conditioning',call.=FALSE)
    coefficients[i,]<-solve(M,c(1,0,0,0))
    .np_udens_beta_check(drop(M%*%coefficients[i,]),c(1,0,0,0),'compact cubic moments',1e-10)
  }
  list(target=target,h=h,kernel=kernel,span=span,s=s,poly=poly,coefficients=coefficients,mass=mass,condition=condition)
}

.np_udens_compact_operator <- function(plan,u,physical=NULL) {
  if(!is.null(physical) && (!identical(plan$kernel,'uniform') ||
     ncol(physical$x)!=1L || nrow(physical$x)!=length(u) ||
     nrow(physical$grid)!=length(plan$target) || length(physical$h)!=1L))
    stop('compact physical operator shape mismatch',call.=FALSE)
  out<-matrix(0,length(plan$target),length(u))
  for(i in seq_along(plan$target)) {
    z<-(u-plan$target[i])/plan$s;c<-plan$coefficients[i,]
    base<-plan$poly[1L]+plan$poly[3L]*z^2
    outside<-if(is.null(physical))abs(u-plan$target[i])>=plan$span*plan$h else
      abs(physical$x[,1L]-physical$grid[i,1L])>=physical$h
    base[outside]<-0
    out[i,]<-base*(((c[4L]*z+c[3L])*z+c[2L])*z+c[1L])/(plan$s*plan$mass[i])
  }
  out
}

.np_udens_compact_rule <- function(theta,order,plan,keep.log=FALSE) {
  a<-theta[1L]+1;b<-theta[2L]+1;budget<-70+abs(lbeta(a,b))
  TL<-max(40,budget/a);TR<-max(40,budget/b)
  if(!is.finite(TL+TR)||max(TL,TR)>1e7)stop('compact plot guide tail domain',call.=FALSE)
  fixed<-c(0,.5,1,2,3,5,10,20,40,80,160,320,640,1280,2560,5120,10240,20480,40960,81920,163840,327680,655360,1310720,2621440,5242880)
  knots<-c(plan$target-plan$span*plan$h,plan$target+plan$span*plan$h)
  cuts<-qlogis(knots[knots>0&knots<1])
  edges<-sort(unique(c(-TL,-fixed[fixed<TL],fixed[fixed<TR],TR,cuts)))
  nodes<-as.double(length(edges)-1L)*order
  if(!is.finite(nodes)||8*(40+as.integer(isTRUE(keep.log)))*nodes>64*1024^2)stop('compact plot quadrature workspace budget',call.=FALSE)
  gl<-.np_udens_beta_rule(c(1,1),order);u<-gl[,1L];w<-gl[,2L]
  z<-as.vector(outer(u,diff(edges),'*')+rep(head(edges,-1L),each=order))
  weight<-as.vector(outer(w,diff(edges),'*'))
  lv<--pmax(-z,0)-log1p(exp(-abs(z)));l1v<--pmax(z,0)-log1p(exp(-abs(z)))
  T<-cbind(lv,l1v);lw<-log(weight)+a*lv+b*l1v;A<-.np_udens_beta_logsum(lw)
  prob<-exp(lw-A);prob<-prob/sum(prob);mean<-colSums(T*prob);S<-sweep(T,2,mean)
  I<-crossprod(S,S*prob)
  tail<-sum(vapply(list(c(a,TL),c(b,TR)),function(v)
    exp(-A-v[1]*v[2])*((v[2]+1)^2/v[1]+2*(v[2]+1)/v[1]^2+2/v[1]^3),numeric(1L)))
  if(!is.finite(tail)||tail>1e-12)stop('compact plot guide tail certificate',call.=FALSE)
  result<-list(A=A,mean=mean,I=I,lv=lv,l1v=l1v,T=T,prob=prob,tail=tail,nodes=length(z))
  # Preserve the discrete normalizer before probabilities can underflow.
  # Existing callers retain their exact probability arithmetic and fields.
  if(isTRUE(keep.log))result$log.prob<-lw-A-log(sum(exp(lw-A)))
  result
}

.np_udens_compact_moments <- function(x,plan,theta,order,full,progress=function(...)NULL,physical=NULL) {
  n<-length(x);m<-length(plan$target);st<-.np_udens_compact_rule(theta,order,plan)
  # All targets share one knot-union measure; preparation is not repeated per pair.
  tile<-max(1L,min(512L,floor(4*1024^2/(8*max(n,m)))))
  blocks<-.np_udens_beta_blocks(st$nodes,tile)
  S<-sweep(st$T,2,st$mean);Ugrid<-solve(st$I,t(S));Td<-cbind(log(x),log1p(-x))
  lgx<-drop(Td%*%theta)-st$A;rho<-n^(-2/5)
  mu<-numeric(m);B<-matrix(0,m,2);lp<-numeric(st$nodes)
  for(ids in blocks) {
    P<-.np_udens_compact_operator(plan,exp(st$lv[ids]))-1
    mu<-mu+drop(P%*%st$prob[ids]);B<-B+P%*%(S[ids,,drop=FALSE]*st$prob[ids])
    v<-exp(st$lv[ids]);a<-v/rho;b<-(1-v)/rho
    logs<-outer(a,log(x))+outer(b,log1p(-x))-lbeta(a+1,b+1)
    logs<-sweep(logs,2,lgx);lp[ids]<-apply(logs,1,.np_udens_beta_logsum)-log(n)
    progress('compact shared first moments',max(ids),st$nodes)
  }
  lw<-log(st$prob)+lp;normal<-.np_udens_beta_logsum(lw);pp<-exp(lw-normal);pp<-pp/sum(pp)
  pilot.mean<-guide.mean<-numeric(m)
  for(ids in blocks) {
    R<-.np_udens_compact_operator(plan,exp(st$lv[ids]))-1-mu-B%*%Ugrid[,ids,drop=FALSE]
    pilot.mean<-pilot.mean+drop(R%*%pp[ids]);guide.mean<-guide.mean+drop(R%*%st$prob[ids])
  }
  if(max(abs(guide.mean))>1e-9)stop('compact plot guide residual is not centered',call.=FALSE)
  covariance<-guide.covariance<-if(full)matrix(0,m,m) else NULL
  pv<-gv<-numeric(m)
  for(ids in blocks) {
    R<-.np_udens_compact_operator(plan,exp(st$lv[ids]))-1-mu-B%*%Ugrid[,ids,drop=FALSE]
    W<-sweep(R-pilot.mean,2,sqrt(pp[ids]),'*');G<-sweep(R-guide.mean,2,sqrt(st$prob[ids]),'*')
    pv<-pv+rowSums(W*W)/n;gv<-gv+rowSums(G*G)/n
    if(full) { covariance<-covariance+tcrossprod(W)/n;guide.covariance<-guide.covariance+tcrossprod(G)/n }
    progress('compact shared covariance',max(ids),st$nodes)
  }
  # The coherent empirical Gram is mandatory, including tiny residual scales.
  canonical<-.np_udens_beta_marginal_rule(theta,64L)
  U<-t(solve(canonical$I,t(sweep(Td,2,colMeans(Td)))))
  Ueta<-sweep(U,2,theta+1,'/');Veta<-crossprod(Ueta)/(n*(n-1))
  observed<-matrix(0,m,n)
  for(ids in .np_udens_beta_blocks(n,tile))observed[,ids]<-.np_udens_compact_operator(plan,x[ids],
    .np_udens_compact_geometry_slice(physical,observations=ids))-1
  rows<-observed-rowMeans(observed)-B%*%t(U)
  ev<-rowSums(rows*rows)/(n*(n-1));empirical<-if(full)tcrossprod(rows)/(n*(n-1)) else NULL
  if(any(!is.finite(c(pv,ev)))||any(pv<=0)||any(ev<=0))stop('compact plot residual scale is not positive',call.=FALSE)
  fitted<-exp(drop(cbind(log(plan$target),log1p(-plan$target))%*%theta)-st$A)
  if(full)for(cv in list(covariance,empirical,guide.covariance))
    if(min(eigen(cv,symmetric=TRUE,only.values=TRUE)$values)< -1e-10)stop('compact plot covariance PSD',call.=FALSE)
  list(center=rowMeans(observed)-mu+fitted,fitted=fitted,mu=mu,B=B,observed=observed,
    pilot_variance=pv,empirical_variance=ev,guide_variance=gv,covariance=covariance,empirical=empirical,
    guide_covariance=guide.covariance,Veta=Veta,ratio_normalization=exp(normal),
    guide_residual_mean=guide.mean,theta=theta,nodes=st$nodes,tail=st$tail)
}

.np_udens_compact_prepare <- function(x,h,grid,kernel,full=TRUE,progress=function(...)NULL,physical=NULL) {
  n<-length(x);m<-length(grid);plan<-.np_udens_compact_operator_plan(grid,h,kernel)
  # Charge observed/influence copies, both refinement results and factor scratch.
  bytes<-8*(8*as.double(m)*n+12*as.double(m)^2+12*512*max(m,n))+64*1024^2
  if(!is.finite(bytes)||bytes>256*1024^2)stop('compact plot inference workspace forecast exceeds 256 MiB',call.=FALSE)
  progress('compact guide',0L,1L);fit<-.np_udens_beta_q1_fit(x);progress('compact guide',1L,1L)
  a<-.np_udens_compact_moments(x,plan,fit$theta,128L,full,progress,physical)
  b<-.np_udens_compact_moments(x,plan,fit$theta,256L,full,progress,physical)
  fields<-c('center','fitted','mu','B','observed','pilot_variance','empirical_variance','guide_variance','Veta','ratio_normalization')
  if(full)fields<-c(fields,'covariance','empirical','guide_covariance')
  a$checks<-setNames(lapply(fields,function(k).np_udens_beta_check(a[[k]],b[[k]],paste('compact',k),if(k=='mu')1e-10 else 1e-9)),fields)
  a$checks$workspace<-list(bytes=bytes,cap=256*1024^2,main.nodes=a$nodes,refined.nodes=b$nodes)
  a$history<-fit$history;a
}
