# Invocation-owned preparation is shared by target blocks and all required laws.
.np_udens_beta_q2_moment_prepare <- function(x,h,grid,theta,inner,outer,state=NULL) {
  n<-nrow(x);st<-if(is.null(state)).np_udens_beta_joint_rule(theta,32L) else state
  axes<-list(.np_udens_beta_marginal_rule(theta[1:2],outer),.np_udens_beta_marginal_rule(theta[3:4],outer))
  nodes<-vapply(axes,`[[`,integer(1),'nodes')
  if(max(nodes)>1600L||prod(nodes)>2560000L)stop("beta plot joint refinement node cap",call.=FALSE)
  cache<-.np_udens_beta_axis_cache(grid,axes,h,n^(-13/54),inner)
  Td<-.np_udens_beta_basis(x);lg<-drop(Td%*%theta)-st$A;rho<-n^(-1/3)
  mats<-vector('list',2);shifts<-numeric(2)
  for(j in 1:2) {
    v<-exp(axes[[j]]$lv);a<-v/rho;b<-(1-v)/rho
    L<-outer(a,log(x[,j]))+outer(b,log1p(-x[,j]))-lbeta(a+1,b+1)
    L<-sweep(L,2,.5*lg);shifts[j]<-max(L);mats[[j]]<-exp(L-shifts[j])
  }
  ratio<-tcrossprod(mats[[1]],mats[[2]])/n
  if(any(!is.finite(ratio))||any(ratio<=0))stop("beta plot joint ratio cache range",call.=FALSE)
  ratio<-log(ratio)+sum(shifts)
  list(x=x,h=h,grid=grid,theta=theta,inner=inner,outer=outer,state=st,axes=axes,
       cache=cache,ratio=ratio,Td=Td,Iinv=solve(st$I),
       logz=st$A-axes[[1]]$A-axes[[2]]$A)
}

.np_udens_beta_q2_moment_block <- function(prep,rows,full=FALSE) {
  if(!is.integer(rows)||!length(rows)||anyNA(rows)||any(rows<1L|rows>nrow(prep$grid))||
     (full&&length(rows)>128L))stop("beta plot moment row/block budget",call.=FALSE)
  n<-nrow(prep$x);m<-length(rows);st<-prep$state;theta<-prep$theta
  cache<-prep$cache;cache$index<-lapply(cache$index,function(ix)ix[rows])
  raw<-.np_udens_beta_joint_raw_native(cache,prep$axes,theta[5],prep$logz,prep$ratio,full)
  guide<-raw[[1]];score<-m+seq_len(5L);ix<-seq_len(m)
  if(abs(exp(guide$logmass)-1)>1e-11)stop("beta plot joint guide mass refinement",call.=FALSE)
  .np_udens_beta_check(guide$mean[score],st$mean,"joint score mean",1e-10)
  .np_udens_beta_check(guide$cross[score,,drop=FALSE],st$I,"joint information")
  mu<-guide$mean[ix];B<-guide$cross[ix,,drop=FALSE];K<-B%*%prep$Iinv
  x<-prep$x;grid<-prep$grid[rows,,drop=FALSE]
  da<-list(list(lv=log(x[,1]),l1v=log1p(-x[,1])),list(lv=log(x[,2]),l1v=log1p(-x[,2])))
  obs<-.np_udens_beta_joint_block_native(.np_udens_beta_axis_cache(grid,da,prep$h,n^(-13/54),prep$inner),seq_len(n),seq_len(n))
  Td<-sweep(prep$Td,2,colMeans(prep$Td));residual<-obs-rowMeans(obs)-K%*%t(Td)
  ve<-rowSums(residual^2)/(n*(n-1));emp<-if(full)tcrossprod(residual)/(n*(n-1)) else NULL
  covs<-variances<-identities<-list();Tg<-.np_udens_beta_basis(grid)
  f<-exp(drop(Tg%*%theta)-st$A);G<-sweep(sweep(Tg,2,st$mean),1,f,'*')%*%prep$Iinv
  for(law in 1:2) {
    name<-c('guide','ratio')[law];r<-raw[[law]];cross<-r$cross[ix,,drop=FALSE];V<-r$cross[score,,drop=FALSE]
    variance<-(r$variance[ix]-2*rowSums(K*cross)+rowSums((K%*%V)*K))/n
    if(any(!is.finite(variance))||any(variance<=0))stop("beta plot residual scale is not positive",call.=FALSE)
    variances[[name]]<-variance
    # The full corrected influence identity is checked without a dense target matrix.
    full.direct<-(r$variance[ix]+2*rowSums((G-K)*cross)+rowSums(((G-K)%*%V)*(G-K)))/n
    expanded<-variance+2*rowSums((cross-K%*%V)*G)/n+rowSums((G%*%V)*G)/n
    identities[[name]]<-.np_udens_beta_check(expanded,full.direct,paste('joint influence diagonal',name))
    if(full) {
      covariance<-(r$cov[ix,ix,drop=FALSE]-tcrossprod(K,cross)-tcrossprod(cross,K)+K%*%V%*%t(K))/n
      .np_udens_beta_check(diag(covariance),variance,paste('joint covariance diagonal',name))
      if(min(eigen(covariance,symmetric=TRUE,only.values=TRUE)$values) < -1e-10)
        stop("beta plot joint covariance PSD",call.=FALSE)
      covs[[name]]<-covariance
    }
  }
  U<-Td%*%prep$Iinv;Ueta<-sweep(U,2,c(1/(theta[1:4]+1),1),'*')
  list(center=rowMeans(obs)-mu+f,fitted=f,mu=mu,B=B,observed=obs,empirical=emp,
    covariance=covs$ratio,guide_covariance=covs$guide,empirical_variance=ve,
    pilot_variance=variances$ratio,guide_variance=variances$guide,Veta=crossprod(Ueta)/(n*(n-1)),
    identities=identities,ratio_normalization=exp(raw[[2]]$logmass),guide_normalization=exp(raw[[1]]$logmass))
}
