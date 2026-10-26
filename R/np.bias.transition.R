# Convex continuation through nested weighted polynomial spaces.
# Geometry is fixed while differentiating the query; no training-square operator.
.np_bias_transition_pilot <- function(state, at, bandwidth, value=NULL,
                     derivatives=integer(), orders=NULL) {
 x<-state$x;Q<-state$Q;P<-ncol(Q);lambda<-state$geometry$lambda
 g<-drop(.np_bias_adaptation_kernel(x,at,bandwidth))
 jj<-order(lambda,seq_along(lambda)); Z<-Q[,jj,drop=FALSE]*sqrt(g)
 norms<-sqrt(colSums(Z^2))
 if(any(!is.finite(norms))||any(norms<=0))stop('continuous adaptation subspace has an unsupported column',call.=FALSE)
 fit<-qr(sweep(Z,2L,norms,'/'),LAPACK=FALSE,tol=1e-11)
 if(fit$rank!=P||!identical(fit$pivot,seq_len(P)))stop('continuous adaptation nested subspace is unidentified',call.=FALSE)
 U<-qr.Q(fit);R<-qr.R(fit)
 A<-backsolve(R,t(U)/(1+lambda[jj]))*rep(sqrt(g),each=P)
 A<-sweep(A,1L,norms,'/'); sorted<-A; A[jj,]<-sorted
 if(is.null(value))value<-.np_bias_adaptation_basis(x,at,state$geometry$q)
 H<-value%*%A;D<-matrix(0,length(derivatives),nrow(x))
 lev<-sort(unique(lambda)); d<-1/(1+lev); mix<-d-c(d[-1L],0)
 for(k in seq_along(derivatives)) {
  j<-derivatives[k];ord<-if(is.null(orders))1L else orders[j]
  values<-do.call(rbind,lapply(0:ord,function(r)
   .np_bias_adaptation_basis(x,at,state$geometry$q,if(r)j else 0L,r)))
  ratios<-.np_bias_gaussian_derivatives((x[,j]-at[1,j])/bandwidth[j],bandwidth[j],ord)
  for(l in seq_along(lev)) {
   cols<-which(lambda<=lev[l])
   D[k,]<-D[k,]+mix[l]*.np_bias_weighted_derivative(Q[,cols,drop=FALSE],g,
    values[,cols,drop=FALSE],ratios,role='continuous adaptation subspace derivative',rank.tolerance=1e-11)
  }
 }
 list(A=A,g=g,H=H,D=D)
}
