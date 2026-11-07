.np_udens_beta_joint_normalizer <- function(eta,order=16L) {
  if(is.null(dim(eta)))eta<-matrix(eta,1L)
  .Call(C_np_plot_udens_beta_normalizer,matrix(as.double(eta),nrow(eta)),as.integer(order))
}


.np_udens_beta_sphere_accept <- function(value,trial,predicted,newnorm,oldnorm) {
  if(trial<=value+1e-4*predicted)return(TRUE)
  roundoff<-32*.Machine$double.eps*(1+abs(value))
  abs(predicted)<=roundoff && abs(trial-value)<=roundoff && newnorm<=.5*oldnorm
}

.np_udens_beta_sphere <- function(objective,start,maxiter=150L,gtol=1e-7) {
  norm<-function(z)sqrt(sum(z*z));z<-start/norm(start);H<-diag(length(z));history<-list()
  out<-objective(z);value<-out$value;g<-out$gradient;tangent<-g-sum(g*z)*z
  # Frozen five-step BFGS phase, followed by the same bounded Newton polishing.
  for(it in seq_len(min(5L,maxiter))) {
    gn<-norm(tangent);history[[length(history)+1L]]<-c(value=value,gradient=gn)
    if(gn<=gtol)return(list(x=z,iterations=it-1L,history=history))
    direction<--drop(H%*%tangent);direction<-direction-z*sum(z*direction)
    length<-norm(direction);if(length>1)direction<-direction/length
    descent<-sum(tangent*direction)
    if(!is.finite(descent)||descent>=0)stop("beta plot sphere BFGS descent",call.=FALSE)
    accepted<-FALSE
    for(line in 0:29) {
      step<-direction*.5^line;length<-norm(step);candidate<-cos(length)*z+sin(length)*step/length
      candidate<-candidate/norm(candidate);trial<-objective(candidate)
      newtangent<-trial$gradient-sum(trial$gradient*candidate)*candidate
      if(.np_udens_beta_sphere_accept(value,trial$value,.5^line*descent,norm(newtangent),gn)) {accepted<-TRUE;break}
    }
    if(!accepted)stop("beta plot sphere BFGS line-search",call.=FALSE)
    K<-tcrossprod(candidate,z)-tcrossprod(z,candidate)
    rotation<-diag(5)+K+K%*%K/(1+sum(candidate*z))
    ss<-drop(rotation%*%step);yy<-newtangent-drop(rotation%*%tangent);H<-rotation%*%H%*%t(rotation)
    sy<-sum(ss*yy)
    if(sy>1e-10*norm(ss)*norm(yy)) {
      A<-diag(5)-tcrossprod(ss,yy)/sy;H<-A%*%H%*%t(A)+tcrossprod(ss)/sy;H<-(H+t(H))/2
    }
    z<-candidate;value<-trial$value;tangent<-newtangent
  }
  for(it in seq_len(maxiter-5L)) {
    out<-objective(z);value<-out$value;g<-out$gradient;tangent<-g-sum(g*z)*z;gn<-norm(tangent)
    history[[length(history)+1L]]<-c(value=value,gradient=gn)
    if(gn<=gtol)return(list(x=z,iterations=it-1L+5L,history=history))
    B<-svd(matrix(z,1L),nu=0L,nv=5L)$v[,2:5,drop=FALSE];grad<-drop(crossprod(B,g));H<-matrix(0,4,4);eps<-1e-4
    for(j in 1:4) {
      zp<-cos(eps)*z+sin(eps)*B[,j];zm<-cos(eps)*z-sin(eps)*B[,j]
      gp<-objective(zp)$gradient;gm<-objective(zm)$gradient
      H[,j]<-drop(crossprod(B,(gp-sum(gp*zp)*zp)-(gm-sum(gm*zm)*zm)))/(2*eps)
    }
    H<-(H+t(H))/2;ev<-eigen(H,symmetric=TRUE);floor<-max(1e-6,.01*max(abs(ev$values)))
    step<--drop(ev$vectors%*%(drop(crossprod(ev$vectors,grad))/pmax(ev$values,floor)))
    length<-norm(step);if(length>1)step<-step/length;descent<-sum(grad*step)
    if(!is.finite(descent)||descent>=0)stop("beta plot sphere Newton descent",call.=FALSE)
    accepted<-FALSE
    for(line in 0:29) {
      s<-step*.5^line;length<-norm(s);candidate<-cos(length)*z+sin(length)*drop(B%*%s)/length
      candidate<-candidate/norm(candidate);trial<-objective(candidate)
      newtangent<-trial$gradient-sum(trial$gradient*candidate)*candidate
      if(.np_udens_beta_sphere_accept(value,trial$value,.5^line*descent,norm(newtangent),gn)) {accepted<-TRUE;break}
    }
    if(!accepted)stop("beta plot sphere Newton line-search",call.=FALSE)
    z<-candidate
  }
  stop("beta plot sphere iteration cap",call.=FALSE)
}

.np_udens_beta_region_certificate <- function(runs,index,target,transform,sign) {
  chosen<-runs[[index]]
  # Certify the returned 32-node quantity, independently of the lower-order
  # objective used to locate the retained sphere candidate.
  eta<-do.call(rbind,lapply(runs,`[[`,'eta'))
  normal<-.np_udens_beta_joint_normalizer(eta,32L)
  ref<-.np_udens_beta_joint_rule(chosen$theta,32L,moments=TRUE)
  actual<-sum(chosen$theta*target)-ref$A
  err<-abs(normal[index]-ref$A)
  if(!is.finite(err)||err>1e-10)
    stop("beta plot joint projection independent-normalizer gate",call.=FALSE)
  gradient<-drop(crossprod(transform,c(chosen$theta[1:4]+1,1)*(target-ref$mean)))
  z<-chosen$direction
  kkt<-sqrt(sum((gradient-sum(gradient*z)*z)^2))
  if(!is.finite(kkt)||kkt>1e-4)
    stop("beta plot joint projection stationarity",call.=FALSE)
  values<-vapply(runs,function(run)sum(run$theta*target),numeric(1))-normal
  ranking<-sign*actual-min(sign*values)
  if(!is.finite(ranking)||ranking>1e-10)
    stop("beta plot joint projection refined-ranking gate",call.=FALSE)
  list(value=actual,refinement=err,kkt=kkt,
       optimization_refinement=abs(actual-chosen$value),ranking_gap=ranking)
}

.np_udens_beta_region_q2 <- function(theta,Veta,grid,alpha) {
  eta<-c(log(theta[1:4]+1),theta[5]);L<-t(chol(Veta));radius<-sqrt(qchisq(1-alpha/2,5))
  transform<-radius*L;T<-.np_udens_beta_basis(grid);starts<-rbind(diag(5),-diag(5))
  lower<-upper<-numeric(nrow(grid));proof<-list()
  for(k in seq_len(nrow(grid)))for(sign in c(1,-1)) {
    objective<-function(z) {
      out<-.Call(C_np_plot_udens_beta_region_objective,eta,transform,as.double(T[k,]),as.double(z),as.double(sign))
      names(out)<-c('value','gradient','eta','theta');out
    }
    runs<-vector('list',10L)
    for(j in seq_len(nrow(starts))) {
      r<-.np_udens_beta_sphere(objective,starts[j,]);out<-objective(r$x)
      if(abs(sum(r$x^2)-1)>1e-12)stop("beta plot sphere norm certificate",call.=FALSE)
      kkt<-sqrt(sum((out$gradient-sum(out$gradient*r$x)*r$x)^2))
      if(kkt>1e-4)stop("beta plot joint projection stationarity",call.=FALSE)
      runs[[j]]<-list(value=sign*out$value,theta=out$theta,eta=out$eta,direction=r$x,kkt=kkt,iterations=r$iterations)
    }
    values<-vapply(runs,`[[`,numeric(1),'value');index<-which.min(sign*values)
    certificate<-.np_udens_beta_region_certificate(runs,index,as.double(T[k,]),transform,sign)
    value<-exp(certificate$value);if(!is.finite(value))stop("beta plot nonfinite joint region",call.=FALSE)
    if(sign==1)lower[k]<-value else upper[k]<-value
    proof[[length(proof)+1L]]<-list(target=k,sign=sign,refinement=certificate$refinement,kkt=certificate$kkt,
      optimization_refinement=certificate$optimization_refinement,ranking_gap=certificate$ranking_gap,
      multistart_log_range=diff(range(values)),iterations=sum(vapply(runs,`[[`,numeric(1),'iterations')))
    # Release completed target temporaries before the next bounded region unit.
    gc(FALSE)
  }
  if(any(lower>upper))stop("beta plot unordered joint region",call.=FALSE)
  list(lower=lower,upper=upper,proof=proof)
}
