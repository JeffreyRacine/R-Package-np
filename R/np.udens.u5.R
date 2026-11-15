# Plot-only U5 support-reproducing correction; native estimators are unchanged.
.np_udens_u5_kernel <- function(name) {
  if (!identical(name, "epanechnikov")) stop("U5 requires Epanechnikov", call.=FALSE)
  list(span=sqrt(5))
}
# Private research-only support-adapted polynomial operator.
.np_udens_u5_poly<-function(z,d=5L) {
 P<-matrix(1,length(z),d+1L)
 if(d>0)P[,2]<-z
 if(d>1)for(j in 2:d)P[,j+1]<-((2*j-1)*z*P[,j]-(j-1)*P[,j-1])/j
 P
}
.np_udens_u5_boundary<-function(name,l,r,order=16L,kernel.order=4L) {
 if (!isTRUE(kernel.order %in% c(4L,6L,8L))) stop("U5 requires kernel order 4, 6 or 8",call.=FALSE)
 degree <- as.integer(kernel.order+1L)
 k<-.np_udens_u5_kernel(name);stopifnot(is.finite(l),is.finite(r),l<r,l<=0,r>=0,l>=-k$span,r<=k$span)
 half<-(r-l)/2;mid<-(r+l)/2;rule<-.np_udens_q6_rule(order);u<-mid+half*rule$x
 weight<-function(v)if(name=='uniform')rep(.5,length(v)) else 3/(4*sqrt(5))*(1-v^2/5)
 P<-.np_udens_u5_poly(rule$x,degree);w<-half*rule$w*weight(u);G<-crossprod(P,P*w);p0<-drop(.np_udens_u5_poly(-mid/half,degree));coef<-solve(G,p0)
 L<-function(v) {keep<-v>=l & v<=r;out<-numeric(length(v));out[keep]<-weight(v[keep])*drop(.np_udens_u5_poly((v[keep]-mid)/half,degree)%*%coef);out}
 oracle<-function(v)drop(.np_udens_u5_poly((v-mid)/half,degree)%*%((2*(0:degree)+1)/2*p0))/half
 list(L=L,analytic_uniform=oracle,condition=kappa(G,exact=TRUE),l=l,r=r)
}

# Complete covariance arithmetic copied exactly from the qualified U5 control.
.np_udens_u5_covariance<-function(z,x,kn,grid,order,kernel.order=4L,gram=FALSE) {
 h<-as.double(z$h);n<-length(x);j<-0:z$k;A<-z$a+j;BB<-z$b+z$k-j;weights<-rowMeans(vapply(x,function(t)dbinom(j,z$k,t),numeric(length(j))));abq<-c(z$aQ,z$bQ)
 logmeans<-cbind(digamma(A)-digamma(A+BB),digamma(BB)-digamma(A+BB));m<-drop(crossprod(weights,logmeans));I<-diag(trigamma(abq))-trigamma(sum(abq));TT<-matrix(0,2,2)
 for(k in seq_along(weights))TT<-TT+weights[k]*(diag(trigamma(c(A[k],BB[k])))-trigamma(A[k]+BB[k])+tcrossprod(logmeans[k,]));TT<-TT-tcrossprod(m)
 span<-.np_udens_u5_kernel(kn)$span;ops<-lapply(grid,function(t).np_udens_u5_boundary(kn,max(-span,-t/h),min(span,(1-t)/h),kernel.order=kernel.order))
 limits<-c(-(100+max(abs(lbeta(A,BB))))/min(c(A,abq[1])),(100+max(abs(lbeta(A,BB))))/min(c(BB,abq[2])))
 physical<-sort(unique(c(grid,grid-span*h,grid+span*h)));physical<-physical[physical>0&physical<1];fixed<-c(-80,-40,-20,-10,-5,0,5,10,20,40,80)
 cuts<-sort(unique(c(limits,fixed[fixed>limits[1]&fixed<limits[2]],qlogis(physical))));stopifnot(cuts[1]==limits[1],tail(cuts,1)==limits[2])
  calc<-function(order) {
   gl<-.np_udens_q6_rule(order);t<-w<-numeric()
   for(k in seq_len(length(cuts)-1L)){mid<-(cuts[k]+cuts[k+1])/2;half<-diff(cuts[k+0:1])/2;t<-c(t,mid+half*gl$x);w<-c(w,half*gl$w)}
   softplus<-function(t)pmax(t,0)+log1p(exp(-abs(t)));lv<--softplus(-t);l1v<--softplus(t);T<-cbind(lv,l1v);u<-exp(lv)
   component<-exp(sweep(outer(lv,A)+outer(l1v,BB),2,lbeta(A,BB),'-'));qprob<-drop(component%*%weights)*w
   gprob<-exp(abq[1]*lv+abq[2]*l1v-lbeta(abq[1],abq[2]))*w
   AX<-t(vapply(seq_along(grid),function(i)ops[[i]]$L((u-grid[i])/h)/h,numeric(length(u))))
   EA<-drop(AX%*%qprob);BG<-AX%*%(sweep(T,2,m)*gprob);gg<-dbeta(grid,abq[1],abq[2]);G<-sweep(sweep(cbind(log(grid),log1p(-grid)),2,m),1,gg,'*');H<-G-BG
   U<-solve(I,t(sweep(T,2,m)));psi<-AX-EA+H%*%U;centered<-psi-drop(psi%*%qprob)
   direct<-tcrossprod(sweep(centered,2,sqrt(qprob),'*'))/n
   AA<-tcrossprod(sweep(AX,2,sqrt(qprob),'*'))-tcrossprod(EA);AT<-AX%*%(sweep(T,2,m)*qprob);HI<-H%*%solve(I)
   expanded<-(AA+AT%*%t(HI)+HI%*%t(AT)+HI%*%TT%*%t(HI))/n
   means<-t(AX%*%(component*w));CQ<-EA-drop(AX%*%gprob)+gg
   factor <- if(gram) .np_udens_u5_gram_factor(sweep(centered,2,sqrt(qprob/n),'*')) else NULL
   list(factor=factor,C=direct,expanded=expanded,center=CQ,EA=EA,B=BG,H=H,HI=HI,component_means=means,AX=AX,T=T,lv=lv,l1v=l1v,w=w,nodes=length(t),mass=sum(qprob),mean=drop(crossprod(qprob,T)),psi_mean=drop(psi%*%qprob),residual_only=AA/n)
  }
 calc(order)
}

.np_udens_u5_check<-function(a,b,label){e<-max(abs(a-b));lim<-1e-12+1e-9*max(abs(b));if(!is.finite(e)||e>lim)stop(paste(label,e,'>',lim));e}
.np_udens_u5_guide<-function(grid,h,ab,kn,budget,relative,kernel.order=4L){
 span<-.np_udens_u5_kernel(kn)$span;low<--(budget+abs(lbeta(ab[1],ab[2])))/ab[1];high<-(budget+abs(lbeta(ab[1],ab[2])))/ab[2]
 softplus<-function(t)pmax(t,0)+log1p(exp(-abs(t)));values<-errors<-numeric(length(grid))
 for(i in seq_along(grid)){
  g<-grid[i];op<-.np_udens_u5_boundary(kn,max(-span,-g/h),min(span,(1-g)/h),kernel.order=kernel.order);phys<-c(g,g-span*h,g+span*h);phys<-phys[phys>0&phys<1];cut<-sort(unique(c(low,0,qlogis(phys),high)));stopifnot(cut[1]==low,tail(cut,1)==high)
  fun<-function(t){lv<--softplus(-t);lu<--softplus(t);op$L((exp(lv)-g)/h)/h*exp(ab[1]*lv+ab[2]*lu-lbeta(ab[1],ab[2]))}
  for(k in seq_len(length(cut)-1L)){v<-integrate(fun,cut[k],cut[k+1L],rel.tol=relative,abs.tol=1e-12,subdivisions=1000L);values[i]<-values[i]+v$value;errors[i]<-errors[i]+v$abs.error}
 }
 list(value=values,error=errors)
}

.np_udens_u5_fitmom <- function(target, start) {
  ab <- start
  for (it in 1:80) {
    grad <- digamma(ab)-digamma(sum(ab))-target
    I <- diag(trigamma(ab))-trigamma(sum(ab)); direction <- solve(I,grad)
    if (max(abs(grad)) < 1e-12) return(list(ab=ab,score=max(abs(grad)),iterations=it))
    objective <- function(z) lbeta(z[1],z[2])-sum((z-1)*target)
    step <- 1
    while (any(ab-step*direction<=0) ||
      objective(ab-step*direction)>objective(ab)-1e-4*step*sum(grad*direction)+1e-14) {
      step <- step/2
      if (step<=2^-40) stop("U5 beta log-moment line search incomplete",call.=FALSE)
    }
    ab <- ab-step*direction
  }
  stop("U5 beta exact-log-moment fit incomplete",call.=FALSE)
}

.np_udens_u5_budget <- function(n, m, block=8L) {
  if (any(!is.finite(c(n,m,block))) || any(c(n,m,block)!=floor(c(n,m,block))) ||
      n<10 || m<1 || block<1 || any(c(n,m,block)>.Machine$integer.max))
    stop("invalid U5 workspace dimensions",call.=FALSE)
  k <- ceiling(2*n^(2/5)); targets <- min(m,block)
  nodes <- 96*(3*targets+15)
  bytes <- 8*(4*n*(k+1)+8*nodes*(k+1)+24*nodes*targets+32*nodes+12*m+16*targets^2)
  if (!is.finite(bytes) || bytes>192*1024^2)
    stop("U5 inference exceeds its 192 MiB workspace bound",call.=FALSE)
  bytes
}

.np_udens_u5_moments <- function(x,h,grid,progress=function(...)NULL,block=8L,kernel.order=4L,joint=FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || any(!is.finite(x)) ||
      any(x<=0|x>=1) || !is.numeric(grid) || !length(grid) ||
      any(!is.finite(grid)) || any(grid<=0|grid>=1) ||
      !is.numeric(h) || length(h)!=1L || !is.finite(h) || h<=0)
    stop("U5 requires finite strictly interior data and targets and positive bandwidth",call.=FALSE)
  workspace <- .np_udens_u5_budget(length(x),length(grid),block)
  targets <- sort(unique(grid)); index <- match(grid,targets)
  pilot <- .np_udens_beta_q1_fit(x); ab <- pilot$theta+1
  if (max(abs(digamma(ab)-digamma(sum(ab))-c(mean(log(x)),mean(log1p(-x))))) >= 2e-8)
    stop("U5 original guide score exceeds its frozen tolerance",call.=FALSE)
  n <- length(x);rho <- n^(-2/5);k <- ceiling(2/rho);jj <- 0:k
  A <- ab[1]+jj;BB <- ab[2]+k-jj
  weights <- rowMeans(vapply(x,function(t)dbinom(jj,k,t),numeric(k+1)))
  if (any(!is.finite(weights)) || any(weights<0)) stop("invalid U5 reference weights",call.=FALSE)
  .np_udens_u5_check(sum(weights),1,"reference mass")
  logQ <- c(sum(weights*(digamma(A)-digamma(sum(ab)+k))),sum(weights*(digamma(BB)-digamma(sum(ab)+k))))
  fq <- .np_udens_u5_fitmom(logQ,ab)
  z <- list(k=k,a=ab[1],b=ab[2],aQ=fq$ab[1],bQ=fq$ab[2],h=h)
  center <- se <- numeric(length(targets)); checks <- list()
  for(first in seq.int(1L,length(targets),by=block)) {
    ix <- first:min(first+block-1L,length(targets)); g <- targets[ix]
    low <- .np_udens_u5_covariance(z,x,"epanechnikov",g,64L,kernel.order)
    ref <- .np_udens_u5_covariance(z,x,"epanechnikov",g,96L,kernel.order)
    errors <- vapply(c("C","center","EA","B","component_means"),function(key)
      .np_udens_u5_check(low[[key]],ref[[key]],paste("refinement",key)),0.0)
    errors <- c(errors,.np_udens_u5_check(ref$mass,1,"mass"),
      .np_udens_u5_check(ref$mean,digamma(fq$ab)-digamma(sum(fq$ab)),"Q log moments"),
      .np_udens_u5_check(ref$C,ref$expanded,"complete covariance"),
      .np_udens_u5_check(ref$psi_mean,rep(0,length(g)),"complete mean"))
    if (min(eigen(ref$C,symmetric=TRUE,only.values=TRUE)$values) < -1e-10 ||
        any(!is.finite(diag(ref$C))) || any(diag(ref$C)<=0))
      stop("U5 complete covariance is not positive finite",call.=FALSE)
    g1 <- .np_udens_u5_guide(g,h,ab,"epanechnikov",80,1e-10,kernel.order)
    g2 <- .np_udens_u5_guide(g,h,ab,"epanechnikov",100,2e-11,kernel.order)
    errors <- c(errors,.np_udens_u5_check(g1$value,g2$value,"original guide refinement"))
    span <- sqrt(5)
    observed <- vapply(g,function(t) {
      op <- .np_udens_u5_boundary("epanechnikov",max(-span,-t/h),min(span,(1-t)/h),kernel.order=kernel.order)
      mean(op$L((x-t)/h)/h)
    },0.0)
    center[ix] <- observed-g2$value+dbeta(g,ab[1],ab[2]);se[ix] <- sqrt(diag(ref$C))
    checks[[length(checks)+1L]] <- errors
    progress("complete marginal moments",max(ix),length(targets))
  }
  if (any(!is.finite(c(center,se)))) stop("nonfinite U5 moments",call.=FALSE)
  out <- list(center=center[index],se=se[index],checks=checks,workspace=workspace,
       original.guide=ab,reference.guide=fq$ab,reference.degree=k)
  if(joint) {
    .np_udens_u5_joint_budget(n,length(targets),1L)
    low <- .np_udens_u5_covariance(z,x,"epanechnikov",targets,64L,kernel.order)
    high <- .np_udens_u5_covariance(z,x,"epanechnikov",targets,96L,kernel.order,gram=TRUE)
    out$joint.checks <- c(.np_udens_u5_check(low$C,high$C,"joint refinement"),
      .np_udens_u5_check(high$C,high$expanded,"joint expanded covariance"),
      .np_udens_u5_check(tcrossprod(high$factor),high$C,"joint factor covariance"),
      .np_udens_u5_check(se^2,diag(high$C),"joint marginal covariance"))
    out$factor <- high$factor;out$joint.se <- se
    progress("complete joint covariance",length(targets),length(targets))
  }
  out
}

.np_udens_u5_request <- function(bws,control,center,method,band,calibration=NULL) {
  inherits(bws,"bandwidth") && is.null(control) && is.null(calibration) &&
    identical(center,"bias-corrected") && identical(method,"asymptotic") && isTRUE(band %in% c("pmzsd","pointwise","bonferroni","simultaneous","all")) &&
    identical(bws[["type",exact=TRUE]],"fixed") && identical(bws[["ckerbound",exact=TRUE]],"fixed") &&
    identical(bws[["ckertype",exact=TRUE]],"epanechnikov") &&
    isTRUE(bws[["ckerorder",exact=TRUE]] %in% c(4L,6L,8L)) && isTRUE(bws[["ndim",exact=TRUE]]==1L) &&
    isTRUE(bws[["ncon",exact=TRUE]]==1L)
}
.np_udens_u5_active <- function(state)
  isTRUE(state[["construction",exact=TRUE]] %in% paste0("density-u5-order",c(4L,6L,8L),"-normal-v1"))

.np_udens_u5_prepare <- function(bws,x,random.seed=42L) {
  npValidateNonNegativeInteger(random.seed,"random.seed")
  if (!.np_udens_u5_request(bws,NULL,"bias-corrected","asymptotic","pmzsd") ||
      !is.data.frame(x) || ncol(x)!=1L || !nzchar(names(x)[1L]) ||
      !is.numeric(x[[1L]]) || is.object(x[[1L]]) || any(!is.finite(x[[1L]])))
    stop("invalid U5 fixed-bound higher-order data",call.=FALSE)
  lower <- bws[["ckerlb",exact=TRUE]];upper <- bws[["ckerub",exact=TRUE]]
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower)!=1L || length(upper)!=1L ||
      any(!is.finite(c(lower,upper))) || upper<=lower)
    stop("U5 requires finite ordered fixed bounds",call.=FALSE)
  width <- upper-lower;h <- as.double(.np_physical_bandwidth(bws))/width
  if (!is.finite(width) || !is.finite(1/width) || width<=0 || length(h)!=1L || !is.finite(h) || h<=0)
    stop("U5 support units or bandwidth are not representable",call.=FALSE)
  normalized <- (x[[1L]]-lower)/width
  if (any(!is.finite(normalized)) || any(normalized<=0|normalized>=1))
    stop("U5 requires strictly interior observations; support endpoints are excluded",call.=FALSE)
  .np_udens_u5_budget(nrow(x),1L)
  list(x=x,normalized=normalized,h=h,lower=lower,width=width,n=nrow(x),q=1L,
       normal.pmzsd=TRUE,plot.law="gaussian",plot.random.seed=random.seed,kernel.order=as.integer(bws[["ckerorder",exact=TRUE]]),
       construction=paste0("density-u5-order",bws[["ckerorder",exact=TRUE]],"-normal-v1"))
}

# Complete positive-weight influence factor: retain all represented directions.
.np_udens_u5_gram_factor <- function(D) {
  if (!is.matrix(D) || !is.numeric(D) || !nrow(D) || ncol(D)<nrow(D) || any(!is.finite(D)))
    stop("invalid U5 complete influence factor",call.=FALSE)
  z <- qr(t(D),LAPACK=TRUE); R <- qr.R(z,complete=FALSE); m <- nrow(D)
  if (!identical(dim(R),c(m,m)) || !identical(sort(z$pivot),seq_len(m)))
    stop("incomplete U5 full QR factor",call.=FALSE)
  F <- matrix(0,m,m); F[z$pivot,] <- t(R)
  .np_udens_u5_check(tcrossprod(F),tcrossprod(D),"complete Gram reconstruction")
  F
}

.np_udens_u5_joint_budget <- function(n,m,B) {
  bytes <- .np_udens_u5_budget(n,m,block=m)+8*(6*as.double(m)^2+256*m+3*B)
  if (!is.finite(bytes) || bytes>192*1024^2)
    stop("U5 joint inference exceeds its 192 MiB workspace bound",call.=FALSE)
  bytes
}

.np_udens_u5_public <- function(state,grid,alpha,band,B=9999L) {
  if (!.np_udens_u5_active(state) || !isTRUE(band %in% c("pmzsd","pointwise","bonferroni","simultaneous","all")) ||
      !is.data.frame(grid) || ncol(grid)!=1L || !nrow(grid) ||
      !identical(names(grid),names(state$x)) || !is.numeric(grid[[1L]]) ||
      any(!is.finite(grid[[1L]])) || !is.numeric(alpha) || length(alpha)!=1L ||
      !is.finite(alpha) || alpha<=0 || alpha>=.5) stop("invalid U5 band request",call.=FALSE)
  target <- (grid[[1L]]-state$lower)/state$width
  if (any(!is.finite(target)) || any(target<=0|target>=1))
    stop("U5 bands require strictly interior evaluation points; support endpoints are excluded",call.=FALSE)
  joint <- band %in% c("simultaneous","all");m <- length(unique(target))
  if (joint) { B <- .np_plot_gaussian_integer(B,"gaussian.draws",.Machine$integer.max)
    .np_udens_u5_joint_budget(state$n,m,B);.np_plot_gaussian_draws(B,alpha,m) }
  progress <- function(stage,done,total) .np_plot_progress_notify(stage=paste("U5 density:",stage),done=done,total=total,unit="work")
  z <- .np_plot_corrected_activity_run("Constructing higher-order density normal intervals",
    .np_udens_u5_moments(state$normalized,state$h,target,progress,kernel.order=state$kernel.order,joint=joint))
  center <- z$center/state$width;se <- z$se/state$width
  critical <- c(pointwise=qnorm(alpha/2,lower.tail=FALSE),bonferroni=qnorm(alpha/(2*m),lower.tail=FALSE))
  maxima <- NULL
  if (joint) {
    maxima <- .np_with_seed(state$plot.random.seed, {
      ans <- numeric(B)
      for(first in seq.int(1L,B,by=128L)) {
        ix <- first:min(first+127L,B);Z <- matrix(rnorm(m*length(ix)),m)
        ans[ix] <- apply(abs((z$factor %*% Z)/z$joint.se),2L,max)
        progress("Gaussian maximum calibration",max(ix),B)
      }
      ans
    })
    if(any(!is.finite(maxima)))stop("nonfinite U5 Gaussian maximum",call.=FALSE)
    critical <- c(critical,simultaneous=unname(quantile(maxima,1-alpha,type=7L)))
  }
  requested <- if(band=="all")c("pointwise","simultaneous","bonferroni") else band
  bands <- setNames(lapply(requested,function(k) {
    width <- critical[[if(k=="pmzsd")"pointwise" else k]]*se
    array(c(center-width,center+width),c(nrow(grid),2L,1L),dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
  }),requested)
  inference <- list(center=center,se=se,bands=bands,intervals=bands,absolute.intervals=bands,
    critical=critical,maximum=maxima,covariance=NULL,
    settings=list(construction=state$construction,calibration=if(joint)"Gaussian maximum and normal quantiles" else "normal-quantile",law="asymptotic",
      kernel="epanechnikov",kernel.order=state$kernel.order,q=1L,n=state$n,alpha=alpha,B=if(joint)B else 0L,available.bands=requested,
      grid.size=m,family="unique supplied interior targets",family.rows=match(target,sort(unique(target))),
      bounds=c(state$lower,state$lower+state$width),bandwidth=state$h*state$width,primary="absolute",quantile.type=7L,
      correction=paste0("degree-",state$kernel.order+1L," support-reproducing operator with fitted beta guide; native fit retained"),
      bonferroni.construction="normal quantile at alpha divided by unique family size; no simulated hull",
      calibration.status="complete positive-reference covariance with guide and cross terms"),
    reference=list(original.guide=z$original.guide,variance.guide=z$reference.guide,degree=z$reference.degree),
    numerical.checks=list(refinement=z$checks,workspace=z$workspace,joint=z$joint.checks))
  if(band=="pmzsd") inference <- .np_plot_bias_normal(inference,alpha)
  widths <- lapply(bands,function(v)cbind(center-v[,1L,1L],v[,2L,1L]-center))
  list(boot.err=unname(cbind(widths[[if(band=="all")"pointwise" else band]],center)),
       boot.all.err=if(band=="all")widths else NULL,inference=inference,bxp=list())
}
