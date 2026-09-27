test_that("scalar conditional ANN trees preserve literal deleted-fold criteria", {
 old<-options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,np.extendednn=FALSE)
 on.exit(options(old),add=TRUE)
 ns<-environment(npcdensbw)
 for(n in c(63L,64L,65L))for(kernel in c("epanechnikov","uniform")){
  set.seed(833+n);x<-data.frame(x=runif(n,-1,1));y<-data.frame(y=sin(x$x)+rnorm(n))
  x$x[2]<-x$x[1];y$y[3]<-y$y[4]
  if(n==64L)x$u<-factor(rep(c("a","b"),length.out=n))
  perm<-sample.int(n);x<-x[perm,,drop=FALSE];y<-y[perm,,drop=FALSE]
  k<-floor(.6*n);grid<-data.frame(y=c(-.4,.3,.8,.3))
  ml<-ls<-cdf<-0
  for(i in seq_len(n)){
   donors<-setdiff(seq_len(n),i)
   radius<-function(z)vapply(donors,function(j)sort(abs(z[-c(i,j)]-z[j]))[k],0.0)
   hx<-radius(x$x);hy<-radius(y$y);u<-(x$x[i]-x$x[donors])/hx
   w<-if(kernel=="uniform").5*(abs(u)<1)/hx else 3/(4*sqrt(5))*pmax(0,1-u*u/5)/hx
   if(ncol(x)>1)w<-w*ifelse(x$u[donors]==x$u[i],.8,.2)
   w<-w/sum(w);f<-sum(w*dnorm(y$y[i],y$y[donors],hy))
   ml<-ml-log(f)
   sd<-sqrt(outer(hy^2,hy^2,"+"))
   conv<-dnorm(outer(y$y[donors],y$y[donors],"-")/sd)/sd
   ls<-ls+sum(outer(w,w)*conv)-2*f
   cdf<-cdf+sum(vapply(grid$y,function(q)
     ((y$y[i]<=q)-sum(w*pnorm((q-y$y[donors])/hy)))^2,0.0))
  }
  # The eval-only density interface returns the maximization orientation
  # (log likelihood and negative CVLS); CDF returns its squared loss.
  expected<-c("cv.ml"=-ml,"cv.ls"=-ls/n,"cdf"=cdf/(n*nrow(grid)))
  for(family in names(expected))for(engine in c("lc","lp"))for(tree in list(FALSE,TRUE,"auto")){
   options(np.tree=tree)
   ctor<-get(if(family=="cdf")"npcdistbw"else"npcdensbw",ns)
   evaluator<-get(if(family=="cdf")".npcdistbw_eval_only"else".npcdensbw_eval_only",ns)
   args<-list(xdat=x,ydat=y,bws=c(k,k,if(ncol(x)>1).2),bwtype="adaptive_nn",
     bwmethod=if(family=="cv.ml")"cv.ml"else"cv.ls",regtype=engine,
     cxkertype=kernel,cykertype="gaussian",uxkertype="aitchisonaitken",bandwidth.compute=FALSE)
   if(engine=="lp")args$degree<-0L
   b<-do.call(ctor,args)
   call<-list(xdat=x,ydat=y,bws=b)
   if(family=="cdf"){call$gydat<-grid;call$do.full.integral<-FALSE}
   value<-do.call(evaluator,call)$objective
   expect_true(abs(value-expected[[family]])<2e-10,
     info=paste(n,kernel,family,engine,tree,format(value,digits=17),format(expected[[family]],digits=17)))
  }
 }
})
