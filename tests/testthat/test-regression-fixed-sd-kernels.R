test_that("fixed SD admission preserves laws and the local empirical default", {
  for(kernel in c("epanechnikov","beta")) {
    b<-list(type="fixed",ckertype=kernel,ckerorder=2L,ckerbound=if(kernel=="beta")"fixed"else"none",ncon=1L)
    expect_identical(.np_bias_regression_kernel_class(b,"absolute"),"local-envelope")
    a<-list(type="fixed",calibration="absolute")
    for(law in c("gaussian","mammen","rademacher","inid","fixed","geom"))
      expect_silent(.np_bias_regression_kernel_law(a,law))
    expect_error(.np_bias_regression_kernel_law(a,"unknown"),"requires Gaussian, Mammen, Rademacher")
    for(type in c("generalized_nn","adaptive_nn")) {
      z<-b;z$type<-type;expect_identical(.np_bias_regression_kernel_class(z,"absolute"),"unsupported")
    }
    expect_identical(.np_bias_regression_calibration(b,"default","bootstrap","simultaneous","wild","mammen"),"basic-local-envelope")
    expect_identical(.np_bias_regression_calibration(b,"default","asymptotic","simultaneous","wild","mammen"),"absolute")
    expect_identical(.np_bias_regression_calibration(b,"default","bootstrap","pmzsd","wild","mammen"),"absolute")
  }
})

test_that("fixed SD pivots match independent full residual-pilot matrices", {
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
tail_bands<-function(R,center,se,map,alpha,calibration) {
 B<-ncol(R);ans<-setNames(lapply(1:3,function(i)array(0,c(nrow(R),2L,length(alpha)))),c('pointwise','simultaneous','bonferroni'))
 for(target in unique(map$target)) {
  ii<-which(map$target==target);rr<-R[ii,,drop=FALSE];m<-length(ii)
  sr<-t(apply(rr,1L,sort));ar<-t(apply(abs(rr),1L,sort));lo<-sort(apply(rr,2L,min));hi<-sort(apply(rr,2L,max));ab<-sort(apply(abs(rr),2L,max))
  for(a in seq_along(alpha)) {
   k<-1+floor(B*alpha[a]/2);kb<-1+floor(B*alpha[a]/(2*m));ka<-ceiling(B*(1-alpha[a]));kab<-ceiling(B*(1-alpha[a]/m))
   for(band in names(ans)) {
    low<-switch(band,pointwise=sr[,k],simultaneous=rep(lo[k],m),bonferroni=pmin(sr[,kb],lo[k]))
    high<-switch(band,pointwise=sr[,B+1-k],simultaneous=rep(hi[B+1-k],m),bonferroni=pmax(sr[,B+1-kb],hi[B+1-k]))
    q<-switch(band,pointwise=ar[,ka],simultaneous=rep(ab[ka],m),bonferroni=pmax(ar[,kab],ab[ka]))
    if(calibration=='absolute'){low<- -q;high<-q}
    if(calibration=='tail-envelope'){low<-pmin(low,-q);high<-pmax(high,q)}
    ans[[band]][ii,1L,a]<-center[ii]-se[ii]*high
    ans[[band]][ii,2L,a]<-center[ii]-se[ii]*low
   }
  }
 };ans
}
# An independent full weighted-design SVD projection, never the production
# residual pilot or its native hat-matrix provider.
full_pilot<-function(st) {
 n<-nrow(st$x);H<-matrix(0,n,n);Q<-st$Q
 groups<-if(isTRUE(st$pooled))list(seq_len(n))else st$groups
 h<-if(isTRUE(st$pooled))st$h$p0 else st$b
 for(ii in groups)for(j in ii) {
  g<-rep(1,length(ii));for(k in seq_along(h))g<-g*exp(-.5*((st$xcon[[k]][ii]-st$xcon[[k]][j])/h[k])^2)
  Z<-Q[ii,,drop=FALSE]*sqrt(g);norm<-sqrt(colSums(Z^2));v<-svd(sweep(Z,2L,norm,'/'))
  stopifnot(all(v$d>0));H[j,ii]<-drop((Q[j,,drop=FALSE]/norm)%*%v$v%*%(t(v$u)/v$d))*sqrt(g)
 }
 H
}
full_smoother<-function(st) {
 n<-nrow(st$x);S<-matrix(0,n,n)
 for(ii in st$groups) {
  K<-matrix(1,length(ii),length(ii))
  for(k in seq_along(st$variance_h))K<-K*exp(-.5*(outer(st$xcon[[k]][ii],st$xcon[[k]][ii],'-')/st$variance_h[k])^2)
  S[ii,ii]<-K/rowSums(K)
 };S
}
  set.seed(112811L);held<-.Random.seed
  for(pooled in c(FALSE,TRUE))for(kernel in c("epanechnikov","beta")) {
    x<-if(pooled)data.frame(x=rep(seq(-.9,.9,length.out=20L),3L),group=factor(rep(letters[1:3],each=20L)))else data.frame(x=seq(-.9,.9,length.out=60L))
    y<-sin(3*x$x)+.15*cos(seq_len(60L))+if(pooled).2*as.integer(x$group)else 0
    args<-list(xdat=x,ydat=y,bws=c(if(kernel=="beta").3 else .5,if(pooled).2),regtype="ll",ckertype=kernel,bandwidth.compute=FALSE)
    if(kernel=="beta")args<-c(args,list(ckerbound="fixed",ckerlb=-1,ckerub=1))
    bw<-do.call(npregbw,args);grid<-x[c(3L,23L,43L,16L,36L,56L),,drop=FALSE]
    targets<-c("mean","d1",if(pooled)"c2");map<-expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE)
    st<-.np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law="mammen",calibration="absolute")
    H<-full_pilot(st);A<-diag(nrow(x))-H;q<-rowSums(A^2);res<-drop(A%*%y)/sqrt(q);mu<-drop(H%*%y);rmu<-drop(A%*%mu)/sqrt(q)
    S<-full_smoother(st);v<-drop(S%*%res^2)
    expect_equal(q,st$residual_row_variance,tolerance=1e-10);expect_equal(res,st$residual,tolerance=1e-10);expect_equal(v,st$variance,tolerance=1e-10)
    rows<-.np_bias_regression_rows(st,grid,map);active<-rows$se>0;B<-63L;alpha<-.2
    expect_equal(sqrt(drop(rows$W^2%*%v)),rows$se,tolerance=1e-10)
    if(pooled)expect_true(any(!active))
    for(law in c("gaussian","mammen"))for(cal in if(law=="gaussian")"absolute"else c("absolute","separate-tails","tail-envelope")) {
      actual.state<-.np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law=law,calibration=cal)
      same<-actual.state;same$settings$kernel.admission$calibration<-"absolute";expect_identical(same,st)
      z<-.np_with_seed(112812L,.np_bias_inference(actual.state,grid,targets,alpha,B,law,calibration=cal,resolution.warn=FALSE))
      draw<-.np_bias_stream(z$rng$start,{
        n<-nrow(x);xi<-if(law=="gaussian")matrix(rnorm(n*B),n)else{p<-(sqrt(5)+1)/(2*sqrt(5));matrix(ifelse(runif(n*B)<p,(1-sqrt(5))/2,(1+sqrt(5))/2),n)}
        ep<-sqrt(v)*xi;num<-rows$W%*%ep;den<-sqrt((rows$W^2%*%S)%*%(A%*%ep/sqrt(q)+rmu)^2)
        stopifnot(all(den[active,,drop=FALSE]>0),all(num[!active,,drop=FALSE]==0));den[!active,]<-1;num/den
      })
      expected<-tail_bands(draw$value,rows$center,rows$se,map,alpha,cal)
      for(band in names(expected))expect_equal(unname(z$bands[[band]]),unname(expected[[band]]),tolerance=1e-10)
      expect_identical(z$rng$end,draw$end);expect_identical(.Random.seed,held)
    }
  }
})
