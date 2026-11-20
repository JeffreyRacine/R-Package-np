test_that("fixed non-Gaussian tails preserve calibration and pairs-law boundaries", {
  for(kernel in c("epanechnikov","beta"))for(cal in c("separate-tails","tail-envelope")) {
    b<-list(type="fixed",ckertype=kernel,ckerorder=2L,ckerbound=if(kernel=="beta")"fixed"else"none",ncon=1L)
    expect_identical(.np_bias_regression_kernel_class(b,cal),"local-envelope")
    admission<-list(type="fixed",calibration=cal)
    for(law in c("inid","fixed","geom"))expect_silent(.np_bias_regression_kernel_law(admission,law))
    for(law in c("gaussian","mammen","rademacher"))
      expect_error(.np_bias_regression_kernel_law(admission,law),"requires inid, fixed or geom")
    for(type in c("generalized_nn","adaptive_nn")) {
      z<-b;z$type<-type;expect_identical(.np_bias_regression_kernel_class(z,cal),"unsupported")
    }
    expect_identical(.np_bias_regression_calibration(b,"default","bootstrap","simultaneous","inid","mammen"),"basic-local-envelope")
  }
})

test_that("fixed tails match literal pooled replicate scales and family inversion", {
  skip_on_cran()
  if(!spawn_mpi_slaves(n=1L)) skip("MPI workers unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
variance<-function(F,L,law){F<-F-rowMeans(F);n<-ncol(F)
 if(L==1L)return(rowSums(F^2))
 if(law=='fixed'){
  part<-function(w)if(w==0)rep(0,nrow(F))else{
   Z<-sapply(seq_len(n),function(i)rowSums(F[,((i-1L+seq_len(w)-1L)%%n)+1L,drop=FALSE]));rowMeans(Z^2)}
  return((n%/%L)*part(L)+part(n%%L))}
 v<-rowSums(F^2)
 for(k in seq_len(n-1L))v<-v+2*(n-k)/n*(1-1/L)^k*rowSums(F*F[,((seq_len(n)-1L+k)%%n)+1L,drop=FALSE])
 v}
# Full-array order statistics/inversion. Shared only: qualified observation rows
# and primitive count generator. No inference reducer/studentizer is called here.
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
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  f<-function(n)get(n,environment(.np_bias_inference))
  x<-data.frame(x=rep(seq(-.9,.9,length.out=20L),3L),group=factor(rep(letters[1:3],each=20L)))
  y<-sin(3*x$x)+.2*as.integer(x$group)+.15*cos(seq_len(nrow(x)))
  grid<-x[c(3L,23L,43L,16L,36L,56L),,drop=FALSE];targets<-c("mean","d1","c2")
  map<-expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE);B<-63L;alpha<-.2
  set.seed(112421L);held<-.Random.seed
  for(kernel in c("epanechnikov","beta")) {
    args<-list(xdat=x,ydat=y,bws=c(if(kernel=="beta").3 else .5,.2),regtype="ll",ckertype=kernel,bandwidth.compute=FALSE)
    if(kernel=="beta")args<-c(args,list(ckerbound="fixed",ckerlb=-1,ckerub=1))
    bw<-do.call(npregbw,args)
    for(law in c("inid","fixed","geom"))for(L in if(law=="inid")1L else c(1L,5L)) {
      local<-.np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law=law,calibration="basic-local-envelope")
      rows<-.np_bias_regression_rows(local,grid,map);F<-sweep(rows$W,2L,local$residual,"*");F<-F-rowMeans(F)
      expect_true(isTRUE(local$pooled));expect_true(any(rows$se==0))
      for(cal in c("separate-tails","tail-envelope")) {
        st<-.np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law=law,calibration=cal)
        same<-st;same$settings$kernel.admission$calibration<-"basic-local-envelope";expect_identical(same,local)
        z<-.np_with_seed(112422L,.np_bias_inference(st,grid,targets,alpha,B,law,if(law=="inid")NULL else L,calibration=cal,refinement=if(law=="inid")"none"else"pilot-block",resolution.warn=FALSE))
        expected.se<-rows$se
        if(law=="inid") {
          draw<-.np_bias_stream(z$rng$start,rmultinom(B,nrow(x),rep(1/nrow(x),nrow(x)))-1)
          raw<-F%*%draw$value;scale<-apply(raw,1L,sd);R<-raw/ifelse(scale==0,1,scale)
        } else {
          expected.se<-sqrt(variance(F,L,law));active<-expected.se>0
          draw<-.np_bias_stream(z$rng$start,{
            if(L==1L) {
              counts<-rmultinom(B,nrow(x),rep(1/nrow(x),nrow(x)));num<-F%*%counts
              V<-(F^2)%*%counts-num^2/nrow(x);stopifnot(all(V[active,,drop=FALSE]>0),all(num[!active,,drop=FALSE]==0))
              R0<-matrix(0,nrow(F),B);R0[active,]<-num[active,,drop=FALSE]/sqrt(V[active,,drop=FALSE]);R0
            } else {
              one<-.np_block_replicate_drawer(nrow(x),nrow(x),L,law,TRUE)
              sapply(seq_len(B),function(i){Z<-F[,one(),drop=FALSE];num<-rowSums(Z);V<-variance(Z,L,law);stopifnot(all(V[active]>0),all(num[!active]==0));v<-numeric(nrow(F));v[active]<-num[active]/sqrt(V[active]);v})
            }
          });R<-draw$value
        }
        expected<-tail_bands(R,rows$center,expected.se,map,alpha,cal)
        expect_equal(z$se,expected.se,tolerance=1e-10);expect_identical(z$rng$end,draw$end)
        for(band in names(expected))expect_equal(unname(z$bands[[band]]),unname(expected[[band]]),tolerance=1e-10)
        expect_identical(.Random.seed,held)
      }
    }
  }
})
