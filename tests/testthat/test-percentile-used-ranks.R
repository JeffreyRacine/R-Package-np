test_that("calibrated diagnostics report the weakest used empirical rank", {
  d <- .np_band_resolution(199L, rep("mean", 3), .05, rep(TRUE, 3),
    c("pointwise", "simultaneous", "bonferroni"), signed=TRUE,
    simultaneous.index=matrix(1L))
  before <- d
  z <- .np_band_resolution_use_indices(d, list(c(3, 2, 4), 1, 2))
  expect_identical(z$order_index, c(2, 1, 2))
  expect_identical(z$resolution_limited, c(FALSE, TRUE, FALSE))
  expect_identical(z$nominal_tail_draws, before$nominal_tail_draws)
  expect_identical(d, before)
  expect_warning(.np_band_resolution_warn(z), "simultaneous",
    class="np_band_resolution_warning")
  expect_error(.np_band_resolution_use_indices(d, list(0, 1, 2)), "outside")
  expect_error(.np_band_resolution_use_indices(d, list(101, 1, 2)), "outside")
  expect_error(.np_band_resolution_use_indices(d, list(NA_real_, 1, 2)), "invalid")
})

test_that("calibrated rank metadata matches independent finite-bank endpoints", {
  en <- new.env(parent=environment(.np_reg_percentile_infer))
  en$.np_reg_percentile_infer <- .np_reg_percentile_infer
  environment(en$.np_reg_percentile_infer) <- en
rows<-list();raw<-list()
for(B in c(3L,19L,199L,200L))for(m in c(1L,3L))for(boundary in c(FALSE,TRUE)){
 set.seed(617);curves<-matrix(round(rnorm(2*m*B),1),2*m,B)
 support<-cbind(apply(curves,1,min)-1,apply(curves,1,max)+1)
 mean.score<-if(boundary)B+1 else ceiling(.83*B)
 marg<-if(boundary)rep(c(B,B+1),length.out=m)else floor(seq(0,.6*B,length.out=m))
 fam<-if(boundary)B+1 else ceiling(.4*B)
 scores<-c(mean.score,marg,fam)
 en$.np_reg_percentile_refresh_geometry<-function(p)NULL
 en$.np_reg_percentile_outer_response_bank<-function(...)list(y=rep(list(rep(0,4)),39),raw=matrix(0,4,39),v=matrix(1,4,39))
 en$.np_reg_percentile_curve_draws<-function(...)list(curves=curves)
 en$.np_reg_percentile_curve_support<-function(...)support
 en$.np_reg_percentile_outer_dispatch<-function(...)matrix(rep(c(scores,scores),each=39),39)
 p<-list(A=matrix(0,2*m,4),y=rep(0,4),mu=rep(0,4),sigma=rep(1,4),center=rep(0,2*m),se=rep(1,2*m),
  map=expand.grid(row=seq_len(m),target=c('mean','d1'),stringsAsFactors=FALSE))
 z<-en$.np_reg_percentile_infer(p,B,10L,1:39,101:139)
 used<-list();expected<-list()
 for(target in c('mean','d1')){
  ii<-which(p$map$target==target);a<-curves[ii,,drop=FALSE];sorted<-t(apply(a,1,sort))
  K<-apply(a,2,function(v)min(vapply(seq_len(m),function(i)min(sum(a[i,]<=v[i]),sum(a[i,]>=v[i])),0)))
  nominal.p<-1+floor(B*.05/2);nominal.b<-1+floor(B*.05/(2*m))
  if(target=='mean'){
   ks<-if(mean.score>B)1 else min(nominal.p,1+floor((B-mean.score)/2),sort(K)[B-mean.score+1])
   kp<-rep(nominal.p,m);ks<-rep(ks,m)
  }else{
   kp<-ifelse(marg>B,1,ifelse(marg==0,ceiling(B/2),1+floor((B-marg)/2)))
   ks<-if(fam>B)1 else if(fam==0)ceiling(B/2)else min(1+floor((B-fam)/2),sort(K)[B-fam+1])
   ks<-pmin(kp,ks)
  }
  kb<-pmin(nominal.b,ks)
  if(m==1)kp<-kb<-ks
  used[[target]]<-list(pointwise=kp,simultaneous=ks,bonferroni=kb)
  for(band in names(used[[target]])){
   kk<-used[[target]][[band]]
   rr<-z$resolution$family==target & z$resolution$band==band
   pass<-identical(as.numeric(z$resolution$order_index[rr]),as.numeric(min(kk))) &&
     identical(z$resolution$resolution_limited[rr],min(kk)==1)
   rows[[length(rows)+1L]]<-data.frame(B,m,boundary,target,band,reported=z$resolution$order_index[rr],expected=min(kk),pass)
   # Independently recover the empirical component of every endpoint.
   got<-matrix(z$bands[[band]][ii,,,drop=FALSE],m,2)
   lo<-sorted[cbind(seq_len(m),kk)];hi<-sorted[cbind(seq_len(m),B+1-kk)]
   extend<-rep(FALSE,m)
   if(target=='mean'&&mean.score>=B&&(band!='pointwise'||m==1))extend[]<-TRUE
   if(target=='d1'){
    extend<-marg>=B
    if((band!='pointwise'||m==1)&&fam>=B)extend[]<-TRUE
   }
   lo[extend]<-pmin(lo[extend],support[ii,1][extend]);hi[extend]<-pmax(hi[extend],support[ii,2][extend])
   stopifnot(identical(as.numeric(got[,1]),as.numeric(lo)),identical(as.numeric(got[,2]),as.numeric(hi)))
  }
  if(target=='mean'){
   band<-matrix(z$bands$simultaneous[ii,,,drop=FALSE],m,2)
   count<-sum(apply(a>=band[,1]&a<=band[,2],2,all))
   rows[[length(rows)+1L]]<-data.frame(B,m,boundary,target,band='content',reported=z$ranks$mean$content,expected=count,
     pass=identical(as.numeric(z$ranks$mean$content),as.numeric(count)))
  }
 }
 raw[[paste(B,m,boundary)]]<-list(curves=curves,used=used,result=z)
}
  tab <- do.call(rbind, rows)
  expect_equal(nrow(tab), 112L)
  expect_true(all(tab$pass))
})
