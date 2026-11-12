test_that("uniform ignored order does not certify empty-window uncertainty", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  x<-c(seq(.01,.19,length.out=25),seq(.81,.99,length.out=25))
  for(q in 1:2) {
    d<-if(q==1)data.frame(x=x) else data.frame(x=x,z=rev(x))
    f2<-npudens(tdat=d,bws=rep(.08,q),ckertype='uniform',ckerorder=2L)
    f4<-suppressWarnings(npudens(tdat=d,bws=rep(.08,q),ckertype='uniform',ckerorder=4L))
    expect_identical(fitted(f4),fitted(f2))
    for(band in c('pointwise','bonferroni')) {
      p2<-suppressWarnings(plot(f2,errors='asymptotic',band=band,neval=9,view='fixed',output='data'))
      p4<-suppressWarnings(plot(f4,errors='asymptotic',band=band,neval=9,view='fixed',output='data'))
      a<-as.data.frame(p2,panel=1);b<-as.data.frame(p4,panel=1)
      grid<-b[paste0('eval_',seq_len(q))]
      empty<-vapply(seq_len(nrow(grid)),function(i){
        positive<-rep(TRUE,nrow(d))
        for(j in seq_len(q))positive<-positive & abs(d[[j]]-grid[[j]][i]) < .08
        !any(positive)
      },TRUE)
      expect_true(any(empty))
      expect_true(all(b$estimate[empty]==0))
      expect_true(all(is.na(b$std_error[empty])))
      expect_true(all(b$interval_status[empty]=='undefined'))
      expect_identical(b[c('estimate','std_error','lower','upper','interval_status')],
                       a[c('estimate','std_error','lower','upper','interval_status')])
    }
  }
})
