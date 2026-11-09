test_that("asymmetric simultaneous roots invert independently ranked family tails", {
  B<-99L;m<-65L;j<-seq_len(B)
  roots<-outer(seq(0,3,length.out=m),seq(-4,2,length.out=B),"+")
  roots[33L,]<-sin(j/9)+seq(-2,3,length.out=B)
  rhat<-seq(-.2,.3,length.out=m);scale<-seq(.5,2,length.out=m)
  region<-list(lower=rep(-.2,m),upper=rep(.8,m));center<-rhat+.1
  store<-.np_udens_beta_store();on.exit(.np_udens_beta_store_close(store),add=TRUE)
  root<-.np_udens_beta_root_store(store,m,B)
  root<-.np_udens_beta_root_write(root,roots,1L)
  for(alpha in c(.1,.2)){
    k<-floor(B*alpha/4)+1L
    high<-sort(vapply(seq_len(B),function(i)max(roots[,i]),0))
    low<-sort(vapply(seq_len(B),function(i)min(roots[,i]),0))
    expected<-list(lower=rhat+region$lower-scale*high[B+1L-k],upper=rhat+region$upper-scale*low[k])
    stored<-.np_udens_beta_stored_reduce(root,rhat,scale,region,center,alpha,"simultaneous","asymmetric-sd")
    direct<-.np_udens_beta_sd_reduce(roots,rhat,scale,region,center,alpha,"simultaneous","asymmetric-sd")
    expect_identical(stored$k_signed,k)
    expect_equal(stored$lower,expected$lower,tolerance=0)
    expect_equal(stored$upper,expected$upper,tolerance=0)
    expect_equal(direct$lower,expected$lower,tolerance=0)
    expect_equal(direct$upper,expected$upper,tolerance=0)
  }
})
