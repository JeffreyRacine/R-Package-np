test_that("opaque custom replay refuses only participating uncached streams", {
  ns<-asNamespace("np");f<-get(".np_bias_replay_check",ns)
  set.seed(8001);before<-.Random.seed;count<-0L
  for(uniform in 0:7)for(normal in 0:5)for(law in c("gaussian","inid","fixed","geom","mammen","rademacher"))for(bl in c(1L,4L)) {
    seed<-as.integer(c(uniform+100L*normal+10000L,rep(123,7)))
    consumes.normal<-law=="gaussian" || (law=="geom" && bl>1L)
    unsafe<-uniform==5L || (consumes.normal && normal%in%c(2L,3L))
    if(unsafe)expect_error(f(law,bl,TRUE,seed=seed),"unsafe RNG replay",fixed=TRUE) else expect_null(f(law,bl,TRUE,seed=seed))
    expect_null(f(law,bl,FALSE,seed=seed));expect_null(f(law,bl,TRUE,retained=TRUE,seed=seed));count<-count+1L
  }
  expect_identical(count,576L);expect_identical(.Random.seed,before)
  # Previously diagnosed absent-uniform and Box-Muller errors retain precedence.
  expect_error(f("inid",1L,TRUE,seed=10405L),"custom uniform state absent",fixed=TRUE)
  expect_error(f("gaussian",1L,TRUE,seed=as.integer(c(10205,1:3))),"Box-Muller normal state",fixed=TRUE)
})
