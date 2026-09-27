test_that("conditional GNN fitting selection has no degree or width cutoff", {
 old<-options(np.messages=FALSE,np.tree=TRUE);on.exit(options(old))
 ns<-asNamespace("np")
 eligible<-get(".np_conditional_gnn_lp_fit_tree_eligible",ns)
 density<-get(".npcdensbw_tree_code",ns)
 distribution<-get(".npcdistbw_tree_code",ns)
 yes<-get("DO_TREE_YES",ns);no<-get("DO_TREE_NO",ns)
 set.seed(759);x<-data.frame(x=runif(32));y<-rnorm(32)
 b<-npcdensbw(xdat=x,ydat=y,bws=c(24,24),bwtype="generalized_nn",
  regtype="lp",degree=2L,bernstein.basis=TRUE,bwmethod="cv.ls",
  cxkertype="epanechnikov",cykertype="epanechnikov",bandwidth.compute=FALSE)
 for(d in c(0L,1L,2L,3L,79L,80L,81L)) {
  q<-b;q$degree.engine<-d
  expect_true(eligible(q))
  expect_identical(density(q,2L,0L,fit.context=TRUE),yes)
  expect_identical(distribution(q,2L,0L,fit.context=TRUE),yes)
  # Whole-support density search belongs to C167's separate qualification.
 }
 for(mode in list(FALSE,"auto")) {
  options(np.tree=mode);expect_false(eligible(b))
 }
 options(np.tree=TRUE)
 for(type in c("fixed","adaptive_nn")) {q<-b;q$type<-type;expect_false(eligible(q))}
 for(kernel in c("gaussian","beta")) {q<-b;q$cxkertype<-kernel;expect_false(eligible(q))}
 q<-b;q$cxkerbound<-"fixed";expect_false(eligible(q))
 q<-b;q$regtype.engine<-"lc";expect_false(eligible(q))
})

test_that("conditional GNN X-tree fitting retains LP0 and donor/query identities", {
 old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old))
 set.seed(758);n<-128L;x<-data.frame(x=runif(n,-1,1));y<-sin(x$x)+rnorm(n,sd=.3)
 permutation<-sample.int(n)
 for(cdf in c(FALSE,TRUE))for(degree in c(0L,1L,3L)) {
  ctor<-if(cdf)npcdistbw else npcdensbw;fit<-if(cdf)npcdist else npcdens
  b<-ctor(xdat=x,ydat=y,bws=c(80,80),bwtype="generalized_nn",
   regtype="lp",degree=degree,bernstein.basis=TRUE,bwmethod="cv.ls",
   cxkertype="epanechnikov",cykertype="epanechnikov",bandwidth.compute=FALSE)
  for(external in c(FALSE,TRUE)) {
   ref<-NULL
   for(perm in list(seq_len(n),permutation))for(tree in c(FALSE,TRUE)) {
    args<-list(bws=b,txdat=x[perm,,drop=FALSE],tydat=y[perm],se=TRUE)
    if(external){args$exdat<-x[c(7,7),,drop=FALSE];args$eydat<-y[c(7,7)]}
    options(np.tree=tree);g<-do.call(fit,args)
    z<-cbind(fitted(g),se(g));if(!external)z<-z[order(perm),,drop=FALSE]
    if(is.null(ref))ref<-z
    expect_identical(is.na(z),is.na(ref))
    expect_true(max(abs(z-ref),na.rm=TRUE)<2e-9)
   }
  }
 }
})
