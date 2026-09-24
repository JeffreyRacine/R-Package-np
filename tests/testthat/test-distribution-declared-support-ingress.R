# C179 support qualification; deleted-NN geometry is separately tracked.
distribution_declared_support_contract <- function(package) {
  ns<-asNamespace(package);bwfun<-get("npudistbw",ns)
  old<-options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-c(-1.31,-.89,-.41,-.17,.24,.52,1.06,1.57);n<-length(x)
  support<-seq(.5,15.5);level<-rep(c(1.5,2.5,4.5,2.5),2)
  dat<-data.frame(x=x,o=ordered(level,levels=support))
  grid<-data.frame(x=c(-.73,.08,.81,1.43),
                   o=ordered(c(.5,3.5,6.5,15.5),levels=support))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(empirical in c(TRUE,FALSE)) {
      evaluation<-if(empirical)dat else grid
      evalue<-as.numeric(as.character(evaluation$o));m<-nrow(evaluation)
      lambda<-.3;h<-.4;k<-4L;expected<-0
      for(i in seq_len(n)) {
        donors<-seq_len(n)[-i]
        for(j in seq_len(m)) {
          if(empirical&&j==i)next
          radius<-if(type=="fixed")h else if(type=="generalized_nn")
            sort(abs(x[donors]-evaluation$x[j]))[k+as.integer(empirical)] else
            vapply(donors,function(d)sort(abs(x[donors]-x[d]))[k+1L],0)
          category<-vapply(donors,function(d)
            sum(lambda^abs(level[d]-support[support<=evalue[j]]))/
              sum(lambda^abs(level[d]-support)),0)
          fitted<-mean(pnorm((evaluation$x[j]-x[donors])/radius)*category)
          target<-as.integer(x[i]<=evaluation$x[j]&&level[i]<=evalue[j])
          expected<-expected+(fitted-target)^2
        }
      }
      expected<-expected/(n*(m-as.integer(empirical)))
      b<-bwfun(dat=dat,bws=c(if(type=="fixed")h else k,lambda),
        bwtype=type,okertype="racineliyan",bandwidth.compute=FALSE)
      for(tree in c(FALSE,TRUE)) {
        options(np.tree=tree)
        args<-list(dat=dat,bws=b,eval.only=TRUE,nmulti=1L,
                   do.full.integral=TRUE,invalid.penalty="dbmax")
        if(!empirical)args$gdat<-grid
        z<-do.call(bwfun,args)
        expect_true(abs(z$fval-expected)<2e-10,
                    info=paste(type,empirical,tree,z$fval,expected))
      }
    }
  options(np.tree=FALSE)
  for(solver in c("powell","mads")) {
    set.seed(81)
    b<-bwfun(dat=dat,bws=c(.4,.3),okertype="racineliyan",
      bwsolver=solver,nmulti=1L,itmax=20L,powell.remin=FALSE,
      nomad.opts=list(MAX_BB_EVAL=20L),do.full.integral=TRUE)
    z<-bwfun(dat=dat,bws=b,eval.only=TRUE,nmulti=1L,do.full.integral=TRUE,
             invalid.penalty="dbmax")
    expect_equal(as.double(b$fval),as.double(z$fval),tolerance=2e-10)
  }
}
test_that("distribution criteria retain declared ordered support", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  distribution_declared_support_contract("npRmpi")
})
