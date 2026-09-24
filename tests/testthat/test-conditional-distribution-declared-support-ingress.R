# C179 support qualification; deleted-NN geometry is separately tracked.
conditional_distribution_declared_support_contract <- function(package) {
  ns<-asNamespace(package);bwfun<-get("npcdistbw",ns);fitfun<-get("npcdist",ns)
  evaluate<-get(".npcdistbw_eval_only",ns)
  old<-options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-data.frame(x=c(.04,.12,.19,.27,.34,.41,.48,.53,.62,.71,.78,.91),
    u=factor(rep(c(1,2,4),4),levels=1:16),
    o=ordered(rep(c(.5,1.5,3.5),4),levels=seq(.5,15.5)))
  n<-nrow(x)
  y<-data.frame(y=ordered(rep(c(1.5,2.5,4.5,2.5),3),levels=seq(.5,15.5)))
  for(type in c("fixed")) for(degree in 0:2)
    for(empirical in c(FALSE,TRUE)) {
      if(empirical && type!="fixed") next
      grid<-if(empirical)y else data.frame(y=ordered(seq(.5,15.5),
                                                    levels=seq(.5,15.5)))
      b<-bwfun(xdat=x,ydat=y,bws=c(.4,if(type=="fixed").4 else 7,.65,.7),
        bwtype=type,regtype="lp",degree=degree,
        oxkertype="racineliyan",oykertype="racineliyan",bandwidth.compute=FALSE)
      expected<-mean(vapply(seq_len(n),function(i) {
        f<-fitted(fitfun(bws=b,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
                        exdat=x[rep(i,nrow(grid)),,drop=FALSE],eydat=grid))
        mean((as.integer(y$y[i]<=grid$y)-f)^2)
      },0))
      for(tree in c(FALSE,TRUE)) {
        options(np.tree=tree)
        observed<-evaluate(x,y,gydat=if(empirical)NULL else grid,
                           bws=b,invalid.penalty="dbmax")$objective
        expect_true(abs(observed-expected)<2e-10,
                    info=paste(type,degree,empirical,tree,observed,expected))
      }
    }
  options(np.tree=FALSE)
  for(solver in c("powell","mads")) {
    set.seed(21)
    b<-bwfun(xdat=x,ydat=y,bws=c(.4,.4,.65,.7),regtype="lc",
      oxkertype="racineliyan",oykertype="racineliyan",bwsolver=solver,
      nmulti=1L,itmax=10L,powell.remin=FALSE,nomad.opts=list(MAX_BB_EVAL=15L))
    observed<-evaluate(x,y,bws=b,invalid.penalty="dbmax")$objective
    expect_equal(as.double(b$fval),as.double(observed),tolerance=2e-10)
  }
}
test_that("conditional CDF preparation preserves declared Y and X support", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  conditional_distribution_declared_support_contract("npRmpi")
})
