# C179 support qualification; deleted-NN geometry is separately tracked.
regression_declared_support_contract <- function(package) {
  ns<-asNamespace(package);bwfun<-get("npregbw",ns);fitfun<-get("npreg",ns)
  evaluate<-get(".npregbw_eval_only",ns);hatfun<-get("npreghat",ns)
  old<-options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-data.frame(x=c(.04,.12,.19,.27,.34,.41,.48,.53,.62,.71,.78,.91),
    u=factor(rep(c(1,2,4),4),levels=1:16),
    o=ordered(rep(c(.5,1.5,3.5),4),levels=seq(.5,15.5)))
  n<-nrow(x);y<-sin(4*x$x)+as.integer(x$u)/8+cos(seq_len(n))/9
  for(type in c("fixed"))
    for(degree in 0:2) {
      options(np.tree=FALSE)
      b<-bwfun(xdat=x,ydat=y,bws=c(if(type=="fixed").4 else 7,.65,.7),
        bwtype=type,regtype="lp",degree=degree,okertype="racineliyan",
        bandwidth.compute=FALSE)
      deleted<-vapply(seq_len(n),function(i)fitted(fitfun(bws=b,
        txdat=x[-i,,drop=FALSE],tydat=y[-i],exdat=x[i,,drop=FALSE])),0)
      for(tree in c(FALSE,TRUE)) {
        options(np.tree=tree)
        observed<-evaluate(x,y,b,invalid.penalty="dbmax")$objective
        expect_true(abs(observed-mean((y-deleted)^2))<2e-10,
                    info=paste(type,degree,tree))
      }
      # delta=.5 has zero normal-score shift: literal LOO check loss is
      # independent of the LSQ objective's shared preparation/optimizer.
      lsq<-get(".nplsqreg_call_fixed_degree_core",ns)
      loss<-function()lsq(x,y,rep(1,n),.5,b,.5,c(.01,.99),
                         list(invalid.penalty="dbmax"),FALSE)$objective
      observed<-if(package=="npRmpi")
        get(".npRmpi_with_local_regression",ns)(loss()) else loss()
      expect_true(abs(observed-mean(abs(y-deleted)/2))<2e-10,
                  info=paste("check",type,degree))
    }
  options(np.tree=FALSE)
  for(degree in 0:2) {
    b<-bwfun(xdat=x,ydat=y,bws=c(.4,.65,.7),bwmethod="cv.aic",
      regtype="lp",degree=degree,okertype="racineliyan",bandwidth.compute=FALSE)
    f<-fitted(fitfun(bws=b,txdat=x,tydat=y))
    tr<-sum(diag(hatfun(bws=b,txdat=x)))
    expected<-log(mean((y-f)^2))+(1+tr/n)/(1-(tr+2)/n)
    observed<-evaluate(x,y,b,invalid.penalty="dbmax")$objective
    expect_true(abs(observed-expected)<2e-10,info=paste("aic",degree))
  }
  for(solver in c("powell","mads")) {
    set.seed(21)
    b<-bwfun(xdat=x,ydat=y,bws=c(.4,.65,.7),regtype="lc",
      okertype="racineliyan",bwsolver=solver,nmulti=1L,itmax=10L,
      powell.remin=FALSE,nomad.opts=list(MAX_BB_EVAL=15L))
    observed<-evaluate(x,y,b,invalid.penalty="dbmax")$objective
    expect_equal(as.double(b$fval),as.double(observed),tolerance=2e-10)
  }
}
test_that("regression and LSQ preparation preserve declared categories", {
  skip_on_cran()
  regression_declared_support_contract("np")
})
