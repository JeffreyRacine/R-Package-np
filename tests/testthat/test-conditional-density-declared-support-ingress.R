# C179 support qualification; deleted-NN geometry is separately tracked.
conditional_declared_support_contract <- function(package) {
  ns<-asNamespace(package);bwfun<-get("npcdensbw",ns);fitfun<-get("npcdens",ns)
  evaluate<-get(".npcdensbw_eval_only",ns)
  old<-options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  n<-12L;x<-data.frame(x=c(.04,.12,.19,.27,.34,.41,.48,.53,.62,.71,.78,.91),
    u=factor(rep(c(1,2,4),4),levels=1:16),
    o=ordered(rep(c(.5,1.5,3.5),4),levels=seq(.5,15.5)))
  y<-data.frame(y=ordered(rep(c(1.5,2.5,4.5,2.5),3),levels=seq(.5,15.5)))
  for(type in c("fixed"))
    for(degree in c(0L,2L)) {
      b<-bwfun(xdat=x,ydat=y,bws=c(.3,if(type=="fixed").3 else 7,.25,.3),
        bwtype=type,bwmethod="cv.ml",regtype="lp",degree=degree,
        oxkertype="racineliyan",oykertype="racineliyan",bandwidth.compute=FALSE)
      fits<-vapply(seq_len(n),function(i)fitted(fitfun(bws=b,
        txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
        exdat=x[i,,drop=FALSE],eydat=y[i,,drop=FALSE])),0)
      # LP densities need not be positive. Independently apply the established
      # guarded CVML criterion, not log() of a potentially negative density.
      expected<-sum(vapply(fits,function(f) {
        if(f>0) log(f)
        else if(f< -.Machine$double.xmin)
          2*log(.Machine$double.xmin)-log(-f)
        else log(.Machine$double.xmin)
      },0))
      for(tree in c(FALSE,TRUE)) {
        options(np.tree=tree)
        observed<-evaluate(x,y,b,invalid.penalty="dbmax")$objective
        expect_true(abs(observed-expected)<2e-9,info=paste(type,degree,tree))
      }
    }
  # Independent finite-support CVLS contraction, no package weight helper.
  options(np.tree=FALSE)
  levels<-seq(.5,15.5);xo<-as.numeric(as.character(x$o))
  yo<-as.numeric(as.character(y$y));xu<-as.integer(x$u)
  w<-dnorm(outer(x$x,x$x,"-"),sd=.3)*
    ifelse(outer(xu,xu,"=="),.75,.25/15)*
    (.3^abs(outer(xo,xo,"-"))/rowSums(.3^abs(outer(xo,levels,"-"))))
  probabilities<-.3^abs(outer(yo,levels,"-"))
  probabilities<-probabilities/rowSums(probabilities)
  expected<-mean(vapply(seq_len(n),function(i) {
    weights<-w[,i];weights[i]<-0;weights<-weights/sum(weights)
    f<-drop(crossprod(weights,probabilities))
    2*f[match(yo[i],levels)]-sum(f^2)
  },0))
  b<-bwfun(xdat=x,ydat=y,bws=c(.3,.3,.25,.3),bwmethod="cv.ls",regtype="lc",
    oxkertype="racineliyan",oykertype="racineliyan",bandwidth.compute=FALSE)
  for(tree in c(FALSE,TRUE)) {
    options(np.tree=tree)
    observed<-evaluate(x,y,b,invalid.penalty="dbmax")$objective
    expect_true(abs(observed-expected)<2e-10,info=paste("cvls",tree,observed,expected))
  }
  options(np.tree=FALSE)
  for(solver in c("powell","mads")) {
    set.seed(21)
    b<-bwfun(xdat=x,ydat=y,bws=c(.3,.3,.25,.3),bwmethod="cv.ml",regtype="lc",
      oxkertype="racineliyan",oykertype="racineliyan",bwsolver=solver,
      nmulti=1L,itmax=10L,powell.remin=FALSE,nomad.opts=list(MAX_BB_EVAL=15L))
    observed<-evaluate(x,y,b,invalid.penalty="dbmax")$objective
    expect_equal(as.double(b$fval),as.double(observed),tolerance=2e-10)
  }
}
test_that("conditional density preparation retains Y and X declared support", {
  skip_on_cran()
  conditional_declared_support_contract("np")
})
