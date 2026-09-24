# C179 support qualification; deleted-NN geometry is separately tracked.
# Private retained contexts must receive support during automatic-degree search.
test_that("automatic-degree searches retain declared training support", {
  skip_on_cran()
  old<-options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(.05,.95,length.out=12),
    o=ordered(rep(c(.5,1.5,3.5),4),levels=seq(.5,15.5)))
  y<-sin(4*x$x)+cos(seq_len(12))/9
  yd<-data.frame(y=ordered(rep(c(1.5,2.5,4.5,2.5),3),levels=seq(.5,15.5)))
  for(family in c("regression","condensity","condistribution")) {
    constructor<-get(switch(family,regression="npregbw",condensity="npcdensbw",
                            condistribution="npcdistbw"),asNamespace("np"))
    a<-if(family=="regression")list(xdat=x,ydat=y,bws=c(.4,.3),okertype="racineliyan") else
      list(xdat=x,ydat=yd,bws=c(.3,.4,.3),oxkertype="racineliyan",oykertype="racineliyan")
    if(family=="condensity")a$bwmethod<-"cv.ml"
    set.seed(921)
    b<-do.call(constructor,c(a,list(regtype="lp",nomad=TRUE,search.engine="nomad",
      degree.select="coordinate",degree.min=0L,degree.max=2L,degree.start=1L,
      degree.verify=FALSE,nmulti=1L,itmax=15L,powell.remin=FALSE,
      nomad.opts=list(MAX_BB_EVAL=20L))))
    evaluate<-get(switch(family,regression=".npregbw_eval_only",
      condensity=".npcdensbw_eval_only",condistribution=".npcdistbw_eval_only"),
      asNamespace("np"))
    result<-evaluate(x,if(family=="regression")y else yd,bws=b,invalid.penalty="dbmax")
    # Match the absolute objective-unit gate of the independent literal-fold
    # support oracles; expect_equal rescales small objectives relatively.
    expect_true(abs(as.numeric(b$fval)-as.numeric(result$objective))<2e-10,
                info=family)
    expect_true(all(b$degree %in% 0:2))
  }
})
