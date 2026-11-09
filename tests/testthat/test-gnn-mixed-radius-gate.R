test_that("GNN shortcut requires every local radius to be large", {
 old<-options(np.messages=FALSE,np.tree=FALSE,np.largeh=TRUE,np.extendednn=TRUE,np.macMseries.accelerate=FALSE);on.exit(options(old),add=TRUE)
 x<-data.frame(x=c(seq(-1,1,length.out=44L),40:43))
 y<-sin(3*pmin(x$x,3))+cos(7*pmin(x$x,3))/4
 force({
  for(k in c(6L,100L,1000L)){
   b<-npregbw(xdat=x,ydat=y,bws=k,bandwidth.compute=FALSE,bwtype="generalized_nn",regtype="lp",degree=2L,bernstein.basis=TRUE)
   options(np.largeh=TRUE);selected<-.npregbw_eval_only(x,y,b,invalid.penalty="dbmax")
   options(np.largeh=FALSE);literal<-.npregbw_eval_only(x,y,b,invalid.penalty="dbmax")
   expect_identical(as.integer(selected$num.feval.fast),0L)
   expect_equal(selected$objective,literal$objective,tolerance=0)
  }
 })
})
