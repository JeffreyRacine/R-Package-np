test_that("fixed categorical plots subtract projected endpoint probabilities", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=60),c=factor(rep(c("a","b"),30)))
  y <- data.frame(y=factor(ifelse(x$x+.7*(x$c=="b")>0,"yes","no")))
  for (proper in c(FALSE,TRUE)) {
    fit <- npconmode(txdat=x,tydat=y,bws=c(0,.25,.03),regtype="ll",
      probabilities=TRUE,gradients=TRUE,level="yes",proper=proper)
    for (q in c(.05,.5,.95)) {
      payload <- plot(fit,gradients=TRUE,view="fixed",neval=5,
                      xq=c(q,.5),output="data",level="yes")
      ex <- data.frame(x=rep(as.double(quantile(x$x,q)),2),
                        c=factor(c("a","b"),levels=c("a","b")))
      direct <- npconmode(fit$bws,txdat=x,tydat=y,exdat=ex,
        probabilities=TRUE,gradients=TRUE,level="yes",proper=proper)
      probability <- direct$probabilities[,"yes"]
      expected <- c(0,probability[2]-probability[1])
      # Absolute oracle tolerance: near-zero raw contrasts can amplify a
      # relative comparison despite only 1e-15 absolute summation differences.
      expect_lte(max(abs(payload$c$effect-expected)),1e-12)
      expect_lte(max(abs(as.vector(gradients(direct)[,"c"])-expected)),1e-12)
    }
  }
})
