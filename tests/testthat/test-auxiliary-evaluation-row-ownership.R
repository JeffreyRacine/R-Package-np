test_that("quantile and class formula outputs use their evaluation row owner", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d <- data.frame(x=seq(-1,1,length.out=20),y=sin(seq_len(20)),
                  cls=factor(rep(letters[1:2],10)))
  d$x[3] <- NA_real_
  grid <- data.frame(x=c(.6,-.7,.2,-.1))
  for(family in c("npqreg","npconmode")) {
    b <- if(family=="npqreg")
      npcdistbw(y~x,data=d,bws=c(.4,.3),bandwidth.compute=FALSE,na.action=na.exclude) else
      npcdensbw(cls~x,data=d,bws=c(.2,.3),bandwidth.compute=FALSE,na.action=na.exclude)
    ff <- get(family)
    extra <- if(family=="npqreg") list(tau=c(.25,.75),gradients=TRUE,se=TRUE) else
      list(probabilities=TRUE,gradients=TRUE,se=TRUE)
    training <- do.call(ff,c(list(bws=b),extra))
    expect_equal(NROW(fitted(training)),nrow(d))
    expect_true(all(is.na(as.matrix(fitted(training))[3,,drop=FALSE])))
    for(missing.row in c(FALSE,TRUE)) {
      e <- grid
      if(missing.row) e$x[2] <- NA_real_
      a <- do.call(ff,c(list(bws=b,exdat=e),extra))
      z <- do.call(ff,c(list(bws=b,newdata=e),extra))
      both <- do.call(ff,c(list(bws=b,exdat=e,newdata=grid[1:2,,drop=FALSE]),extra))
      expect_equal(fitted(a),fitted(z),tolerance=2e-11)
      expect_equal(fitted(a),fitted(both),tolerance=0)
      expect_equal(NROW(fitted(a)),4L)
      if(family=="npqreg") {
        expect_equal(se(a),se(z),tolerance=2e-11)
        expect_equal(gradients(a),gradients(z),tolerance=2e-11)
      } else {
        expect_equal(a$probabilities,z$probabilities,tolerance=2e-11)
      }
    }
  }
})
