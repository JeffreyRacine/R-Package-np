test_that("single-index reentry preserves original response and row metadata", {
  withr::local_options(np.messages=FALSE, np.plot.progress=FALSE)
  for (method in c("ichimura", "kleinspady"))
    for (type in c("fixed", "generalized_nn", "adaptive_nn"))
      for (named in c(FALSE,TRUE)) {
        set.seed(731); n<-40L
        d<-data.frame(x=rnorm(n),z=rnorm(n))
        d$y<-if(method=="kleinspady")rbinom(n,1,plogis(d$x+.5*d$z))else rpois(n,2)
        if(named)row.names(d)<-paste0("subject-",seq_len(n))
        formula<-stats::as.formula("y~x+z",env=baseenv())
        b<-npindexbw(formula,data=d,bws=c(1,.5,if(type=="fixed").8 else 20),
          bwtype=type,method=method,bandwidth.compute=FALSE)
        original<-b
        f<-npindex(b,se=FALSE)
        expect_identical(f$bws$.np.formula.training,b$.np.formula.training)
        expect_identical(b,original)
        expect_type(f$bws$.np.formula.training$frame$y,"integer")
        expect_identical(fitted(npindex(unserialize(serialize(b,NULL)),se=FALSE)),fitted(f))
        # Explicit replacement/reordering must still retain the actual new sample.
        replacement<-d[n:1,]; replacement$y<-rev(d$y)
        g<-npindex(b,data=replacement,se=FALSE)
        expect_identical(g$bws$.np.formula.training$frame$y,replacement$y)
        expect_identical(row.names(g$bws$.np.formula.training$frame),row.names(replacement))
        changed<-d; changed$y<-if(method=="kleinspady")1L-d$y else d$y+2L
        g<-npindex(b,data=changed,se=FALSE)
        expect_identical(g$bws$.np.formula.training$frame$y,changed$y)
        expect_false(identical(fitted(f),fitted(g)))
      }
})

test_that("single-index retained omissions survive conversion and explicit refit", {
  withr::local_options(np.messages=FALSE,np.plot.progress=FALSE)
  set.seed(52);d<-data.frame(y=rpois(45,2),x=rnorm(45),z=rnorm(45));d$x[3]<-NA
  formula<-stats::as.formula("y~x+z",env=baseenv())
        b<-npindexbw(formula,data=d,na.action=na.exclude,bws=c(1,.5,.8),bandwidth.compute=FALSE)
  f<-npindex(b,se=TRUE)
  expect_identical(f$bws$.np.formula.training,b$.np.formula.training)
  expect_length(fitted(f),45L);expect_true(is.na(fitted(f)[3L]))
  expect_true(all(is.finite(fitted(f)[-3L])))
  expect_identical(predict(f,newdata=d[4:9,c("x","z")]),predict(f,exdat=d[4:9,c("x","z")]))
})
