# C179 support qualification; deleted-NN geometry is separately tracked.
declared_density_category <- function(train, evaluation, lambda, type, support) {
  gap <- abs(outer(train,evaluation,"-"))
  if(type=="aitchisonaitken")
    return(ifelse(gap==0,1-lambda,lambda/(length(support)-1)))
  if(type=="liracine")
    return(ifelse(gap==0,1,lambda)/(1+(length(support)-1)*lambda))
  lambda^gap/rowSums(lambda^abs(outer(train,support,"-")))
}

test_that("density criteria retain declared categorical support", {
  skip_on_cran()
  old<-options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-c(-1.31,-.89,-.41,-.17,.24,.52,1.06,1.57);n<-length(x)
  for(type in c("aitchisonaitken","liracine","racineliyan"))
    for(levels in list(0:3,0:15)) {
      category<-factor(rep(c(0,1,3,1),2),levels=levels,ordered=type=="racineliyan")
      dat<-data.frame(x=x,category=category);h<-.4;lambda<-.3
      codes<-if(is.ordered(category))as.numeric(as.character(category)) else as.integer(category)
      support<-if(is.ordered(category))as.numeric(levels(category)) else seq_along(levels(category))
      weights<-declared_density_category(codes,support,lambda,type,support)
      integral<-mean(dnorm(outer(x,x,"-"),sd=sqrt(2)*h)*tcrossprod(weights))
      cross<-dnorm(outer(x,x,"-"),sd=h)*declared_density_category(codes,codes,lambda,type,support)
      diag(cross)<-0
      expected<-2*sum(cross)/(n*(n-1))-integral
      args<-list(dat=dat,bws=c(h,lambda),bwmethod="cv.ls",bandwidth.compute=FALSE)
      if(is.ordered(category))args$okertype<-type else args$ukertype<-type
      b<-do.call(npudensbw,args)
      for(tree in c(FALSE,TRUE)) {
        options(np.tree=tree)
        got<-npudensbw(dat=dat,bws=b,eval.only=TRUE,nmulti=1,bwsolver="powell")
        expect_true(abs(got$fval-expected)<2e-10,info=paste(type,length(levels),tree))
      }
    }
})

test_that("native density support validation precedes preparation", {
  skip_on_cran()
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  dat<-data.frame(x=seq(.1,.9,length.out=8),u=factor(rep(1:2,4),levels=1:4))
  b<-npudensbw(dat=dat,bws=c(.2,.3),bandwidth.compute=FALSE)
  p<-getFromNamespace(".npudensbw_nomad_native_prepare_args","np")(dat,b,invalid.penalty="dbmax")
  native<-function(support) .Call("C_np_density_bw_eval",
    p$duno,p$dord,p$dcon,p$mysd,p$myopti,p$myoptd,as.double(b$bw),1L,
    p$penalty_mode,p$penalty_multiplier,p$ckerlb,p$ckerub,support,PACKAGE="np")
  expect_error(native(list(numeric())),"nonempty numeric levels",fixed=TRUE)
  expect_error(native(list(1:4)),"nonempty numeric levels",fixed=TRUE)
  expect_error(native(list(c(1,2,3,4),c(1,2))),"categorical dimensions",fixed=TRUE)
  expect_error(native(list()),"categorical dimensions",fixed=TRUE)
  expect_error(native(list(c(1,2,2,4))),"finite and increasing",fixed=TRUE)
  expect_error(native(list(c(1,3,4))),"outside declared categorical support",fixed=TRUE)
  expect_error(native(list(c(1,2,Inf))),"finite and increasing",fixed=TRUE)
  expect_true(is.finite(native(p$declared.support)$fval[1]))
  expect_identical(unname(getFromNamespace(".np_native_categorical_support","np")(b)),list(as.double(1:4)))
})

test_that("density CVML and native MADS retain declared support", {
  skip_on_cran()
  old<-options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  dat<-data.frame(x=c(-1.31,-.89,-.41,-.17,.24,.52,1.06,1.57),
                 category=factor(rep(c(0,1,3,1),2),levels=0:15))
  n<-nrow(dat);h<-.4;lambda<-.3
  for (type in c("aitchisonaitken","liracine","racineliyan")) {
    if(type=="racineliyan") dat$category<-ordered(dat$category)
    codes<-if(is.ordered(dat$category))as.numeric(as.character(dat$category)) else as.integer(dat$category)
    support<-if(is.ordered(dat$category))as.numeric(levels(dat$category)) else seq_along(levels(dat$category))
    category<-declared_density_category(codes,codes,lambda,type,support)
    for (bwtype in c("fixed")) {
      args<-list(dat=dat,bws=c(if(bwtype=="fixed")h else 4,lambda),
                 bwmethod="cv.ml",bwtype=bwtype,bandwidth.compute=FALSE)
      if(type=="racineliyan")args$okertype<-type else args$ukertype<-type
      b<-do.call(npudensbw,args)
      expected<-sum(vapply(seq_len(n),function(i) {
        donors<-seq_len(n)[-i]
        radii<-if(bwtype=="fixed")h else if(bwtype=="generalized_nn")
          sort(abs(dat$x[donors]-dat$x[i]))[4] else
          vapply(donors,function(j)sort(abs(dat$x[donors]-dat$x[j]))[5],0)
        log(mean(dnorm(dat$x[donors]-dat$x[i],sd=radii)*category[donors,i]))
      },0))
      got<-npudensbw(dat=dat,bws=b,eval.only=TRUE,nmulti=1,
                     invalid.penalty="dbmax")$fval
      expect_true(abs(got-expected)<2e-9,info=paste(type,bwtype))
    }
    # Both search owners must hand the same support to the final objective.
    args$bwtype<-"fixed";args$bws<-c(h,lambda);args$bandwidth.compute<-TRUE
    for(solver in c("powell","mads")) {
      set.seed(19)
      fit<-do.call(npudensbw,c(args,list(bwsolver=solver,nmulti=1L,
        itmax=20L,powell.remin=FALSE,nomad.opts=list(MAX_BB_EVAL=20L))))
      check<-npudensbw(dat=dat,bws=fit,eval.only=TRUE,nmulti=1L,
                       invalid.penalty="dbmax")
      expect_equal(as.double(fit$fval),as.double(check$fval),tolerance=2e-10)
    }
  }
})
