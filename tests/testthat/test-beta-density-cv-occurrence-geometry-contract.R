test_that("beta density CVML uses occurrence-excluded nearest neighbours", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(.02,.06,.13,.23,.38,.51,.67,.8,.89,.98))
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    bw <- npudensbw(dat=x,bws=if(type=="fixed") .3 else 3,
      bwtype=type,ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1,
      bandwidth.compute=FALSE,bwmethod="cv.ml")
    oracle <- vapply(seq_len(nrow(x)),function(i)
      fitted(npudens(bws=bw,tdat=x[-i,,drop=FALSE],edat=x[i,,drop=FALSE])),numeric(1))
    actual <- npudensbw.bandwidth(dat=x,bws=bw,bandwidth.compute=TRUE,
                                 eval.only=TRUE,nmulti=1L)$fval
    expect_equal(actual,sum(log(oracle)),tolerance=1e-11,info=type)
  }
})

test_that("beta density CVLS preserves its full-sample integral and deletes its cross term", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(.02,.06,.13,.23,.38,.51,.67,.8,.89,.98)
  n <- length(x); k <- 3L; q <- 81L
  grid <- seq(0,1,length.out=q)
  weights <- rep(1/(q-1),q); weights[c(1,q)] <- weights[c(1,q)]/2
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    bw <- npudensbw(dat=data.frame(x=x),bws=if(type=="fixed") .3 else k,
      bwtype=type,ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1,
      bandwidth.compute=FALSE,bwmethod="cv.ls")
    estimate <- function(target,omit=integer()) {
      donors <- setdiff(seq_len(n),omit)
      h <- switch(type,fixed=rep(.3,length(donors)),
        generalized_nn=rep(sort(abs(x[donors]-target))[k],length(donors)),
        adaptive_nn=vapply(donors,function(j)
          sort(abs(x[setdiff(donors,j)]-x[j]))[k],numeric(1)))
      mean(dbeta(x[donors],1+target/h^2,1+(1-target)/h^2))
    }
    full <- vapply(grid,estimate,numeric(1))
    cross <- vapply(seq_len(n),function(i) estimate(x[i],i),numeric(1))
    # The public fval reports the maximized negative CVLS loss.
    expected <- 2*mean(cross)-sum(weights*full^2)
    actual <- npudensbw.bandwidth(dat=data.frame(x=x),bws=bw,
      bandwidth.compute=TRUE,eval.only=TRUE,nmulti=1L)$fval
    expect_equal(actual,expected,tolerance=1e-11,info=type)
  }
})

test_that("conditional beta CVML selects matching X and Y folds", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920401)
  x <- data.frame(x=runif(21,.02,.98),c=factor(rep(letters[1:3],7)))
  y <- data.frame(y=runif(21,.02,.98))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(cx in c("beta","gaussian")) for(cy in c("beta","gaussian")) {
      if(cx=="gaussian" && cy=="gaussian") next
      for(reg in c("lc","ll","lp")) {
        args <- list(xdat=x,ydat=y,bws=c(rep(if(type=="fixed") .3 else 9,2),.2),
          bwtype=type,bandwidth.compute=FALSE,bwmethod="cv.ml",regtype=reg,
          cxkertype=cx,cykertype=cy)
        if(cx=="beta") args <- c(args,list(cxkerbound="fixed",cxkerlb=0,cxkerub=1))
        if(cy=="beta") args <- c(args,list(cykerbound="fixed",cykerlb=0,cykerub=1))
        if(reg=="lp") args$degree <- 2L
        bw <- do.call(npcdensbw,args)
        oracle <- vapply(seq_len(nrow(x)),function(i)
          fitted(npcdens(bws=bw,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
            exdat=x[i,,drop=FALSE],eydat=y[i,,drop=FALSE],se=FALSE)),numeric(1))
        # Signed LP densities retain the documented guarded-CVML objective.
        # Geometry changes neither its sign policy nor its underflow policy.
        log_terms <- log(abs(oracle))
        negative <- oracle < 0
        log_terms[negative] <- 2*log(.Machine$double.xmin)-log_terms[negative]
        log_terms[oracle==0] <- log(.Machine$double.xmin)
        expect_true(all(is.finite(oracle)),info=paste(type,cx,cy,reg))
        expect_equal(.npcdensbw_eval_only(x,y,bw)$objective,sum(log_terms),
                     tolerance=1e-9,info=paste(type,cx,cy,reg))
      }
    }
})
