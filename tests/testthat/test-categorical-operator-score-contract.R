r21_score_call <- function(ordered, kernel, lambda, operator, score = FALSE,
                           weights = NULL, response = NULL, loo = FALSE) {
  s <- c(0,1,3,6)
  ids <- rep(1:4,c(2,4,3,3))
  d <- data.frame(x = factor(s[ids], levels=s, ordered=ordered))
  e <- if(loo) d else data.frame(x=factor(s,levels=s,ordered=ordered))
  args <- list(txdat=d,bws=lambda,operator=operator,compute.score=score,
    return.kernel.weights=TRUE,return.derivative.kernel.weights=score,leave.one.out=loo)
  if(!loo) args$exdat <- e
  args[[if(ordered) "okertype" else "ukertype"]] <- kernel
  if(!is.null(weights)) args$weights <- weights
  if(!is.null(response)) args$tydat <- response
  do.call(npksum,args)
}

test_that("categorical scores differentiate each requested operator", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  for(ordered in c(FALSE,TRUE))
    for(kernel in if(ordered) c("liracine","nliracine","wangvanryzin","racineliyan")
        else c("aitchisonaitken","liracine"))
      for(op in c("normal","convolution","integral")) for(lambda in c(.17,.63)) {
        f <- r21_score_call(ordered,kernel,lambda,op,TRUE)
        a <- r21_score_call(ordered,kernel,lambda,op,FALSE)
        eps <- 1e-6
        plus <- r21_score_call(ordered,kernel,lambda+eps,op,FALSE)
        minus <- r21_score_call(ordered,kernel,lambda-eps,op,FALSE)
        expect_equal(f$ksum,a$ksum,tolerance=1e-13)
        expect_equal(f$kw,a$kw,tolerance=1e-13)
        expect_equal(as.numeric(f$p.ksum),as.numeric((plus$ksum-minus$ksum)/(2*eps)),
                     tolerance=2e-8)
        expect_equal(as.numeric(f$p.kw),as.numeric((plus$kw-minus$kw)/(2*eps)),
                     tolerance=2e-8)
      }
})

test_that("categorical operator scores retain endpoint limits", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  for(ordered in c(FALSE,TRUE))
    for(kernel in if(ordered) c("liracine","nliracine","wangvanryzin","racineliyan")
        else c("aitchisonaitken","liracine"))
      for(op in c("normal","convolution","integral"))
        for(lambda in c(0,if(kernel=="aitchisonaitken") .75 else 1)) {
          f <- r21_score_call(ordered,kernel,lambda,op,TRUE)
          step <- if(lambda==0) 1e-6 else -1e-6
          a <- r21_score_call(ordered,kernel,lambda,op,FALSE)
          b <- r21_score_call(ordered,kernel,lambda+step,op,FALSE)
          c <- r21_score_call(ordered,kernel,lambda+2*step,op,FALSE)
          expected <- (-3*a$ksum+4*b$ksum-c$ksum)/(2*step)
          exact <- ordered && kernel=="liracine" && op=="integral" && lambda==1
          if(exact) {
            # The old CDF closed form loses precision close to one. Its
            # derivative at one is the exact finite sum of integer distances.
            s <- c(0,1,3,6)
            ids <- rep(1:4,c(2,4,3,3))
            expected <- vapply(s,function(q)sum(abs(outer(s[ids],0:q,"-"))),0.0)
          }
          expect_true(all(is.finite(f$p.ksum)))
          expect_equal(as.numeric(f$p.ksum),as.numeric(expected),
                       tolerance=if(exact) 0 else 2e-7)
        }
})

test_that("operator scores preserve weighted multiresponse packing and LOO", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  W <- cbind(seq_len(12)/13,rev(seq_len(12))/17)
  Y <- cbind(sin(seq_len(12)),cos(seq_len(12)))
  for(op in c("normal","convolution","integral")) for(loo in c(FALSE,TRUE)) {
    f <- r21_score_call(TRUE,"racineliyan",.4,op,TRUE,W,Y,loo)
    a <- r21_score_call(TRUE,"racineliyan",.4,op,FALSE,W,Y,loo)
    plus <- r21_score_call(TRUE,"racineliyan",.400001,op,FALSE,W,Y,loo)
    minus <- r21_score_call(TRUE,"racineliyan",.399999,op,FALSE,W,Y,loo)
    expect_equal(f$ksum,a$ksum,tolerance=1e-13)
    expect_equal(as.numeric(f$p.ksum),as.numeric((plus$ksum-minus$ksum)/2e-6),
                 tolerance=2e-8)
    # LOO applies to sums; exported weights retain the raw diagonal.
    if(loo) expect_equal(as.numeric(f$p.kw),
      as.numeric((plus$kw-minus$kw)/2e-6),tolerance=2e-8)
  }
})

test_that("finite ordered support uses literal convolution and CDF derivatives", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  s <- c(0,1,3,6)
  ids <- rep(1:4,c(2,4,3,3))
  D <- abs(outer(s,s,"-"))
  lambda <- .4
  K <- lambda^D
  G <- ifelse(D==0,0,D*lambda^(D-1))
  expected <- (G %*% t(K)+K %*% t(G))[ids,,drop=FALSE]
  actual <- r21_score_call(TRUE,"liracine",lambda,"convolution",TRUE)
  expect_equal(as.numeric(actual$p.kw),as.numeric(expected),tolerance=2e-13)
  q <- data.frame(x=ordered(c(0,10000),levels=c(0,10000)))
  actual <- npksum(txdat=q,exdat=q,bws=1,okertype="liracine",
    operator="integral",compute.score=TRUE,return.kernel.weights=TRUE,
    return.derivative.kernel.weights=TRUE)
  expect_equal(as.numeric(actual$p.ksum),c(10000,100010000),tolerance=0)
  expect_equal(as.numeric(actual$p.kw),c(0,10000,50005000,50005000),tolerance=0)
})

test_that("cached LR convolution values share the retained-support owner", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  a <- c(0,1,3,6)
  b <- c(0,2,4,7)
  for(n in c(64L,128L,144L)) {
    d <- data.frame(a=ordered(rep(a,length.out=n),levels=a),
                    b=ordered(rep(b,each=4,length.out=n),levels=b))
    e <- d[1:16,,drop=FALSE]
    K <- .3^abs(outer(a,a,"-"))
    conv <- tcrossprod(K)
    expected <- conv[as.integer(d$a),as.integer(e$a)]*
      .4^abs(outer(b[as.integer(d$b)],b[as.integer(e$b)],"-"))
    f <- npksum(txdat=d,exdat=e,bws=c(.3,.4),okertype="liracine",
      operator=c("convolution","normal"),return.kernel.weights=TRUE)
    expect_equal(unname(f$kw),expected,tolerance=2e-13)
    expect_equal(as.numeric(f$ksum),colSums(expected),tolerance=2e-13)
  }
})

test_that("categorical operator scores retain mixed tree and NN weight ownership", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  s <- c(0,1,3,6)
  ids <- rep(1:4,12)
  d <- data.frame(x=seq(-3,3,length.out=48),o=ordered(s[ids],levels=s))
  e <- d[c(2,9,15,21,26,33,39,47),,drop=FALSE]
  D <- abs(outer(s,s,"-"))
  N <- .4^D
  G <- ifelse(D==0,0,D*.4^(D-1))
  den <- rowSums(N)
  G <- (G*den-N*rowSums(G))/den^2
  N <- N/den
  for(bwtype in c("fixed","generalized_nn","adaptive_nn"))
    for(tree in c(FALSE,TRUE)) for(op in c("convolution","integral")) {
      options(np.tree=tree)
      h <- if(bwtype=="fixed").2 else 5
      ck <- if(bwtype=="fixed")"epanechnikov" else "gaussian"
      cw <- npksum(txdat=d["x"],exdat=e["x"],bws=h,bwtype=bwtype,
                    ckertype=ck,return.kernel.weights=TRUE)$kw
      if(op=="convolution") {
        value <- tcrossprod(N)
        score <- G%*%t(N)+N%*%t(G)
      } else {
        value <- t(apply(N,1,cumsum))
        score <- t(apply(G,1,cumsum))
      }
      value <- cw*value[ids,as.integer(e$o)]
      score <- cw*score[ids,as.integer(e$o)]
      f <- npksum(txdat=d,exdat=e,bws=c(h,.4),bwtype=bwtype,
        ckertype=ck,okertype="racineliyan",operator=c("normal",op),
        compute.score=TRUE,return.kernel.weights=TRUE,
        return.derivative.kernel.weights=TRUE)
      expect_equal(as.numeric(f$kw),as.numeric(value),tolerance=2e-12)
      expect_equal(as.numeric(f$ksum),colSums(value),tolerance=2e-12)
      expect_equal(as.numeric(f$p.kw),as.numeric(score),tolerance=2e-12)
      expect_equal(as.numeric(f$p.ksum),colSums(score),tolerance=2e-12)
    }
})
