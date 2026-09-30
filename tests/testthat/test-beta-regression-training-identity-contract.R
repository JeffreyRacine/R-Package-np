test_that("beta GNN fits and hats exclude only training neighbour identity", {
  set.seed(4040)
  x <- runif(60L,.03,.97)
  x[8L] <- x[7L]
  training <- data.frame(x=x)
  response <- sin(4*x)+rnorm(length(x),sd=.15)
  for(degree in 0:2) {
    bw <- npregbw(xdat=training,ydat=response,bws=12,
      bandwidth.compute=FALSE,bwtype="generalized_nn",bwmethod="cv.aic",
      regtype=if(degree==0L) "lc" else if(degree==1L) "ll" else "lp",
      degree=degree,bernstein.basis=FALSE,
      ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1)
    weights <- workhorse_nn_training_beta_weights(training,12,"generalized_nn")
    design <- outer(x,0:degree,`^`)
    expected <- t(vapply(seq_along(x),function(i) {
      w <- weights[,i]/max(weights[,i])
      drop(design %*% solve(crossprod(design,design*w),design[i,]))*w
    },numeric(length(x))))
    fit <- npreg(bws=bw,txdat=training,tydat=response,se=FALSE)
    observed <- npreghat(bws=bw,txdat=training,output="matrix")
    expect_equal(unname(fitted(fit)),drop(expected %*% response),tolerance=2e-9)
    expect_equal(as.vector(observed),as.vector(expected),tolerance=2e-9)

    # Equal coordinates explicitly supplied as external queries do not become
    # training identities. A genuine fix must distinguish these two calls.
    external <- npreg(bws=bw,txdat=training,tydat=response,
                      exdat=training,se=FALSE)
    expect_gt(max(abs(fitted(fit)-fitted(external))),1e-5)
  }
})
