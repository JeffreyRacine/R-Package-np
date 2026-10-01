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

test_that("beta GNN LC derivative hats retain implicit observation identity", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- c(.03, .10, .10, .24, .43, .67, .67, .82, .96)
  X <- data.frame(x = x)
  y <- sin(2*x) + x
  for (order in c(2L, 4L, 6L, 8L)) {
    bw <- npregbw(xdat = X, ydat = y, bws = 3,
      bandwidth.compute = FALSE, bwtype = "generalized_nn", regtype = "lc",
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1)
    for (external in c(FALSE, TRUE)) {
      expected <- t(vapply(seq_along(x), function(i) {
        # Holding this radius fixed defines the documented kernel derivative.
        donors <- if (external) x else x[-i]
        h <- sort(abs(donors-x[i]))[3L]
        w <- dw <- numeric(length(x))
        for (scale in seq_len(order/2L)) {
          tau <- 1/(scale*h*h)
          a <- 1+x[i]*tau; b <- 1+(1-x[i])*tau
          v <- dbeta(x, a, b)
          coefficient <- (-1)^(scale+1L)*choose(order/2L, scale)
          w <- w+coefficient*v
          dw <- dw+coefficient*v*tau*(log(x)-log1p(-x)-digamma(a)+digamma(b))
        }
        (dw*sum(w)-w*sum(dw))/sum(w)^2
      }, numeric(length(x))))
      args <- list(bws = bw, txdat = X, s = 1L)
      if (external) args$exdat <- X
      H <- do.call(npreghat, args)
      applied <- do.call(npreghat, c(args, list(output = "apply", y = y)))
      expect_equal(as.vector(H), as.vector(expected), tolerance = 8e-13)
      expect_equal(as.vector(applied), drop(expected %*% y), tolerance = 8e-13)
    }
  }
})
